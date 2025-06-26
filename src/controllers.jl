import ControlSystemIdentification: PredictionStateSpace
import ControlSystemsBase: StateSpace, Discrete
import LinearAlgebra: pinv
using DataStructures: Queue, enqueue!, dequeue!
# import ControlSystemsBase: lsim, ss, pid, d2c_exact, d2c, StateSpace

export lsim_step, model_step, state_prediction_matrices, Actuator, DelayedActuator,
    Controller, LinearController, PVLC, run_closed_loop_sim

"""
    lsim_step(
        sys::Union{PredictionStateSpace, StateSpace}, u::Vector{Float64};
        x0::Vector{Float64}=zeros(size(sys.A, 1)),
    )::Tuple{Vector{Float64}, Vector{Float64}}

Single step version of [ControlSystemsBase.lsim](https://juliacontrol.github.io/ControlSystems.jl/dev/lib/timefreqresponse/#ControlSystemsBase.lsim-Tuple%7BAbstractStateSpace,%20AbstractVecOrMat,%20AbstractVector%7D)
"""
function lsim_step(
    sys::Union{PredictionStateSpace, StateSpace}, u::Vector{Float64};
    x0::Vector{Float64}=zeros(size(sys.A, 1)),
)::Tuple{Vector{Float64}, Vector{Float64}}
    x = sys.A * x0 + sys.B * u
    y = sys.C * x0 + sys.D * u
    return y, x
end

"""
    model_step(
        sys::Union{PredictionStateSpace, StateSpace},
        inp::Vector{Float64};
        x0::Vector{Float64}=zeros(size(sys.A, 1)),
        inp_cond::Union{Nothing, Function}=nothing,
        inp_offset::Float64=0.0, inp_factor::Float64=1.0,
        out_offset::Float64=0.0, out_factor::Float64=1.0,
        inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
    )::Tuple{Vector{Float64}, Vector{Float64}} where {T}

Single step version of [`model_evolve()`](@ref).
"""
function model_step(
    sys::Union{PredictionStateSpace, StateSpace},
    inp::Vector{Float64};
    x0::Vector{Float64}=zeros(size(sys.A, 1)),
    inp_cond::Union{Nothing, Function}=nothing,
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
)::Tuple{Vector{Float64}, Vector{Float64}} where {T}
    inp_so = offset_scale(inp; offset=inp_offset, factor=inp_factor)
    if !isnothing(inp_cond)
        inp_so = inp_cond(inp_so; inp_cond_kwargs...)
    end
    modeled_out_so, x0 = lsim_step(sys, inp_so; x0)
    modeled_out = unscale_unoffset(modeled_out_so; offset=out_offset, factor=out_factor)
    return modeled_out, x0
end

function _calculate_LMNO(
    sys::Union{PredictionStateSpace, StateSpace},
    h::Int=size(sys.A, 1),
)
    A = sys.A
    B = sys.B
    C = sys.C
    D = sys.D
    order = size(A, 1)
    ninps = size(B, 2)
    nouts = size(C, 1)

    L = zeros((h * nouts, order))
    M = zeros((h * nouts, h * ninps))
    O = zeros((order, h * ninps))
    Atoim1 = Matrix{Float64}(I, (order, order))

    for j ∈ 1:h
        # println("Filling M with D: ", (nouts*(j-1)+1):(nouts*j), ", ", (nouts*(j-1)+1):(nouts*j))
        M[(nouts*(j-1)+1):(nouts*j), (ninps*(j-1)+1):(ninps*j)] = D
    end
    for i ∈ 1:h
        # println("i = ", i)
        if i > 1
            Atoim1 = Atoim1 * A
        end
        CAtoim1 = C * Atoim1
        # println("Filling L[$((nouts*(i-1)+1):(nouts*i)), :] = CA^$(i-1)")
        L[(nouts*(i-1)+1):(nouts*i), :] = CAtoim1

        for j ∈ (i+1):h
            # println("Filling M[$((nouts*(j-1)+1):(nouts*j)), $((ninps*(j-1-i)+1):(ninps*(j-i)))] = CA^$(i-1)B")
            M[(nouts*(j-1)+1):(nouts*j), (ninps*(j-1-i)+1):(ninps*(j-i))] = CAtoim1 * B
        end
        # println("Filling O[:, $((ninps*(h-i)+1):(ninps*(h-i+1)))] = A^$(i-1)B")
        O[:, (ninps*(h-i)+1):(ninps*(h-i+1))] = Atoim1 * B
    end
    N = Atoim1 * A
    return L, M, N, O
end

"""
    state_prediction_matrices(
        sys::Union{PredictionStateSpace, StateSpace},
        h::Int=size(sys.A, 1),
    )

Calculate state prediction matrices for a linear system `sys` given `h` steps of input
and output history. This function returns two matrices, `Y2x` and `U2x` that can be
used to calculate least square fitted estimate of the current state vector of the
system.
`Y2x` has size: (State size of `sys`) x (No. of outputs times h)
`U2x` has size: (State size of `sys`) x (No. of inputs times h)
The estimated state vector is obtained by:
Y2x * Y + U2x * U
where Y is all the `h` outputs of system stacked into a single vector.
and U is all the `h` inputs to the system stacked into a single vector.
"""
function state_prediction_matrices(
    sys::Union{PredictionStateSpace, StateSpace},
    h::Int=size(sys.A, 1),
)
    L, M, N, O = _calculate_LMNO(sys, h)
    Linv = pinv(L)
    Y2x = N * Linv
    U2x = (O - Y2x * M)
    return Y2x, U2x
end

"""
    Actuator

Abstract parent type for all actuators. Whenever a user defined actuator is created,
it must be subtype of `Actuator`.

To create a new actuator with custom function, it must be defined as a mutable stucture
which is daughter of Actuator that contains all settings and state information for
the actuator and the instance itself should be callable to take as input a
`Vector{Float64}` and should return an output of `Vector{Float64}`.

```julia
mutable struct CustomActuator <: Actuator
    settings
    state
    # ... Anything more
end

function (ca::CustomActuator)(inp::Vector{Float64})::Vector{Float64}
    # perform the actuation calcualtions with inp
    # update ca.state if required
    return output
end
```

**NOTE:** If you need to add a delay in the actuator, add a [`DelayedActuator`](@ref)
instance in the attributes of your `CustomActuator` and just call the DelayedActuator
inside your function call.
"""
abstract type Actuator end

"""
    DelayedActuator{U}

Implementation of delayed actuation. It stores `delay::Int` for number of time steps of
delay and `buffer::Queue{U}` which stores the delayed actuations in a queue.
Constructor:

    DelayedActuator(
        delay::Int;
        default_output::T=[0.0],
    ) where {T}

Creates a DelayedActuator{T} instance with `delay` and initializes the `buffer`
pre-filled upto brim with the default_output.
"""
mutable struct DelayedActuator{U} <: Actuator
    delay::Int
    buffer::Queue{U}
    function DelayedActuator(
        delay::Int;
        default_output::T=[0.0],
    ) where {T}
        buffer = Queue{T}(delay)
        for i ∈ 1:delay
            enqueue!(buffer, default_output)
        end
        return new{T}(delay, buffer)
    end
end

function (da::DelayedActuator)(inp::Union{Float64, Vector{Float64}})::Vector{Float64}
    enqueue!(da.buffer, inp)
    return dequeue!(da.buffer)
end

function get_future_inp(act::Actuator)
    for fieldname ∈ fieldnames(typeof(act))
        if fieldtype(typeof(act), fieldname) == DelayedActuator
            return get_future_inp(getfield(act, fieldname))
        end
    end
    return Matrix{Float64}(undef, 0, 0)
end

get_future_inp(act::DelayedActuator) = stack(act.buffer)

"""
    Controller

Abstract parent type for all controllers. Whenever a user defined controller is created
it must be a subtype of `Controller`.

To create a new controller algorithm, it should be defined as a mutable structure that
is a daughter of `Controller` and should contain all settings and state information to
be stored. It must be a callable structure that can use any of the following keyword
arguments:

  - `ii::Int`: Iteration index
  - `target::Matrix{Float64}`: Target waveform (No. of signals x Time steps)
  - `plant_inp::Matrix{Float64}`: Inputs to plant (No. of inputs x Time steps)
  - `plant_out::Matrix{Float64}`: Outputs from plant (No. of outputs x Time steps)
  - `future_inp::Matrix{Float64}`: Upcoming future known inputs to plant (for delayed actuators) (No. of outputs x Time steps)
  - `inp_offset::Float64`: Input offset for plant input
  - `inp_factor::Float64`: Input factor for plant input
  - `out_offset::Float64`: Output offset for plant output
  - `out_factor::Float64`: Output factor for plant output
  - `kwargs..`: Required to ignore unused keyword arguments
"""
abstract type Controller end

"""
    LinearController

Implementation for using any linear controller. It stores
`ctrl_ss::StateSpace{TE} where {TE <: Discrete}` for storing any linear controller as a
discrete state space model using [ControlSystemsBase.ss](https://juliacontrol.github.io/ControlSystems.jl/dev/man/creating_systems/#State-Space-Systems).
It also stores the state vector for the state space model as `ctrl_x0::Vector{Float64}`.
It's call signature is:

    (lc::LinearController)(;
        ii::Int,
        target::Matrix{Float64},
        plant_out::Matrix{Float64},
        kwargs...,
    )::Vector{Float64}

Calcualtes error as `target[:, ii] .- plant_out[:, ii]` and runs it through
[`lsim_step()`](@ref).
"""
mutable struct LinearController <: Controller
    ctrl_ss::StateSpace{TE} where {TE <: Discrete}
    ctrl_x0::Vector{Float64}
end

function (lc::LinearController)(;
    ii::Int,
    target::Matrix{Float64},
    plant_out::Matrix{Float64},
    kwargs...,
)::Vector{Float64}
    err = target[:, ii] .- plant_out[:, ii]
    ctrl_out, lc.ctrl_x0 = lsim_step(lc.ctrl_ss, err; x0=lc.ctrl_x0)
    return ctrl_out
end

"""
    PVLC

Implementation of Predicted Variable Linear Controller (PVLC). It stores
`ctrl_ss::StateSpace{TE} where {TE <: Discrete}` for storing any linear controller as a
discrete state space model using [ControlSystemsBase.ss](https://juliacontrol.github.io/ControlSystems.jl/dev/man/creating_systems/#State-Space-Systems).
It also stores the state vector for the state space model as `ctrl_x0::Vector{Float64}`.
Additionally, it stores the `plant_model`, the number of steps of history `h` used for
state tracking and state prediction matrices `Y2x` and `U2x` from
[`state_prediction_matrices()`](@ref).
There is a convinience contructor:

    PVLC(
        ctrl_ss::StateSpace{TE},
        ctrl_x0::Vector{Float64},
        plant_model::Union{PredictionStateSpace, StateSpace},
        h::Int,
    ) where {TE <: Discrete}

This controller has a call signature:

    (pvlc::PVLC)(;
        ii::Int,
        target::Matrix{Float64},
        plant_inp::Matrix{Float64},
        plant_out::Matrix{Float64},
        future_inp::Matrix{Float64},
        inp_offset::Float64=0.0, inp_factor::Float64=1.0,
        out_offset::Float64=0.0, out_factor::Float64=1.0,
        kwargs...,
    )::Vector{Float64}

Tracks the state vector of the plant using `h` steps of history from `plant_inp` and
`plant_out` and uses it to calculate future output of the plant. It compares it with a
future `target` value and applies the linear controller `ctrl_ss` there.
"""
mutable struct PVLC <: Controller
    ctrl_ss::StateSpace{TE} where {TE <: Discrete}
    ctrl_x0::Vector{Float64}
    plant_model::Union{PredictionStateSpace, StateSpace}
    h::Int
    Y2x::Matrix{Float64}
    U2x::Matrix{Float64}
    function PVLC(
        ctrl_ss::StateSpace{TE},
        ctrl_x0::Vector{Float64},
        plant_model::Union{PredictionStateSpace, StateSpace},
        h::Int,
    ) where {TE <: Discrete}
        Y2x, U2x = state_prediction_matrices(plant_model, h)
        return new(ctrl_ss, ctrl_x0, plant_model, h, Y2x, U2x)
    end
end

function (pvlc::PVLC)(;
    ii::Int,
    target::Matrix{Float64},
    plant_inp::Matrix{Float64},
    plant_out::Matrix{Float64},
    future_inp::Matrix{Float64},
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    kwargs...,
)::Vector{Float64}
    ctrl_out = zeros(Float64, size(pvlc.ctrl_ss.C, 1))
    lat = size(future_inp, 2)
    if ii - pvlc.h + 1 > 0 && ii + lat <= size(target, 2)
        # Gather history of inputs and outputs
        U = vec(plant_inp[:, (ii-pvlc.h+1):ii])
        Y = vec(plant_out[:, (ii-pvlc.h+1):ii])
        # Estimate current state vector of plant model
        est_x = pvlc.Y2x * Y .+ pvlc.U2x * U
        # Evolve into the future for the delay in actuator
        future_out = model_evolve(
            pvlc.plant_model, future_inp;
            x0=est_x,
            inp_offset, inp_factor, out_offset, out_factor,
        )
        # Estimate the error signal in future
        err = target[:, ii+lat] - future_out[:, end]
        # Apply the linear controller in future
        ctrl_out, pvlc.ctrl_x0 = lsim_step(pvlc.ctrl_ss, err; x0=pvlc.ctrl_x0)
    end
    return ctrl_out
end

"""
    run_closed_loop_sim(
        plant_model::Union{PredictionStateSpace{TE}, StateSpace{TE}},
        act::Actuator,
        ctrl::Controller,
        target::Matrix{Float64};
        inp_cond::Union{Nothing, Function}=nothing,
        inp_offset::Float64=0.0, inp_factor::Float64=1.0,
        out_offset::Float64=0.0, out_factor::Float64=1.0,
        inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
        inp_feedforward::Matrix{Float64}=zeros(
            Float64,
            (size(plant_model.B, 2), length(target)),
        ),
        ctrl_start_ind::Int=1,
    ) where {T, TE <: Discrete}

Generic function to run closed loop simulations with provided `plant_model`, actuator
`act`, controller `ctrl`, and target waveform `target`. The length of simulation is
determined by provided `target`. Keyword arguments are possible for providing
adjustments to inputs and outputs of the plant model as explained in
[`model_evolve()`](@ref). Additionally, `ctrl_start_ind` can be provided to start
control loop at an arbitrary point in the loop.
"""
function run_closed_loop_sim(
    plant_model::Union{PredictionStateSpace{TE}, StateSpace{TE}},
    act::Actuator,
    ctrl::Controller,
    target::Matrix{Float64};
    inp_cond::Union{Nothing, Function}=nothing,
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
    inp_feedforward::Matrix{Float64}=zeros(
        Float64,
        (size(plant_model.B, 2), size(target, 2)),
    ),
    ctrl_start_ind::Int=1,
    noise_plant_inp::Matrix{Float64}=zeros(
        Float64,
        (size(plant_model.B, 2), size(target, 2)),
    ),
    noise_plant_out::Matrix{Float64}=zeros(
        Float64,
        (size(plant_model.C, 1), size(target, 2)),
    ),
    noise_ctrl_out::Matrix{Float64}=zeros(
        Float64,
        (size(plant_model.B, 2), size(target, 2)),
    ),
) where {T, TE <: Discrete}
    ctrl_out = zeros(Float64, size(inp_feedforward))
    plant_inp = zeros(Float64, (size(plant_model.B, 2), length(target)))
    plant_out = zeros(Float64, (size(plant_model.C, 1), length(target)))
    x0 =
        inv(Matrix{Float64}(I, size(plant_model.A)...) - plant_model.A) *
        plant_model.B * inp_feedforward[:, 1]
    for ii ∈ eachindex(target)
        plant_inp[:, ii] =
            act(ctrl_out[:, ii] .+ inp_feedforward[:, ii]) .+ noise_plant_inp[:, ii]
        plant_out[:, ii], x0 = model_step(
            plant_model, plant_inp[:, ii];
            x0,
            inp_cond, inp_offset, inp_factor, out_offset, out_factor,
            inp_cond_kwargs,
        )
        plant_out[:, ii] .+= noise_plant_out[:, ii]
        future_inp = get_future_inp(act)
        if ii >= ctrl_start_ind && ii < length(target)
            ctrl_out[:, ii+1] =
                ctrl(;
                    ii, target, plant_inp, plant_out, future_inp,
                    inp_offset, inp_factor, out_offset, out_factor,
                ) .+ noise_ctrl_out[:, ii+1]
        end
    end
    return Dict(
        :plant_model => plant_model,
        :act => act,
        :ctrl => ctrl,
        :target => target,
        :inp_cond => inp_cond,
        :inp_offset => inp_offset,
        :inp_factor => inp_factor,
        :out_offset => out_offset,
        :out_factor => out_factor,
        :inp_cond_kwargs => inp_cond_kwargs,
        :inp_feedforward => inp_feedforward,
        :ctrl_start_ind => ctrl_start_ind,
        :plant_inp => plant_inp,
        :plant_out => plant_out,
        :ctrl_out => ctrl_out,
        :noise_plant_inp => noise_plant_inp,
        :noise_plant_out => noise_plant_out,
        :noise_ctrl_out => noise_ctrl_out,
    )
end
