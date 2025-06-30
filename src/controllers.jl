import ControlSystemIdentification: PredictionStateSpace
import ControlSystemsBase: StateSpace, Discrete
import LinearAlgebra: pinv
import DataStructures: Queue, enqueue!, dequeue!, empty!, first
import LsqFit: curve_fit, coef

export lsim_step, model_step, state_prediction_matrices, Actuator, DelayedActuator,
    Controller, LinearController, PVLC, MPC, run_closed_loop_sim

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

function get_min_delay(act::Actuator)::Int
    min_delay = typemax(Int)
    has_delay = false
    for fieldname ∈ fieldnames(typeof(act))
        if fieldtype(typeof(act), fieldname) == DelayedActuator
            min_delay = min(min_delay, getfield(act, fieldname).delay)
            has_delay = true
        end
    end
    if has_delay
        return min_delay
    else
        return 0
    end
end

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
  - `act::Actuator`: Actuator model in the loop
  - `inp_cond::Union{Nothing, Function}=nothing`: Input conditioning on plant model
  - `inp_cond_kwargs::Dict{Symbol, Any}`: Keyword arguments for `inp_cond`
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
Additionally, it stores the `plant`, the number of steps of history `h` used for
state tracking and state prediction matrices `Y2x` and `U2x` from
[`state_prediction_matrices()`](@ref).
There is a convinience contructor:

    PVLC(
        ctrl_ss::StateSpace{TE},
        ctrl_x0::Vector{Float64},
        plant::Union{PredictionStateSpace, StateSpace},
        h::Int,
    ) where {TE <: Discrete}

This controller has a call signature:

    (pvlc::PVLC)(;
        ii::Int,
        target::Matrix{Float64},
        plant_inp::Matrix{Float64},
        plant_out::Matrix{Float64},
        act::Actuator,
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
    plant::Union{PredictionStateSpace, StateSpace}
    h::Int
    Y2x::Matrix{Float64}
    U2x::Matrix{Float64}
    function PVLC(
        ctrl_ss::StateSpace{TE},
        ctrl_x0::Vector{Float64},
        plant::Union{PredictionStateSpace, StateSpace},
        h::Int,
    ) where {TE <: Discrete}
        Y2x, U2x = state_prediction_matrices(plant, h)
        return new(ctrl_ss, ctrl_x0, plant, h, Y2x, U2x)
    end
end

function (pvlc::PVLC)(;
    ii::Int,
    target::Matrix{Float64},
    plant_inp::Matrix{Float64},
    plant_out::Matrix{Float64},
    act::Actuator,
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    kwargs...,
)::Vector{Float64}
    ctrl_out = zeros(Float64, size(pvlc.ctrl_ss.C, 1))
    future_inp = get_future_inp(act)
    lat = size(future_inp, 2)
    if ii - pvlc.h + 1 > 0 && ii + lat <= size(target, 2)
        # Gather history of inputs and outputs
        U = vec(plant_inp[:, (ii-pvlc.h+1):ii])
        Y = vec(plant_out[:, (ii-pvlc.h+1):ii])
        # Estimate current state vector of plant model
        est_x = pvlc.Y2x * Y .+ pvlc.U2x * U
        # Evolve into the future for the delay in actuator
        future_out = model_evolve(
            pvlc.plant, future_inp;
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

function expand(
    steps::StepRange{Int, Int},
    red_vec::Vector{Float64},
    ninps::Int,
    horizon::Int,
)::Vector{Vector{Float64}}
    red_mat = reshape(red_vec, ninps, length(steps))
    mat = zeros(Float64, ninps, horizon + 1)
    for i ∈ 1:ninps
        mat[i, :] =
            Interpolations.linear_interpolation(
                steps,
                red_mat[i, :];
                extrapolation_bc=Interpolations.Line(),
            ).(
                1:(horizon+1),
            )
    end
    return [mat[:, i] for i ∈ 1:(horizon+1)]
end

"""
    MPC

Implementation of simple leaqt square optimized Model Predictive Controller (MPC). It
stores the `plant`, the number of steps of history `h` used for state tracking and
state prediction matrices `Y2x` and `U2x` from [`state_prediction_matrices()`](@ref).
It stores a current deep copy of actuator instance  in `act` to try it during
optimization and it stores setting for least square optimization. It also stores a
`future_evolve` that is created based on `plant`, `act`, and optimization setting and
is the model function that is used for optimization later. This controller stores a
buffer for control outputs, `ctrl_out_buffer` so that it can be called less often and
it can reuse it's previous optimization results.

There is a convinience contructor:

    MPC(
        plant::Union{PredictionStateSpace, StateSpace},
        h::Int,
        act::Actuator,               # Actuator model without delay
        horizon::Int,                # Number of steps in future after latency
        nopt::Int,                   # Number of optimization points in horizon window
        opt_every::Int,              # Run cost optimization every opt_every steps
        ctrl_out_bounds::Tuple{Vector{Float64}, Vector{Float64}}=(
            Array{Float64}(undef, 0),
            Array{Float64}(undef, 0),
        ),
    )

This contructor takes in minimum required information to create a self-consistent MPC
instance. It sets the other dependent quantities in MPC such as `Y2x`, `U2x`,
`min_delay`, and create a `future_evolve` function and initializes the `ctrl_out_buffer`.
`act` needs to be a deepcopy of the actuator instance. `horizon` is the number of steps
after the `min_delay` among all acturators for which the optimization is carried out.
`nopt` is the number of optimization points taken along the `horizon` which are lineary
distributed. The gaps between this optimzation points are interpolated linearly in the
output. `opt_every` defines the frequency of optimization, i.e. at every `opt_every`
call of this controller, the optimization is carried out. This avoids unnecessary
over-calculations and thus results in a faster controller.

This controller has a call signature:

    (mpc::MPC)(;
        ii::Int,
        target::Matrix{Float64},
        plant_inp::Matrix{Float64},
        plant_out::Matrix{Float64},
        act::Actuator,
        inp_cond::Union{Nothing, Function}=nothing,
        inp_offset::Float64=0.0, inp_factor::Float64=1.0,
        out_offset::Float64=0.0, out_factor::Float64=1.0,
        inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
        inp_feedforward::Matrix{Float64}=zeros(
            Float64,
            (size(plant.B, 2), length(target)),
        ),
        kwargs...,
    )::Vector{Float64} where {T}

Tracks the state vector of the plant using `h` steps of history from `plant_inp` and
`plant_out` and uses it along with `act` to run an optimization to match `target` in
future after the minimum delay from all the actuators. This function uses
[`curve_fit()`](https://julianlsolvers.github.io/LsqFit.jl/latest/api/#LsqFit.curve_fit)
for the non-linear least squares fitting which uses Levenberg-Marquardt algorithm. This
function performs the optimization every `opt_every` call and uses the stored control
output in `ctrl_out_buffer` meanwhile.
"""
mutable struct MPC <: Controller
    plant::Union{PredictionStateSpace, StateSpace}
    h::Int
    Y2x::Matrix{Float64}
    U2x::Matrix{Float64}
    act::Actuator               # Actuator model without delay
    horizon::Int                # Number of steps in future after latency
    nopt::Int                   # Number of optimization points in horizon window
    opt_every::Int              # Run cost optimization every opt_every steps
    last_opt::Int               # Step index when optimization ran last time
    min_delay::Int              # Minimum delay in all actuators
    ctrl_out_bounds::Tuple{Vector{Float64}, Vector{Float64}}
    future_evolve::Function
    ctrl_out_buffer::Queue{Vector{Float64}}
    function MPC(
        plant::Union{PredictionStateSpace, StateSpace},
        h::Int,
        act::Actuator,               # Actuator model without delay
        horizon::Int,                # Number of steps in future after latency
        nopt::Int,                   # Number of optimization points in horizon window
        opt_every::Int,              # Run cost optimization every opt_every steps
        ctrl_out_bounds::Tuple{Vector{Float64}, Vector{Float64}}=(
            Array{Float64}(undef, 0),
            Array{Float64}(undef, 0),
        ),
    )
        Y2x, U2x = state_prediction_matrices(plant, h)
        min_delay = get_min_delay(act)
        function fe(
            red_steps::StepRange{Int, Int}, red_ctrl_out::Vector{Float64};
            mpc::MPC,
            inp_ff::Matrix{Float64}=zeros(
                Float64,
                size(mpc.plant.B, 2),
                min_delay + horizon + 1,
            ),
            x0::Vector{Float64}=zeros(Float64, size(mpc.plant.A, 2)),
            inp_cond::Union{Nothing, Function}=nothing,
            inp_offset::Float64=0.0, inp_factor::Float64=1.0,
            out_offset::Float64=0.0, out_factor::Float64=1.0,
            inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
        ) where {T}
            act = deepcopy(mpc.act)
            ninps = size(mpc.plant.B, 2)
            nouts = size(mpc.plant.C, 1)
            ctrl_out = expand(red_steps, red_ctrl_out, ninps, mpc.horizon)
            # Padding far in future extra zero inputs, these won't make any effect
            # because of the minium delay
            ctrl_out = [ctrl_out; repeat(zeros(Float64, ninps), mpc.min_delay)]
            plant_inp = zeros(Float64, ninps)
            plant_out = zeros(Float64, nouts, mpc.min_delay + mpc.horizon + 1)
            for i ∈ 1:(mpc.min_delay+mpc.horizon+1)
                plant_inp = act(ctrl_out[i] .+ inp_ff[:, i])
                plant_out[:, i], x0 = model_step(
                    mpc.plant, plant_inp;
                    x0,
                    inp_cond, inp_offset, inp_factor, out_offset, out_factor,
                    inp_cond_kwargs,
                )
            end
            return vec(plant_out)
        end
        last_opt = 0
        ctrl_out_buffer = Queue{Vector{Float64}}(horizon)
        return new(
            plant, h, Y2x, U2x, act, horizon, nopt, opt_every, last_opt,
            min_delay, ctrl_out_bounds, fe, ctrl_out_buffer,
        )
    end
end

function (mpc::MPC)(;
    ii::Int,
    target::Matrix{Float64},
    plant_inp::Matrix{Float64},
    plant_out::Matrix{Float64},
    act::Actuator,
    inp_cond::Union{Nothing, Function}=nothing,
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
    inp_feedforward::Matrix{Float64}=zeros(
        Float64,
        (size(plant.B, 2), length(target)),
    ),
    kwargs...,
)::Vector{Float64} where {T}
    fut_lookup = mpc.min_delay + mpc.horizon + 1

    enough_history = ii - mpc.h + 1 > 0
    enough_future = ii + fut_lookup <= size(target, 2)
    opt_step = ii > mpc.opt_every + mpc.last_opt
    buf_empty = isempty(mpc.ctrl_out_buffer)

    if buf_empty
        prev_ctrl_out = zeros(Float64, size(mpc.plant.B, 2))
    else
        prev_ctrl_out = first(mpc.ctrl_out_buffer)
    end

    if enough_history && enough_future && (opt_step || buf_empty)
        red_steps = 1:div(mpc.horizon, mpc.nopt-1):(mpc.horizon+1)
        lower = repeat(mpc.ctrl_out_bounds[1], length(red_steps))
        upper = repeat(mpc.ctrl_out_bounds[2], length(red_steps))
        ninps = size(mpc.plant.B, 2)

        # Gather history of inputs and outputs
        U = vec(plant_inp[:, (ii-mpc.h+1):ii])
        Y = vec(plant_out[:, (ii-mpc.h+1):ii])
        # Estimate current state vector of plant model
        x0 = mpc.Y2x * Y .+ mpc.U2x * U

        # Prepare arguments for cost optimization
        mpc.act = deepcopy(act)
        target_vec = vec(target[:, (ii+1):(ii+fut_lookup)])
        guess = repeat(prev_ctrl_out, length(red_steps))
        inp_ff = inp_feedforward[:, (ii+1):(ii+fut_lookup)]

        # Optimize
        fit_res = curve_fit(
            (x, y) -> mpc.future_evolve(
                x, y; mpc, inp_ff, x0,
                inp_cond, inp_offset, inp_factor, out_offset, out_factor,
                inp_cond_kwargs,
            ),
            red_steps, target_vec, guess;
            lower, upper,
        )

        # Empty buffer and fill with updated control outputs
        empty!(mpc.ctrl_out_buffer)
        for co ∈ expand(red_steps, coef(fit_res), ninps, mpc.horizon)
            enqueue!(mpc.ctrl_out_buffer, co)
        end
    end
    if opt_step
        mpc.last_opt = ii
    end
    if isempty(mpc.ctrl_out_buffer)
        return zeros(Float64, size(mpc.plant.B, 2))
    else
        return dequeue!(mpc.ctrl_out_buffer)
    end
end

"""
    run_closed_loop_sim(
        plant::Union{PredictionStateSpace{TE}, StateSpace{TE}},
        act::Actuator,
        ctrl::Controller,
        target::Matrix{Float64};
        inp_cond::Union{Nothing, Function}=nothing,
        inp_offset::Float64=0.0, inp_factor::Float64=1.0,
        out_offset::Float64=0.0, out_factor::Float64=1.0,
        inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
        inp_feedforward::Matrix{Float64}=zeros(
            Float64,
            (size(plant.B, 2), size(target, 2)),
        ),
        ctrl_start_ind::Int=1,
        noise_plant_inp::Matrix{Float64}=zeros(
            Float64,
            (size(plant.B, 2), size(target, 2)),
        ),
        noise_plant_out::Matrix{Float64}=zeros(
            Float64,
            (size(plant.C, 1), size(target, 2)),
        ),
        noise_ctrl_out::Matrix{Float64}=zeros(
            Float64,
            (size(plant.B, 2), size(target, 2)),
        ),
    ) where {T, TE <: Discrete}

Generic function to run closed loop simulations with provided `plant`, actuator
`act`, controller `ctrl`, and target waveform `target`. The length of simulation is
determined by provided `target`. Keyword arguments are possible for providing
adjustments to inputs and outputs of the plant model as explained in
[`model_evolve()`](@ref). Additionally, `ctrl_start_ind` can be provided to start
control loop at an arbitrary point in the loop. `noise_plant_inp`, `noise_plant_out`,
and `noise_ctrl_out` allow addition of predefined noise waveforms at the input of plant,
output of plant, and the output of controller respectively.
"""
function run_closed_loop_sim(
    plant::Union{PredictionStateSpace{TE}, StateSpace{TE}},
    act::Actuator,
    ctrl::Controller,
    target::Matrix{Float64};
    inp_cond::Union{Nothing, Function}=nothing,
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
    inp_feedforward::Matrix{Float64}=zeros(
        Float64,
        (size(plant.B, 2), size(target, 2)),
    ),
    ctrl_start_ind::Int=1,
    noise_plant_inp::Matrix{Float64}=zeros(
        Float64,
        (size(plant.B, 2), size(target, 2)),
    ),
    noise_plant_out::Matrix{Float64}=zeros(
        Float64,
        (size(plant.C, 1), size(target, 2)),
    ),
    noise_ctrl_out::Matrix{Float64}=zeros(
        Float64,
        (size(plant.B, 2), size(target, 2)),
    ),
) where {T, TE <: Discrete}
    # Initialization
    ctrl_out = zeros(Float64, size(inp_feedforward))
    plant_inp = zeros(Float64, (size(plant.B, 2), length(target)))
    plant_out = zeros(Float64, (size(plant.C, 1), length(target)))
    x0 =
        inv(Matrix{Float64}(I, size(plant.A)...) - plant.A) *
        plant.B * inp_feedforward[:, 1]

    # Closed loop simulation
    for ii ∈ eachindex(target)
        # Actuation
        plant_inp[:, ii] =
            act(ctrl_out[:, ii] .+ inp_feedforward[:, ii]) .+ noise_plant_inp[:, ii]

        # Plant response
        plant_out[:, ii], x0 = model_step(
            plant, plant_inp[:, ii];
            x0,
            inp_cond, inp_offset, inp_factor, out_offset, out_factor,
            inp_cond_kwargs,
        )
        plant_out[:, ii] .+= noise_plant_out[:, ii]
        # future_inp = get_future_inp(act)

        # Controller
        if ii >= ctrl_start_ind && ii < length(target)
            ctrl_out[:, ii+1] =
                ctrl(;
                    ii, target, plant_inp, plant_out, act, inp_feedforward,
                    inp_cond, inp_cond_kwargs,
                    inp_offset, inp_factor, out_offset, out_factor,
                ) .+ noise_ctrl_out[:, ii+1]
        end
    end
    return Dict(
        :plant => plant,
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
