import ControlSystemIdentification: PredictionStateSpace
import ControlSystemsBase: StateSpace, Discrete
import LinearAlgebra: pinv
import DataStructures: Queue, enqueue!, dequeue!, empty!
import LsqFit: curve_fit, coef

export Controller,
    LinearController, PVLC, state_prediction_matrices, MPC, run_closed_loop_sim

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
  - `inp_feedforward::Matrix{Float64}`: Feedforward input to the plant (No. of inputs x Time steps)
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
    PVLC

    # Constructor
    PVLC(
        ctrl_ss::StateSpace{TE},
        ctrl_x0::Vector{Float64},
        plant::LinearPlant,
        h::Int,
    ) where {TE <: Discrete}

Implementation of Predicted Variable Linear Controller (PVLC). It stores
`ctrl_ss::StateSpace{TE} where {TE <: Discrete}` for storing any linear controller as a
discrete state space model using [ControlSystemsBase.ss](https://juliacontrol.github.io/ControlSystems.jl/dev/man/creating_systems/#State-Space-Systems).
It also stores the state vector for the state space model as `ctrl_x0::Vector{Float64}`.
Additionally, it stores the `plant`, the number of steps of history `h` used for
state tracking and state prediction matrices `Y2x` and `U2x` from
[`state_prediction_matrices()`](@ref).

This controller has a call signature:

    (pvlc::PVLC)(;
        ii::Int,
        target::Matrix{Float64},
        plant_inp::Matrix{Float64},
        plant_out::Matrix{Float64},
        act::Actuator,
        kwargs...,
    )::Vector{Float64}

Tracks the state vector of the plant using `h` steps of history from `plant_inp` and
`plant_out` and uses it to calculate future output of the plant. It compares it with a
future `target` value and applies the linear controller `ctrl_ss` there.
"""
mutable struct PVLC <: Controller
    ctrl_ss::StateSpace{TE} where {TE <: Discrete}
    ctrl_x0::Vector{Float64}
    plant::LinearPlant
    h::Int
    Y2x::Matrix{Float64}
    U2x::Matrix{Float64}
    function PVLC(
        ctrl_ss::StateSpace{TE},
        ctrl_x0::Vector{Float64},
        plant::LinearPlant,
        h::Int,
    ) where {TE <: Discrete}
        # Take the plant by value
        plant = deepcopy(plant)
        Y2x, U2x = state_prediction_matrices(plant.sys, h)
        return new(ctrl_ss, ctrl_x0, plant, h, Y2x, U2x)
    end
end

function (pvlc::PVLC)(;
    ii::Int,
    target::Matrix{Float64},
    plant_inp::Matrix{Float64},
    plant_out::Matrix{Float64},
    act::Actuator,
    kwargs...,
)::Vector{Float64}
    ctrl_out = zeros(Float64, size(pvlc.ctrl_ss.C, 1))
    future_inp = get_future_inp(act)
    lat = size(future_inp, 2)
    if ii - pvlc.h + 1 > 0 && ii + lat <= size(target, 2)
        # Gather history of inputs and outputs
        U = offset_scale(
            vec(plant_inp[:, (ii-pvlc.h+1):ii]);
            offset=pvlc.plant.inp_offset, factor=pvlc.plant.inp_factor,
        )
        Y = offset_scale(
            vec(plant_out[:, (ii-pvlc.h+1):ii]);
            offset=pvlc.plant.out_offset, factor=pvlc.plant.out_factor,
        )
        # Estimate current state vector of plant model
        pvlc.plant.x = pvlc.Y2x * Y .+ pvlc.U2x * U
        # Evolve into the future for the delay in actuator
        future_out = model_evolve(pvlc.plant, future_inp)
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

    # Constructor
    MPC(
        plant::Plant,
        h::Int,
        act::Actuator,               # Actuator model without delay
        horizon::Int,                # Number of steps in future after latency
        nopt::Int,                   # Number of optimization points in horizon window
        opt_every::Int;              # Run cost optimization every opt_every steps
        ctrl_out_bounds::Tuple{Vector{Float64}, Vector{Float64}}=(
            Array{Float64}(undef, 0),
            Array{Float64}(undef, 0),
        ),
        guess::Union{Symbol, Vector{Float64}}=:zeros,
        curve_fit_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
    )

Implementation of simple leaqt square optimized Model Predictive Controller (MPC). It
stores the `plant`, the number of steps of history `h` used for state tracking and
state prediction matrices `Y2x` and `U2x` from [`state_prediction_matrices()`](@ref).
It stores a current deep copy of actuator instance  in `act` to try it during
optimization and it stores setting for least square optimization. It also stores a
`future_evolve` that is created based on `plant`, `act`, and optimization setting and
is the model function that is used for optimization later. This controller stores a
buffer for control outputs, `ctrl_out_buffer` so that it can be called less often and
it can reuse it's previous optimization results.

This contructor takes in minimum required information to create a self-consistent MPC
instance. It sets the other dependent quantities in MPC such as `Y2x`, `U2x`,
`min_delay`, and create a `future_evolve` function and initializes the `ctrl_out_buffer`.
`act` needs to be a deepcopy of the actuator instance. `horizon` is the number of steps
after the `min_delay` among all acturators for which the optimization is carried out.
`nopt` is the number of optimization points taken along the `horizon` which are lineary
distributed. The gaps between this optimzation points are interpolated linearly in the
output. `opt_every` defines the frequency of optimization, i.e. at every `opt_every`
call of this controller, the optimization is carried out. This avoids unnecessary
over-calculations and thus results in a faster controller. `ctrl_out_bounds` is a tuple
of lower bounds and upper bounds for the control output. `guess` is used to create the
initial guess during least square optimization. If `:last`, it would use the last
controller output as initial setting. If it is a `Vector`, each initialization starts
with this Vector. In all other cases, initial guess is zeros. `curve_fit_kwargs` can be
used to provide keyword arguments that go to `curve_fit`.

This controller has a call signature:

    function (mpc::MPC)(;
        ii::Int,
        target::Matrix{Float64},
        plant_inp::Matrix{Float64},
        plant_out::Matrix{Float64},
        act::Actuator,
        inp_feedforward::Matrix{Float64}=zeros(
            Float64,
            (size(get_sys(plant).B, 2), length(target)),
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
mutable struct MPC{U} <: Controller
    plant::LinearPlant
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
    guess::Union{Symbol, Vector{Float64}}
    curve_fit_kwargs::Dict{Symbol, U}
    function MPC(
        plant::LinearPlant,
        h::Int,
        act::Actuator,               # Actuator model without delay
        horizon::Int,                # Number of steps in future after latency
        nopt::Int,                   # Number of optimization points in horizon window
        opt_every::Int;              # Run cost optimization every opt_every steps
        ctrl_out_bounds::Tuple{Vector{Float64}, Vector{Float64}}=(
            Array{Float64}(undef, 0),
            Array{Float64}(undef, 0),
        ),
        guess::Union{Symbol, Vector{Float64}}=:zeros,
        curve_fit_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
    ) where {T}
        # Take the plant and actuator by value
        plant = deepcopy(plant)
        act = deepcopy(act)
        Y2x, U2x = state_prediction_matrices(plant.sys, h)
        min_delay = get_min_delay(act)
        function fe(
            red_steps::StepRange{Int, Int}, red_ctrl_out::Vector{Float64};
            mpc::MPC,
            inp_ff::Matrix{Float64}=zeros(
                Float64,
                size(mpc.plant.sys.B, 2),
                min_delay + horizon + 1,
            ),
        )
            act = deepcopy(mpc.act)
            ninps = size(mpc.plant.sys.B, 2)
            nouts = size(mpc.plant.sys, 1)
            ctrl_out = expand(red_steps, red_ctrl_out, ninps, mpc.horizon)
            # Padding far in future extra zero inputs, these won't make any effect
            # because of the minium delay
            ctrl_out = [ctrl_out; repeat(zeros(Float64, ninps), mpc.min_delay)]
            plant_inp = zeros(Float64, ninps)
            plant_out = zeros(Float64, nouts, mpc.min_delay + mpc.horizon + 1)
            for i ∈ 1:(mpc.min_delay+mpc.horizon+1)
                plant_inp = act(ctrl_out[i] .+ inp_ff[:, i])
                plant_out[:, i] = mpc.plant(plant_inp)
            end
            return vec(plant_out[:, (1+mpc.min_delay):end])
        end
        last_opt = -opt_every
        ctrl_out_buffer = Queue{Vector{Float64}}(horizon)
        return new{T}(
            plant, h, Y2x, U2x, act, horizon, nopt, opt_every, last_opt,
            min_delay, ctrl_out_bounds, fe, ctrl_out_buffer, guess, curve_fit_kwargs,
        )
    end
end

function (mpc::MPC)(;
    ii::Int,
    target::Matrix{Float64},
    plant_inp::Matrix{Float64},
    plant_out::Matrix{Float64},
    act::Actuator,
    inp_feedforward::Matrix{Float64}=zeros(
        Float64,
        (size(get_sys(plant).B, 2), length(target)),
    ),
    kwargs...,
)::Vector{Float64}
    fut_lookup = mpc.min_delay + mpc.horizon + 1

    enough_history = ii - mpc.h + 1 > 0
    enough_future = ii + fut_lookup <= size(target, 2)
    opt_step = ii > mpc.opt_every + mpc.last_opt
    buf_empty = isempty(mpc.ctrl_out_buffer)

    if enough_history && enough_future && (opt_step || buf_empty)
        red_steps = 1:div(mpc.horizon, mpc.nopt-1):(mpc.horizon+1)
        lower = repeat(mpc.ctrl_out_bounds[1], length(red_steps))
        upper = repeat(mpc.ctrl_out_bounds[2], length(red_steps))
        ninps = size(mpc.plant.sys.B, 2)

        # Gather history of inputs and outputs
        U = offset_scale(
            vec(plant_inp[:, (ii-mpc.h+1):ii]);
            offset=mpc.plant.inp_offset, factor=mpc.plant.inp_factor,
        )
        Y = offset_scale(
            vec(plant_out[:, (ii-mpc.h+1):ii]);
            offset=mpc.plant.out_offset, factor=mpc.plant.out_factor,
        )
        # Estimate current state vector of plant model
        mpc.plant.x = mpc.Y2x * Y .+ mpc.U2x * U

        # Prepare arguments for cost optimization
        mpc.act = deepcopy(act)
        target_vec = vec(target[:, (ii+1+mpc.min_delay):(ii+fut_lookup)])

        if mpc.guess == :last && !isempty(mpc.ctrl_out_buffer)
            guess = repeat(first(mpc.ctrl_out_buffer), length(red_steps))
        elseif isa(mpc.guess, Vector)
            guess = repeat(mpc.guess, length(red_steps))
        else
            guess = repeat(zeros(Float64, ninps), length(red_steps))
        end
        inp_ff = inp_feedforward[:, (ii+1):(ii+fut_lookup)]

        # Optimize
        fit_res = curve_fit(
            (x, y) -> mpc.future_evolve(x, y; mpc, inp_ff),
            red_steps, target_vec, guess;
            lower, upper,
            mpc.curve_fit_kwargs...,
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
        return zeros(Float64, size(mpc.plant.sys.B, 2))
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
        inp_feedforward::Matrix{Float64}=zeros(
            Float64,
            (size(get_sys(plant).B, 2), size(target, 2)),
        ),
        ctrl_start_ind::Int=1,
        noise_plant_inp::Matrix{Float64}=zeros(
            Float64,
            (size(get_sys(plant).B, 2), size(target, 2)),
        ),
        noise_plant_out::Matrix{Float64}=zeros(
            Float64,
            (size(get_sys(plant).C, 1), size(target, 2)),
        ),
        noise_ctrl_out::Matrix{Float64}=zeros(
            Float64,
            (size(get_sys(plant).B, 2), size(target, 2)),
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
    plant::Plant,
    act::Actuator,
    ctrl::Controller,
    target::Matrix{Float64};
    inp_feedforward::Matrix{Float64}=zeros(
        Float64,
        (size(get_sys(plant).B, 2), size(target, 2)),
    ),
    ctrl_start_ind::Int=1,
    noise_plant_inp::Matrix{Float64}=zeros(
        Float64,
        (size(get_sys(plant).B, 2), size(target, 2)),
    ),
    noise_plant_out::Matrix{Float64}=zeros(
        Float64,
        (size(get_sys(plant).C, 1), size(target, 2)),
    ),
    noise_ctrl_out::Matrix{Float64}=zeros(
        Float64,
        (size(get_sys(plant).B, 2), size(target, 2)),
    ),
)
    # Take the plant, actuator, and controller by value
    plant = deepcopy(plant)
    act = deepcopy(act)
    ctrl = deepcopy(ctrl)

    # Initialization
    ctrl_out = zeros(Float64, size(inp_feedforward))
    plant_sys = get_sys(plant)
    plant_inp = zeros(Float64, (size(plant_sys.B, 2), length(target)))
    plant_out = zeros(Float64, (size(plant_sys.C, 1), length(target)))

    # Closed loop simulation
    for ii ∈ axes(target, 2)
        # Actuation
        plant_inp[:, ii] =
            act(ctrl_out[:, ii] .+ inp_feedforward[:, ii]) .+ noise_plant_inp[:, ii]

        # Plant response
        plant_out[:, ii] = plant(plant_inp[:, ii])
        plant_out[:, ii] .+= noise_plant_out[:, ii]

        # Controller
        if ii >= ctrl_start_ind && ii < length(target)
            ctrl_out[:, ii+1] =
                ctrl(; ii, target, plant_inp, plant_out, act, inp_feedforward
                ) .+ noise_ctrl_out[:, ii+1]
        end
    end
    return Dict(
        :plant => plant,
        :act => act,
        :ctrl => ctrl,
        :target => target,
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
