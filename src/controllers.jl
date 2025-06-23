import ControlSystemIdentification: PredictionStateSpace
import ControlSystemsBase: StateSpace, Discrete
import LinearAlgebra: pinv
using DataStructures: Queue, enqueue!, dequeue!
# import ControlSystemsBase: lsim, ss, pid, d2c_exact, d2c, StateSpace

export lsim_step, state_prediction_matrices

function lsim_step(
    sys::Union{PredictionStateSpace, StateSpace}, u::Vector{Float64};
    x0::Vector{Float64}=zeros(size(sys.A, 1)),
)::Tuple{Vector{Float64}, Vector{Float64}}
    x = sys.A * x0 + sys.B * u
    y = sys.C * x0 + sys.D * u
    return y, x
end

function model_step(
    sys::Union{PredictionStateSpace, StateSpace},
    inp::Vector{Float64};
    x0::Vector{Float64}=zeros(size(sys.A, 1)),
    inp_cond::Union{Nothing, Function}=nothing,
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
)::Vector{Float64} where {T}
    inp_so = offset_scale(inp; offset=inp_offset, factor=inp_factor)
    if !isnothing(inp_cond)
        inp_so = inp_cond(inp_so; inp_cond_kwargs...)
    end
    modeled_out_so, _, x0, _ = lsim_step(sys, u; x0)
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

# function latency_jump_part(
#     future_flow_rates::Vector{Float64}, ssm::PredictionStateSpace,
#     fut_x0::Vector{Float64};
#     gas_offset::Float64=gas_offset, gas_factor::Float64=gas_factor,
#     boltzmann_k::Float64=1.380649e-23, Temp::Float64=300.0,
# )::Tuple{Float64, Vector{Float64}}
#     fut_gas_atomps = gas_SI_to_aps(future_flow_rates; boltzmann_k, Temp)
#     fut_gas_atomsps_scaled = zeros((1, length(fut_gas_atomps)))
#     fut_gas_atomsps_scaled[1, :] =
#         gas_offset_scale(fut_gas_atomps; gas_offset, gas_factor)
#     y_int, _, post_lat_x0, _ = lsim(ssm, fut_gas_atomsps_scaled; x0=fut_x0)
#     fut_meas_hf = y_int[1, end]
#     return fut_meas_hf, post_lat_x0[:, end]
# end

abstract type Actuator end

mutable struct DelayedActuator{U} <: Actuator
    delay::Int
    default_output::U
    buffer::Queue{U}
    function DelayedActuator(
        delay::Int;
        default_output::T=0.0,
    ) where {T <: Union{Float64, Vector{Float64}}}
        return new{T}(delay, default_output, Queue{T}(delay))
    end
end

function (da::DelayedActuator)(inp::Union{Float64, Vector{Float64}})
    enqueue!(da.buffer, inp)
    if length(da.buffer) <= da.delay
        return da.default_output
    else
        return dequeue!(da.buffer)
    end
end

get_T(::DelayedActuator{T}) where {T} = T

function has_delay(act::Actuator)
    for fieldname ∈ fieldnames(act)
        if isa(fieldtype(act, fieldname), DelayedActuator)
            return true
        end
    end
    return false
end

function get_future_inp(act::Actuator)
    for fieldname ∈ fieldnames(act)
        if isa(fieldtype(act, fieldname), DelayedActuator)
            return collect(getfield(act, fieldname).buffer)
        end
    end
    return Matrix{Float64}(undef, 0, 0)
end

abstract type Controller end

mutable struct linear_controller <: Controller
    ctrl_ss::StateSpace{Discrete}
    ctrl_x0::Vector{Float64}
end

function (lc::linear_controller)(;
    ii::Int,
    target::Matrix{Float64},
    plant_out::Matrix{Float64},
    kwargs...,
)::Vector{Float64}
    err = target[:, ii] .- plant_out[:, ii]
    ctrl_out, lc.ctrl_x0 = lsim_step(lc.ctrl_ss, err; x0=lc.ctrl_x0)
    return ctrl_out
end

mutable struct PVLC <: Controller
    ctrl_ss::StateSpace{Discrete}
    ctrl_x0::Vector{Float64}
    h::Int
    plant_model::Union{PredictionStateSpace, StateSpace}
    Y2x::Matrix{Float64}
    U2x::Matrix{Float64}
    function PVLC(
        ctrl_ss::StateSpace{Discrete},
        ctrl_x0::Vector{Float64},
        plant_model::Union{PredictionStateSpace, StateSpace},
        h::Int=size(sys.A, 1),
    )
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
    ctrl_out = zeros(Float64, size(ctrl_ss.C, 1))
    if length(plant_inp) >= size(pvlc.U2x, 2)
        # Gather history of inputs and outputs
        U = vec(plant_inp[:, (end-pvlc.h+1):pvlc.h])
        Y = vec(plant_out[:, (end-pvlc.h+1):pvlc.h])
        # Estimate current state vector of plant model
        est_x = pvlc.Y2x * Y .+ U2x * U
        # Evolve into the future for the delay in actuator
        future_out, _ = model_evolve(
            plant_model, future_inp;
            x0=est_x,
            inp_offset, inp_factor, out_offset, out_factor,
        )
        # Estimate the error signal in future
        err = target[:, ii+size(future_inp, 2)] - future_out[:, end]
        # Apply the linear controller in future
        ctrl_out, pvlc.ctrl_x0 = lsim_step(pvlc.ctrl_ss, err; x0=pvlc.ctrl_x0)
    end
    return ctrl_out
end

function run_closed_loop_sim(
    plant_model::Union{PredictionStateSpace, StateSpace},
    act::Actuator,
    ctrl::Controller,
    target::Union{Vector{Float64}, Matrix{Float64}};
    inp_cond::Union{Nothing, Function}=nothing,
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
    inp_feedforward::Matrix{Float64}=zeros(
        Float64,
        (size(plant_model.B, 2), length(target)),
    ),
    ctrl_start_ind::Int=1,
) where {T}
    if isa(target, Vector)
        if size(plant_model.C, 1) == 1
            target = Matrix(target')
        else
            error("Plant has multiple outputs, target must be Matrix")
        end
    end
    ctrl_out = zeros(Float64, size(inp_feedforward, 1))
    plant_inp = zeros(Float64, (size(plant_model.B, 2), length(target)))
    plant_out = zeros(Float64, (size(plant_model.C, 1), length(target)))
    x0 =
        inv(Matrix{Float64}(I, size(plant_model.A)...) - plant_model.A) *
        plant_model.B * inp_feedforward[:, 1]
    for ii ∈ eachindex(target)
        plant_inp[:, ii] = act(ctrl_out .+ inp_feedforward[:, ii])
        plant_out[:, ii], _ = model_step(
            plant_model, plant_inp[:, ii];
            x0,
            inp_cond, inp_offset, inp_factor, out_offset, out_factor,
            inp_cond_kwargs,
        )
        # err = target[:, ii] .- plant_out[:, ii]
        future_inp = get_future_inp(act)
        if ii >= ctrl_start_ind
            ctrl_out = ctrl(;
                ii, target, plant_inp, plant_out, future_inp,
                inp_offset, inp_factor, out_offset, out_factor,
            )
        end
    end
end
