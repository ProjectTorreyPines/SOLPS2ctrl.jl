import ControlSystemIdentification: PredictionStateSpace
import ControlSystemsBase: StateSpace

export lsim_step, offset_scale, unscale_unoffset

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
    offset_scale(
        val::Union{Float64, Vector{Float64}, Matrix{Float64}};
        offset::Union{Float64, Vector{Float64}}=0.0,
        factor::Union{Float64, Vector{Float64}}=1.0,
    )::typeof(val)

Subtract an `offset` and multiply by a `factor`, the `val` to make it nominally in the
range of -1 to 1 (not strictly) for easy identification of system.

    (val .- offset) .* factor
"""
function offset_scale(
    val::Union{Float64, Vector{Float64}, Matrix{Float64}};
    offset::Union{Float64, Vector{Float64}}=0.0,
    factor::Union{Float64, Vector{Float64}}=1.0,
)::typeof(val)
    return (val .- offset) .* factor
end

"""
    unscale_unoffset(
        offset_scaled::Union{Float64, Vector{Float64}, Matrix{Float64}};
        offset::Union{Float64, Vector{Float64}}=0.0,
        factor::Union{Float64, Vector{Float64}}=1.0,
    )::typeof(offset_scaled)

Undo previously applied offset and scaling.

    offset_scaled ./ factor .+ offset
"""
function unscale_unoffset(
    offset_scaled::Union{Float64, Vector{Float64}, Matrix{Float64}};
    offset::Union{Float64, Vector{Float64}}=0.0,
    factor::Union{Float64, Vector{Float64}}=1.0,
)::typeof(offset_scaled)
    return offset_scaled ./ factor .+ offset
end
