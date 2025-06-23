import ControlSystemsBase: lsim, StateSpace
import ControlSystemIdentification: iddata, newpem, PredictionStateSpace
import LinearAlgebra: I, inv
import LsqFit: curve_fit, coef

export offset_scale,
    unscale_unoffset, model_evolve, system_id,
    system_id_optimal_inp_cond

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

"""
    system_id(
        inp::Union{Vector{Float64}, Matrix{Float64}},
        out::Union{Vector{Float64}, Matrix{Float64}},
        tt::Vector{Float64},
        order::Int;
        inp_cond::Union{Nothing, Function}=nothing,
        inp_offset::Float64=0.0, inp_factor::Float64=1.0,
        out_offset::Float64=0.0, out_factor::Float64=1.0,
        inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
        newpem_kwargs::Dict{Symbol, U}=Dict{Symbol, Any}(),
        verbose::Bool=false,
    ) where {T, U}

Perform system identification for a set on input data `inp`, output data `out`, and
time series vector `tt`. If there are more than one inputs or outputs, provide them as
Matrix with first dimension for ports (input or output) and second dimension for time.

If `inp_cond` is provided, it is applied to offseted and scaled input before
performing system identification with keywords for this function provided in
`inp_cond_kwargs`.

This function uses [ControlSystemIdentification.newpem](https://baggepinnen.github.io/ControlSystemIdentification.jl/stable/ss/#ControlSystemIdentification.newpem)
to perform the system identification. Any additional keywords for this function should
be passed as dictionary in `newpem_kwargs`. For advanced use, it is recommended to do
system identification directly with `newpem` instead of using this function.

Returns a linear state space model of the `order`.
"""
function system_id(
    inp::Union{Vector{Float64}, Matrix{Float64}},
    out::Union{Vector{Float64}, Matrix{Float64}},
    tt::Vector{Float64},
    order::Int;
    inp_cond::Union{Nothing, Function}=nothing,
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
    newpem_kwargs::Dict{Symbol, U}=Dict{Symbol, Any}(),
    verbose::Bool=false,
) where {T, U}
    @assert size(inp)[end] == size(out)[end] == length(tt)
    inp_so = offset_scale(inp; offset=inp_offset, factor=inp_factor)
    if !isnothing(inp_cond)
        inp_so = inp_cond(inp_so; inp_cond_kwargs...)
    end
    u = Matrix{Float64}[]
    if isa(inp_so, Vector)
        u = Matrix(inp_so')
    else
        u = inp_so
    end

    out_so = offset_scale(out; offset=out_offset, factor=out_factor)
    v = Matrix{Float64}[]
    if isa(out_so, Vector)
        v = Matrix(out_so')
    else
        v = out_so
    end

    id_data = iddata(v, u, tt[2] - tt[1])

    # newpem generates a lot of output when identifying the system and it can be a lot
    # when running this function iteratievely, verbose=false helps in that case.
    original_stdout = stdout
    if !verbose
        redirect_stdout(devnull)
    end
    sys = newpem(id_data, order; newpem_kwargs...).sys
    if !verbose
        redirect_stdout(original_stdout)
    end
    return sys
end

"""
    model_evolve(
        sys::Union{PredictionStateSpace, StateSpace},
        inp::Union{Vector{Float64}, Matrix{Float64}};
        x0::Union{Nothing, Vector{Float64}}=nothing,
        inp_cond::Union{Nothing, Function}=nothing,
        inp_offset::Float64=0.0, inp_factor::Float64=1.0,
        out_offset::Float64=0.0, out_factor::Float64=1.0,
        inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
    )::Union{Vector{Float64}, Matrix{Float64}} where {T}

Evolve a state space model `sys` with the input steps `inp`. The input is offseted and
scaled with `inp_offset` and `inp_factor` and the output is unscaled and unoffseted
with `out_offset` and `out_factor`. If a function is provided as `inp_cond` it is
applied to the input after scaling and offseting along with any keywords passed for it.
"""
function model_evolve(
    sys::Union{PredictionStateSpace, StateSpace},
    inp::Union{Vector{Float64}, Matrix{Float64}};
    x0::Union{Nothing, Vector{Float64}}=nothing,
    inp_cond::Union{Nothing, Function}=nothing,
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
)::typeof(inp) where {T}
    inp_so = offset_scale(inp; offset=inp_offset, factor=inp_factor)
    if !isnothing(inp_cond)
        inp_so = inp_cond(inp_so; inp_cond_kwargs...)
    end
    u = Matrix{Float64}[]
    if isa(inp_so, Vector)
        u = Matrix(inp_so')
    else
        u = inp_so
    end

    if isnothing(x0)
        x0 = inv(Matrix{Float64}(I, size(sys.A)...) - sys.A) * sys.B * u[:, 1]
    end
    modeled_out_so, _, x0, _ = lsim(sys, u; x0)
    modeled_out = unscale_unoffset(modeled_out_so; offset=out_offset, factor=out_factor)
    if size(modeled_out)[1] == 1
        return modeled_out[1, :]
    else
        return modeled_out
    end
end

"""
    system_id_optimal_inp_cond(
        inp::Union{Vector{Float64}, Matrix{Float64}},
        out::Union{Vector{Float64}, Matrix{Float64}},
        tt::Vector{Float64},
        order::Int,
        inp_cond::Union{Nothing, Function},
        inp_cond_args_guess::Dict{Symbol, T};
        inp_offset::Float64=0.0, inp_factor::Float64=1.0,
        out_offset::Float64=0.0, out_factor::Float64=1.0,
        inp_cond_args_lower::Dict{Symbol, V}=Dict{Symbol, Any}(),
        inp_cond_args_upper::Dict{Symbol, W}=Dict{Symbol, Any}(),
        newpem_kwargs::Dict{Symbol, U}=Dict{Symbol, Any}(),
        verbose::Bool=false,
    ) where {T, U, V, W}

Perform system identification for a set on input data `inp`, output data `out`, and
time series vector `tt`. If there are more than one inputs or outputs, provide them as
Matrix with first dimension for ports (input or output) and second dimension for time.

The `inp_cond` is applied to offseted and scaled input before performing system
identification. The `inp_cond_args_guess` is used as initial keyword arguments that
provide the parameters of the `inp_cond`. These arguments are then used to find the
best fit while iteratively performing [`system_id`](@ref) in each step.

This function uses [ControlSystemIdentification.newpem](https://baggepinnen.github.io/ControlSystemIdentification.jl/stable/ss/#ControlSystemIdentification.newpem)
to perform the system identification. Any additional keywords for this function should
be passed as dictionary in `newpem_kwargs`.

This function uses [LsqFit.curve_fit](https://julianlsolvers.github.io/LsqFit.jl/latest/api/#LsqFit.curve_fit)
to fit the parameters of input conditions along with performing the system
identification.

For advanced use, it is recommended to do system identification directly with `newpem`
and optimize using your favorite fitting method instead of using this function.

Returns a linear state space model of the `order` and the keyword argument dictionary
containing optimal parameters for `inp_cond` function.
"""
function system_id_optimal_inp_cond(
    inp::Union{Vector{Float64}, Matrix{Float64}},
    out::Union{Vector{Float64}, Matrix{Float64}},
    tt::Vector{Float64},
    order::Int,
    inp_cond::Function,
    inp_cond_args_guess::Dict{Symbol, T};
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    inp_cond_args_lower::Dict{Symbol, V}=Dict{Symbol, Any}(),
    inp_cond_args_upper::Dict{Symbol, W}=Dict{Symbol, Any}(),
    newpem_kwargs::Dict{Symbol, U}=Dict{Symbol, Any}(),
    verbose::Bool=false,
) where {T, U, V, W}
    key_list = collect(keys(inp_cond_args_guess))
    function model(inp, param)
        inp_cond_kwargs = Dict(key => param[ii] for (ii, key) ∈ enumerate(key_list))
        sys = system_id(
            inp, out, tt, order;
            inp_cond, inp_offset, inp_factor, out_offset, out_factor,
            inp_cond_kwargs,
            newpem_kwargs, verbose,
        )
        return model_evolve(
            sys, inp;
            inp_cond, inp_offset, inp_factor, out_offset, out_factor, inp_cond_kwargs,
        )
    end
    guess = [inp_cond_args_guess[key] for key ∈ key_list]
    lower = [
        inp_cond_args_lower[key] for
        key ∈ key_list if key in keys(inp_cond_args_lower)
    ]
    lower = [
        inp_cond_args_upper[key] for
        key ∈ key_list if key in keys(inp_cond_args_upper)
    ]

    fit_result = coef(curve_fit(model, inp, out, guess))
    opt = Dict{Symbol, T}(key => fit_result[ii] for (ii, key) ∈ enumerate(key_list))

    sys = system_id(
        inp, out, tt, order;
        inp_cond, inp_offset, inp_factor, out_offset, out_factor,
        inp_cond_kwargs=opt,
        newpem_kwargs, verbose,
    )
    return sys, opt
end
