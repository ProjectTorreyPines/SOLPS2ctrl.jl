import ControlSystemsBase: lsim, StateSpace
import ControlSystemIdentification: iddata, newpem, PredictionStateSpace
import LinearAlgebra: I, inv
import LsqFit: curve_fit, coef

export model_evolve, system_id, system_id_optimal_inp_cond

"""
    system_id(
        inp::Union{Vector{Float64}, Matrix{Float64}},
        out::Union{Vector{Float64}, Matrix{Float64}},
        tt::Vector{Float64},
        order::Int;
        inp_offset::Float64=0.0, inp_factor::Float64=1.0,
        out_offset::Float64=0.0, out_factor::Float64=1.0,
        inp_cond::Union{Nothing, InputConditioning}=nothing,
        inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
        newpem_kwargs::Dict{Symbol, U}=Dict{Symbol, Any}(),
        verbose::Bool=false,
    ) where {T, U}

Perform system identification for a set on input data `inp`, output data `out`, and
time series vector `tt`. If there are more than one inputs or outputs, provide them as
Matrix with first dimension for ports (input or output) and second dimension for time.

If `inp_cond` is provided, it is applied before offsetting and scaling for
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
    inp_offset::Float64=0.0, inp_factor::Float64=1.0,
    out_offset::Float64=0.0, out_factor::Float64=1.0,
    inp_cond::Union{Nothing, InputConditioning}=nothing,
    inp_cond_kwargs::Dict{Symbol, T}=Dict{Symbol, Any}(),
    newpem_kwargs::Dict{Symbol, U}=Dict{Symbol, Any}(),
    verbose::Bool=false,
) where {T, U}
    @assert size(inp)[end] == size(out)[end] == length(tt)
    if isa(inp, Vector)
        inp_m = Matrix(inp')
    else
        inp_m = inp
    end
    inp_co = zeros(Float64, size(inp_m))
    if !isnothing(inp_cond)
        temp_inst = deepcopy(inp_cond)
        for ii ∈ eachindex(tt)
            inp_co[:, ii] = temp_inst(inp_m[:, ii]; inp_cond_kwargs...)
        end
    end
    u = offset_scale(inp_co; offset=inp_offset, factor=inp_factor)

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
    linear_plant = LinearPlant(sys; inp_offset, inp_factor, out_offset, out_factor)
    if isnothing(inp_cond)
        return linear_plant
    else
        return InpCondLinPlant(linear_plant, inp_cond, inp_cond_kwargs)
    end
end

"""
    model_evolve(
        plant::Plant,
        inp::Union{Vector{Float64}, Matrix{Float64}};
        x0::Union{Nothing, Vector{Float64}}=nothing,
        initialize_x::Bool=false,
    )

Evolve a plant model through a time series of inputs. If `inp` is a Matrix, the first
dimension will hold different inputs and second dimension will be along time. If `inp`
is a vector, it would be assumed that the model has a single input. If the model also
happens to have a single output and `inp` is a vector, the returned outputs will be a
vector as well. If `x0` is not provided, the stored state vector of the plant model will
be used for initialization. If `initialize_x` is true, then the state vector of the
plant model would be initialized to match first input value so that no sudden jumps
happen.
"""
function model_evolve(
    plant::Plant,
    inp::Union{Vector{Float64}, Matrix{Float64}};
    x0::Union{Nothing, Vector{Float64}}=nothing,
    initialize_x::Bool=false,
)
    if isa(inp, Vector)
        inp_m = Matrix(inp')
    else
        inp_m = inp
    end
    temp_plant = deepcopy(plant)

    lp = get_linear_plant(temp_plant)

    if initialize_x
        sys = lp.sys
        initial_inp =
            offset_scale(inp_m[:, 1]; offset=lp.inp_offset, factor=lp.inp_factor)
        x0 = inv(Matrix{Float64}(I, size(sys.A)...) - sys.A) * sys.B * initial_inp
    end

    if !isnothing(x0)
        lp.x = x0
    end

    out_m = zeros(Float64, size(lp.sys.C, 1), size(inp_m, 2))
    for ii ∈ axes(inp_m, 2)
        out_m[:, ii] = temp_plant(inp_m[:, ii])
    end
    if isa(inp, Vector) && size(out_m, 1) == 1
        return out_m[1, :]
    else
        return out_m
    end
end

"""
    system_id_optimal_inp_cond(
        inp::Union{Vector{Float64}, Matrix{Float64}},
        out::Union{Vector{Float64}, Matrix{Float64}},
        tt::Vector{Float64},
        order::Int,
        inp_cond::InputConditioning,
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

The `inp_cond` is applied before offsetting and scaling for performing system
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

Returns a NonLinInpLinearPlant that has stored best fit parameters for the input
non-linearity and a best fit linear model.
"""
function system_id_optimal_inp_cond(
    inp::Union{Vector{Float64}, Matrix{Float64}},
    out::Union{Vector{Float64}, Matrix{Float64}},
    tt::Vector{Float64},
    order::Int,
    inp_cond::InputConditioning,
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
        plant = system_id(
            inp, out, tt, order;
            inp_offset, inp_factor, out_offset, out_factor,
            inp_cond, inp_cond_kwargs,
            newpem_kwargs, verbose,
        )
        return model_evolve(plant, inp; initialize_x=true)
    end
    guess = [inp_cond_args_guess[key] for key ∈ key_list]
    lower = Vector{Float64}(
        [
        inp_cond_args_lower[key] for
        key ∈ key_list if key in keys(inp_cond_args_lower)
    ]
    )
    upper = Vector{Float64}(
        [
        inp_cond_args_upper[key] for
        key ∈ key_list if key in keys(inp_cond_args_upper)
    ]
    )

    fit_result = coef(curve_fit(model, inp, out, guess; lower, upper))
    opt = Dict{Symbol, T}(key => fit_result[ii] for (ii, key) ∈ enumerate(key_list))

    return system_id(
        inp, out, tt, order;
        inp_cond, inp_offset, inp_factor, out_offset, out_factor,
        inp_cond_kwargs=opt,
        newpem_kwargs, verbose,
    )
end
