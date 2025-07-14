import ControlSystemIdentification: PredictionStateSpace
import ControlSystemsBase: StateSpace

export Plant,
    LinearPlant, InputConditioning, InpCondLinPlant, get_linear_plant, get_sys, get_x,
    set_x!

"""
    Plant

Abstract parent type for all plants. Whenever a user defined plant is created,
it must be subtype of `Plant`.

To create a new plant with custom feautres, it must be defined as a mutable stucture
which is daughter of `Plant` that contains all settings and state information for
the plant and the instance itself should be callable to take as input a
`Vector{Float64}` and should return an output of `Vector{Float64}`.

```julia
mutable struct CustomPlant <: Plant
    settings
    state
    # ... Anything more
end

function (cp::CustomPlant)(inp::Vector{Float64})::Vector{Float64}
    # perform the plant single step forward calcualtion with inp
    # update cp.state if required
    return output
end
```

**NOTE:** If you need to add a linear system in your plant model, add a
[`LinearPlant`](@ref) instance in the attributes of your `CustomPlant` and just call
the LinearPlant inside your function call.
"""
abstract type Plant end

"""
    LinearPlant

Implementation of linear system with option of scaling and offseting inputs and outputs.
It stores the system in `sys` and state of the system in `x`.
Constructor:

    LinearPlant(
        sys::Union{PredictionStateSpace, StateSpace},
        x=zeros(Float64, size(sys.A, 1));
        inp_offset::Float64=0.0, inp_factor::Float64=1.0,
        out_offset::Float64=0.0, out_factor::Float64=1.0,
    )

Creates a LinearPlant instance with state space system `sys` and state vector `x`. It
also defined input offsetting and scaling with  `inp_offset` and `inp_factor`, and
similarly output unscaling and unoffseting with `out_offset` and `out_factor`.
"""
mutable struct LinearPlant <: Plant
    sys::Union{PredictionStateSpace, StateSpace}
    x::Vector{Float64}
    inp_offset::Float64
    inp_factor::Float64
    out_offset::Float64
    out_factor::Float64
    function LinearPlant(
        sys::Union{PredictionStateSpace, StateSpace},
        x=zeros(Float64, size(sys.A, 1));
        inp_offset::Float64=0.0, inp_factor::Float64=1.0,
        out_offset::Float64=0.0, out_factor::Float64=1.0,
    )
        return new(sys, x, inp_offset, inp_factor, out_offset, out_factor)
    end
end

function (lp::LinearPlant)(inp::Vector{Float64})::Vector{Float64}
    u = offset_scale(inp; offset=lp.inp_offset, factor=lp.inp_factor)
    modeled_out_so, lp.x = lsim_step(lp.sys, u; x0=lp.x)
    return unscale_unoffset(modeled_out_so; offset=lp.out_offset, factor=lp.out_factor)
end

"""
    InputConditioning

Abstract parent type for creating input conditioning to plants. Whenever a user
defined input coniditioning is created it must be subtype of `InputConditioning`.

To create a new input conditioning with custom feautres, it must be defined as a
mutable stucture which is daughter of `InputConditioning` that contains all settings
and state information for the function and the instance itself should be callable to
take as input a `Vector{Float64}` and kwargs for any parameters used and should return
an output of `Vector{Float64}` which is same size as the input.

```julia
mutable struct CustomInpCond <: InputConditioning
    settings
    state
    # ... Anything more
end

function (inp_cond::CustomInpCond)(inp::Vector{Float64}; kwargs...)::Vector{Float64}
    # perform the input coniditioning single step forward calcualtion with inp
    # Use parameters from kwargs
    # update inp_cond.state if required
    return output
end
```

Note that `inp_cond` must have a call signature
of `inp_cond(inp::Vector{Float64}; kwargs...)::Vector{Float64}` to take input as a
vector for multiple inputs at a particular time instance
"""
abstract type InputConditioning end

"""
    InpCondLinPlant

Implementation of a LinearPlant with an added input conditioning which could be used to
make it non-linear.
It stores the LinearPlant in `linear_plant`, the input coniditioning function in
`inp_cond` and any keyword arguments for the `inp_cond` in `inp_cond_kwargs`.
Constructor:

    InpCondLinPlant{T}(
        linear_plant::LinearPlant,
        inp_cond::InputConditioning
        inp_cond_kwargs::Dict{Symbol, T}
    ) where {T}

Creates a InpCondLinPlant instance. `inp_cond_kwargs` can be used to have changeable
parameters in the input conditioning which can be used for optimization and fiting.
"""
mutable struct InpCondLinPlant{T} <: Plant
    linear_plant::LinearPlant
    inp_cond::InputConditioning
    inp_cond_kwargs::Dict{Symbol, T}
end

function (nlp::InpCondLinPlant)(inp::Vector{Float64})::Vector{Float64}
    conditioned_inp = nlp.inp_cond(inp; nlp.inp_cond_kwargs...)
    return nlp.linear_plant(conditioned_inp)
end

function get_linear_plant(p::Plant)
    for fieldname ∈ fieldnames(typeof(p))
        if fieldtype(typeof(p), fieldname) == LinearPlant
            return getfield(p, fieldname)
        end
    end
    return error(
        "Plant must use an attribute of type LinearPlant to store the linear plant model",
    )
end

get_linear_plant(p::LinearPlant) = p

get_sys(p::Plant) = get_linear_plant(p).sys

get_x(p::Plant) = get_linear_plant(p).x

set_x!(p::Plant, x::Vector{Float64}) = setfield!(get_linear_plant(p), :x, x)
