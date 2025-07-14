import DataStructures: Queue, enqueue!, dequeue!

export Actuator, DelayedActuator

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
