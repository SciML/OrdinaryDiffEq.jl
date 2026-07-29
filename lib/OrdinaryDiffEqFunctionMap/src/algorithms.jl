@doc raw"""
    FunctionMap(; scale_by_time = false)

A fixed timestep method for when the ODE is a discrete dynamical system. In the 
operator setting, this is equivalent to operator splitting for additive operators.

When `scale_by_time = true`, the method becomes `u_{n+1} = u_n + dt*f(u_n,p,t_n)`, 
otherwise it's `u_{n+1} = f(u_n,p,t_n)`.

!!! note
    This method requires a fixed timestep dt and is not adaptive.
"""
struct FunctionMap{scale_by_time, StepLimiter} <: OrdinaryDiffEqAlgorithm
    step_limiter!::StepLimiter
end
FunctionMap{scale_by_time}() where {scale_by_time} = FunctionMap{scale_by_time, typeof(trivial_limiter!)}(
    trivial_limiter!
)
function FunctionMap(; scale_by_time = false, step_limiter = trivial_limiter!, kwargs...)
    old_kw = Symbol("step_limiter!")
    if haskey(kwargs, old_kw)
        if step_limiter === trivial_limiter!
            step_limiter = get(kwargs, old_kw, trivial_limiter!)
        end
    end
    unsupported_kwargs = filter(!=(old_kw), keys(kwargs))
    if !isempty(unsupported_kwargs)
        throw(ArgumentError("Unsupported keyword argument(s): $(unsupported_kwargs)"))
    end

    return FunctionMap{scale_by_time, typeof(step_limiter)}(step_limiter)
end
