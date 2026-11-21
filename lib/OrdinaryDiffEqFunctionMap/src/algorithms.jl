"""
    FunctionMap(; scale_by_time = false, step_limiter = trivial_limiter!)

Fixed-step algorithm for discrete dynamical systems. By default, each step applies the
problem function directly, ``u_{n + 1} = f(u_n, p, t_n)``. With `scale_by_time = true`,
it instead applies an explicit-Euler update, ``u_{n + 1} = u_n + dt * f(u_n, p, t_n)``.

`FunctionMap` is non-adaptive. Supply `dt` when using `scale_by_time = true` with an ODE
problem.

# Fields

  - `step_limiter!`: function applied after each completed step. It receives the ordinary
    DifferentialEquations limiter arguments and can enforce problem-specific constraints.

# Keywords

  - `scale_by_time = false`: choose the direct map update. Set to `true` for the
    explicit-Euler update.
  - `step_limiter = trivial_limiter!`: post-step limiter. The legacy `step_limiter!`
    keyword is also accepted for compatibility; prefer `step_limiter` in new code.

# Throws

  - `ArgumentError`: unsupported keyword arguments are supplied.

# Examples

```julia
using OrdinaryDiffEqFunctionMap
using SciMLBase: DiscreteProblem, solve

prob = DiscreteProblem((u, p, t) -> 0.5u, 1.0, (0.0, 3.0))
sol = solve(prob, FunctionMap())
```
"""
struct FunctionMap{scale_by_time, StepLimiter} <: OrdinaryDiffEqAlgorithm
    step_limiter!::StepLimiter
end
FunctionMap{scale_by_time}() where {scale_by_time} = FunctionMap{scale_by_time, typeof(trivial_limiter!)}(
    trivial_limiter!
)
function FunctionMap(; scale_by_time = false, step_limiter = trivial_limiter!, kwargs...)
    kwargs_nt = values(kwargs)
    old_kw = Symbol("step_limiter!")
    if haskey(kwargs_nt, old_kw)
        if step_limiter === trivial_limiter!
            step_limiter = get(kwargs_nt, old_kw, trivial_limiter!)
        end
    end
    extra_kwargs = Base.structdiff(kwargs_nt, NamedTuple{(old_kw,)})
    if !isempty(extra_kwargs)
        throw(ArgumentError("Unsupported keyword argument(s): $(keys(extra_kwargs))"))
    end

    return FunctionMap{scale_by_time, typeof(step_limiter)}(step_limiter)
end
