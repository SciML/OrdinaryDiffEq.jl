"""
    QPRK98(; stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!,
        thread = Serial())

Ninth-order explicit embedded Runge-Kutta pair with an eighth-order error estimate,
designed for quadruple-precision computations. Its parallel tableau is most useful when
function evaluations are expensive and extended-precision accuracy is required.

# Fields

  - `stage_limiter!`: function applied after every Runge-Kutta stage.
  - `step_limiter!`: function applied after every accepted step.
  - `thread`: OrdinaryDiffEqCore threading option for independent stage work.

# Keywords

  - `stage_limiter! = trivial_limiter!`: stage limiter with the OrdinaryDiffEq limiter
    calling convention. Use this to enforce stage-wise constraints.
  - `step_limiter! = trivial_limiter!`: limiter applied to the completed step. Use this
    to enforce solution constraints after an accepted update.
  - `thread = Serial()`: threading configuration. Use `OrdinaryDiffEqCore.BaseThreads()`
    or `OrdinaryDiffEqCore.PolyesterThreads()` to enable the corresponding
    OrdinaryDiffEqCore backend. These types are public but not exported, so they need
    either the `OrdinaryDiffEqCore.` qualification shown here or a `using
    OrdinaryDiffEqCore: BaseThreads, PolyesterThreads`.

# Examples

```julia
using OrdinaryDiffEqQPRK
using SciMLBase: ODEProblem, solve

prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
sol = solve(prob, QPRK98())
```

# References

V. N. Kovalnogov, R. V. Fedorov, T. V. Karpukhina, T. E. Simos, and C. Tsitouras,
"Runge-Kutta pairs of orders 9(8) for use in quadruple precision computations,"
*Numerical Algorithms* (2023). https://doi.org/10.1007/s11075-023-01632-8
"""
Base.@kwdef struct QPRK98{StageLimiter, StepLimiter, Thread} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = Serial()
end

OrdinaryDiffEqCore.has_stage_limiter(::QPRK98) = true
