"""
    KuttaPRK2p5(; threading = true)

Fifth-order explicit parallel Runge-Kutta method with five stages arranged for two
processors. It is useful when each right-hand-side evaluation is sufficiently expensive
to benefit from the available stage parallelism.

# Fields

  - `threading`: threading configuration used for the parallel stages. It must be a
    `Bool` or an OrdinaryDiffEqCore threading option such as
    `OrdinaryDiffEqPRK.BaseThreads()`.

# Keywords

  - `threading = true`: threading configuration for the parallel stages. `true` uses
    `Threads.@threads`, `false` is serial, and
    [`BaseThreads()`](@ref OrdinaryDiffEqCore.BaseThreads) or
    [`PolyesterThreads()`](@ref OrdinaryDiffEqCore.PolyesterThreads) select the
    corresponding OrdinaryDiffEqCore backend. These are public but not exported, so
    write them qualified (`OrdinaryDiffEqPRK.BaseThreads()`) or bring them into
    scope with `using OrdinaryDiffEqCore: BaseThreads, PolyesterThreads`.
    `PolyesterThreads` additionally requires `using Polyester`.

# Examples

```julia
using OrdinaryDiffEq

prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
sol = solve(prob, KuttaPRK2p5())
```

# References

K. R. Jackson and S. P. Norsett, "The potential for parallelism in Runge--Kutta
methods. Part 1: RK formulas in standard form," *SIAM Journal on Numerical Analysis*,
32(1), 49-82 (1995).
"""
Base.@kwdef struct KuttaPRK2p5{TO} <: OrdinaryDiffEqAlgorithm
    threading::TO = true
end
