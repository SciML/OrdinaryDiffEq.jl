"""
    Feagin10(; step_limiter! = trivial_limiter!) -> Feagin10

Tenth-order explicit Runge-Kutta method for high-accuracy integration of non-stiff
ordinary differential equations. Its embedded error estimate has order eight.

# Fields

  - `step_limiter!`: deprecated algorithm-level step limiter. It is called as
    `step_limiter!(u, integrator, p, t)` after each accepted step when no solve-level
    limiter is supplied.

# Keywords

  - `step_limiter! = trivial_limiter!`: step limiter stored in the algorithm. Prefer the
    `step_limiter` keyword of `solve` or `init`; the constructor keyword is retained for
    compatibility.

# Returns

  - A `Feagin10` algorithm instance for use with `solve` or `init`.

# Examples

```julia
using OrdinaryDiffEqFeagin: Feagin10
using SciMLBase: ODEProblem, solve

prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
sol = solve(prob, Feagin10(); abstol = 1.0e-12, reltol = 1.0e-12)
```

# References

T. Feagin, "High-order explicit Runge-Kutta methods using m-symmetry,"
*Neural, Parallel & Scientific Computations* (2012).
"""
Base.@kwdef struct Feagin10{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
    step_limiter!::StepLimiter = trivial_limiter!
end

"""
    Feagin12(; step_limiter! = trivial_limiter!) -> Feagin12

Twelfth-order explicit Runge-Kutta method for high-accuracy integration of non-stiff
ordinary differential equations.

# Fields

  - `step_limiter!`: deprecated algorithm-level step limiter. It is called as
    `step_limiter!(u, integrator, p, t)` after each accepted step when no solve-level
    limiter is supplied.

# Keywords

  - `step_limiter! = trivial_limiter!`: step limiter stored in the algorithm. Prefer the
    `step_limiter` keyword of `solve` or `init`; the constructor keyword is retained for
    compatibility.

# Returns

  - A `Feagin12` algorithm instance for use with `solve` or `init`.

# Examples

```julia
using OrdinaryDiffEqFeagin: Feagin12
using SciMLBase: ODEProblem, solve

prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
sol = solve(prob, Feagin12(); abstol = 1.0e-12, reltol = 1.0e-12)
```

# References

T. Feagin, "High-order explicit Runge-Kutta methods using m-symmetry,"
*Neural, Parallel & Scientific Computations* (2012).
"""
Base.@kwdef struct Feagin12{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
    step_limiter!::StepLimiter = trivial_limiter!
end

"""
    Feagin14(; step_limiter! = trivial_limiter!) -> Feagin14

Fourteenth-order explicit Runge-Kutta method for high-accuracy integration of non-stiff
ordinary differential equations. Its embedded error estimate has order twelve.

# Fields

  - `step_limiter!`: deprecated algorithm-level step limiter. It is called as
    `step_limiter!(u, integrator, p, t)` after each accepted step when no solve-level
    limiter is supplied.

# Keywords

  - `step_limiter! = trivial_limiter!`: step limiter stored in the algorithm. Prefer the
    `step_limiter` keyword of `solve` or `init`; the constructor keyword is retained for
    compatibility.

# Returns

  - A `Feagin14` algorithm instance for use with `solve` or `init`.

# Examples

```julia
using OrdinaryDiffEqFeagin: Feagin14
using SciMLBase: ODEProblem, solve

prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
sol = solve(prob, Feagin14(); abstol = 1.0e-12, reltol = 1.0e-12)
```

# References

T. Feagin, "An explicit Runge-Kutta method of order fourteen," *Numerical Algorithms*
(2009).
"""
Base.@kwdef struct Feagin14{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
    step_limiter!::StepLimiter = trivial_limiter!
end
