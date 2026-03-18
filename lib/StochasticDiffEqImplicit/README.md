# StochasticDiffEqImplicit.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

StochasticDiffEqImplicit.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo providing drift-implicit and split-step implicit methods for stiff stochastic differential equations. While completely independent and usable on its own, users wanting all SDE solvers should use [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

## Solvers

### Drift-Implicit Methods
- `ImplicitEM` - Implicit Euler-Maruyama (strong order 0.5)
- `ImplicitEulerHeun` - Implicit Euler-Heun (Stratonovich interpretation)
- `ImplicitRKMil` - Implicit Runge-Kutta Milstein (strong order 1.0, diagonal noise)

### Split-Step Methods
- `STrapezoid` - Stochastic Trapezoidal rule
- `SImplicitMidpoint` - Stochastic Implicit Midpoint rule
- `ISSEM` - Implicit Split-Step Euler-Maruyama
- `ISSEulerHeun` - Implicit Split-Step Euler-Heun (Stratonovich)

### Stochastic SDIRK
- `SKenCarp` - Stochastic Kennedy-Carpenter ESDIRK method (strong order 1.5 for additive noise)

## Example

```julia
using StochasticDiffEqImplicit
f(u, p, t) = -100.0 * u  # stiff drift
g(u, p, t) = 0.1 * u
u0 = 1.0
tspan = (0.0, 1.0)
prob = SDEProblem(f, g, u0, tspan)
sol = solve(prob, ImplicitEM(), dt = 0.01)
```

## Available Solvers

For the full list of available SDE solvers, refer to the [DifferentialEquations.jl SDE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/) documentation.
