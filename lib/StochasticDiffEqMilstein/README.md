# StochasticDiffEqMilstein.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

StochasticDiffEqMilstein.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo providing generalized Milstein methods for stochastic differential equations that support non-diagonal noise via Levy area approximation. While completely independent and usable on its own, users wanting all SDE solvers should use [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

## Solvers

- `RKMilGeneral` - Generalized Runge-Kutta Milstein (strong order 1.0, supports non-commutative noise via Levy area)
- `WangLi3SMil_A` through `WangLi3SMil_F` - Wang-Li three-stage Milstein methods

## Example

```julia
using StochasticDiffEqMilstein
f(u, p, t) = 1.01 * u
g(u, p, t) = 0.87 * u
u0 = 0.5
tspan = (0.0, 1.0)
prob = SDEProblem(f, g, u0, tspan)
sol = solve(prob, RKMilGeneral(), dt = 0.01)
```

## Available Solvers

For the full list of available SDE solvers, refer to the [DifferentialEquations.jl SDE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/) documentation.
