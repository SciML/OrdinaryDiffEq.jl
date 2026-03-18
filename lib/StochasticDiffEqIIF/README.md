# StochasticDiffEqIIF.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

StochasticDiffEqIIF.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo providing integrating factor methods for stochastic differential equations with semi-linear structure. While completely independent and usable on its own, users wanting all SDE solvers should use [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

## Solvers

- `IIF1M` - Implicit Integrating Factor Euler-Maruyama
- `IIF2M` - Implicit Integrating Factor order 2
- `IIF1Mil` - Implicit Integrating Factor Milstein

These methods are designed for SDEs of the form `du = (Au + f(u,p,t))dt + g(u,p,t)dW` where `A` is a linear operator, leveraging the integrating factor `exp(At)` to handle the linear part exactly.

## Example

```julia
using StochasticDiffEqIIF, LinearAlgebra
A = [-1.0 0.0; 0.0 -2.0]
f = SplitFunction(DiffEqArrayOperator(A), (du, u, p, t) -> du .= 0.1 .* u)
g(du, u, p, t) = du .= 0.1 .* u
u0 = [1.0, 1.0]
tspan = (0.0, 1.0)
prob = SDEProblem(f, g, u0, tspan)
sol = solve(prob, IIF1M(), dt = 0.01)
```

## Available Solvers

For the full list of available SDE solvers, refer to the [DifferentialEquations.jl SDE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/) documentation.
