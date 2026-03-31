# StochasticDiffEqHighOrder.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

StochasticDiffEqHighOrder.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo providing high-order adaptive stochastic Runge-Kutta methods based on the Roessler SRI/SRA families. While completely independent and usable on its own, users wanting all SDE solvers should use [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

## Solvers

### SRI Family (Stochastic Runge-Kutta for Ito SDEs)
- `SRI` - General SRI framework (user-supplied tableau)
- `SRIW1` - SRI with weak order 2 (strong order 1.5 for scalar/diagonal noise)
- `SRIW2` - SRI variant
- `SOSRI` - Stability-optimized SRI
- `SOSRI2` - Stability-optimized SRI variant

### SRA Family (Stochastic Runge-Kutta for additive noise)
- `SRA` - General SRA framework (user-supplied tableau)
- `SRA1` - SRA with strong order 1.5 for additive noise
- `SRA2` - SRA variant
- `SRA3` - SRA variant
- `SOSRA` - Stability-optimized SRA
- `SOSRA2` - Stability-optimized SRA variant

## Example

```julia
using StochasticDiffEqHighOrder
f(u, p, t) = 1.01 * u
g(u, p, t) = 0.87 * u
u0 = 0.5
tspan = (0.0, 1.0)
prob = SDEProblem(f, g, u0, tspan)
sol = solve(prob, SOSRI())
```

## Available Solvers

For the full list of available SDE solvers, refer to the [DifferentialEquations.jl SDE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/) documentation.
