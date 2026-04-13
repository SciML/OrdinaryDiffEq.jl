# StochasticDiffEqRODE.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/rode_solve/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

StochasticDiffEqRODE.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo providing solvers for random ordinary differential equations (RODEs) and Langevin dynamics. While completely independent and usable on its own, users wanting all SDE/RODE solvers should use [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

## Solvers

### RODE Solvers
- `RandomEM` - Random Euler-Maruyama for RODEProblems
- `RandomHeun` - Random Heun method for RODEProblems
- `RandomTamedEM` - Tamed Euler-Maruyama for RODEProblems with superlinear growth

### Langevin Dynamics
- `BAOAB` - BAOAB splitting integrator for Langevin dynamics (second-order ODEs with friction and noise)

## Example

```julia
using StochasticDiffEqRODE
f(u, p, t, W) = 2u * sin(W)
u0 = 1.0
tspan = (0.0, 1.0)
prob = RODEProblem(f, u0, tspan)
sol = solve(prob, RandomEM(), dt = 0.01)
```

## Available Solvers

For the full list of available RODE solvers, refer to the [DifferentialEquations.jl RODE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/rode_solve/) documentation.
