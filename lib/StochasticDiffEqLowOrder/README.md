# StochasticDiffEqLowOrder.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

StochasticDiffEqLowOrder.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo providing low-order stochastic differential equation solvers. While completely independent and usable on its own, users wanting all SDE solvers should use [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

## Solvers

- `EM` - Euler-Maruyama (strong order 0.5)
- `EulerHeun` - Euler-Heun (Stratonovich interpretation, strong order 0.5)
- `LambaEM` - Adaptive Euler-Maruyama
- `LambaEulerHeun` - Adaptive Euler-Heun (Stratonovich)
- `SimplifiedEM` - Simplified Euler-Maruyama
- `SplitEM` - Split-step Euler-Maruyama for split problems
- `RKMil` - Runge-Kutta Milstein (strong order 1.0, diagonal noise)
- `RKMilCommute` - Milstein for commutative noise
- `PCEuler` - Predictor-corrector Euler

## Example

```julia
using StochasticDiffEqLowOrder
f(u, p, t) = 1.01 * u
g(u, p, t) = 0.87 * u
u0 = 0.5
tspan = (0.0, 1.0)
prob = SDEProblem(f, g, u0, tspan)
sol = solve(prob, EM(), dt = 0.01)
```

## Available Solvers

For the full list of available SDE solvers, refer to the [DifferentialEquations.jl SDE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/) documentation.
