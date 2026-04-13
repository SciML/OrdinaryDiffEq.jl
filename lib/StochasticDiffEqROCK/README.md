# StochasticDiffEqROCK.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

StochasticDiffEqROCK.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo providing stabilized explicit stochastic Runge-Kutta-Chebyshev (ROCK) methods for stiff SDEs. While completely independent and usable on its own, users wanting all SDE solvers should use [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

## Solvers

- `SROCK1` - Stochastic ROCK with strong order 0.5
- `SROCK2` - Stochastic ROCK with weak order 2.0
- `KomBurSROCK2` - Komori-Burrage variant of SROCK2
- `SROCKC2` - SROCK with weak order 2.0 (Chebyshev variant)
- `SROCKEM` - SROCK Euler-Maruyama variant
- `SKSROCK` - SK-SROCK stabilized method
- `TangXiaoSROCK2` - Tang-Xiao variant of SROCK2

These methods are designed for stiff SDEs where standard explicit methods would require prohibitively small time steps for stability. The Chebyshev polynomial basis enables an extended stability region along the negative real axis.

## Example

```julia
using StochasticDiffEqROCK
f(u, p, t) = -100.0 * u  # stiff drift
g(u, p, t) = 0.1 * u
u0 = 1.0
tspan = (0.0, 1.0)
prob = SDEProblem(f, g, u0, tspan)
sol = solve(prob, SROCK1(), dt = 0.01)
```

## Available Solvers

For the full list of available SDE solvers, refer to the [DifferentialEquations.jl SDE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/) documentation.
