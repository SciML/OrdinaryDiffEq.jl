# OrdinaryDiffEqVerner.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/nonstiff/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

OrdinaryDiffEqVerner.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo. Verner's efficient explicit Runge-Kutta methods of orders 6 through 9.
While completely independent and usable on its own, users wanting the full ODE solver suite should use [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).

## Solvers

- `Vern6`
- `Vern7`
- `Vern8`
- `Vern9`
- `AutoVern6`
- `AutoVern7`
- `AutoVern8`
- `AutoVern9`

For a full description of the solvers and their options, see the [ODE solver documentation](https://docs.sciml.ai/DiffEqDocs/stable/).
