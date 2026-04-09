# OrdinaryDiffEqSDIRK.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/stiff/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

OrdinaryDiffEqSDIRK.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo. Singly Diagonally Implicit Runge-Kutta (SDIRK) and ESDIRK methods for stiff problems.
While completely independent and usable on its own, users wanting the full ODE solver suite should use [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).

## Solvers

- `ImplicitEuler`
- `Trapezoid`
- `TRBDF2`
- `Kvaerno3`
- `KenCarp3`
- `Kvaerno4`
- `KenCarp4`
- `KenCarp5`
- `ESDIRK54I8L2SA`

For a full description of the solvers and their options, see the [ODE solver documentation](https://docs.sciml.ai/DiffEqDocs/stable/).
