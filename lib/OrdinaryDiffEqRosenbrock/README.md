# OrdinaryDiffEqRosenbrock.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/stiff/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

OrdinaryDiffEqRosenbrock.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo. Rosenbrock and Rosenbrock-Wanner (ROW) methods for stiff ODEs.
While completely independent and usable on its own, users wanting the full ODE solver suite should use [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).

## Solvers

- `Rosenbrock23`
- `Rosenbrock32`
- `RosShamp4`
- `Rodas3`
- `Rodas4`
- `Rodas4P`
- `Rodas5`
- `Rodas5P`
- `Rodas5Pe`

For a full description of the solvers and their options, see the [ODE solver documentation](https://docs.sciml.ai/DiffEqDocs/stable/).
