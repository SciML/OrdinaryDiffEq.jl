# StochasticDiffEqLeaping.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/jump_solve/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

StochasticDiffEqLeaping.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo providing tau-leaping methods for jump process simulation. While completely independent and usable on its own, users wanting all SDE solvers should use [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

## Solvers

- `TauLeaping` - Explicit tau-leaping
- `CaoTauLeaping` - Cao's adaptive tau-leaping with step size selection
- `ImplicitTauLeaping` - Implicit tau-leaping for stiff systems
- `ThetaTrapezoidalTauLeaping` - Theta-method trapezoidal tau-leaping

These methods approximate the stochastic simulation algorithm (SSA) by allowing multiple reactions to fire in each time step, providing significant speedups for systems with many reactions.

## Example

```julia
using StochasticDiffEqLeaping, JumpProcesses
# Define a simple birth-death process via JumpProblem
# then solve with tau-leaping:
sol = solve(jump_prob, TauLeaping(), dt = 0.01)
```

## Available Solvers

For the full list of available jump solvers, refer to the [DifferentialEquations.jl Jump Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/jump_solve/) documentation.
