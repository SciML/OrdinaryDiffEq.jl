# StochasticDiffEq.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

StochasticDiffEq.jl is a component package in the DifferentialEquations ecosystem. It holds the stochastic differential equation solvers and utilities. While completely independent and usable on its own, users interested in using this functionality should check out [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

This is an umbrella package that re-exports all SDE solver subpackages from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo.

## Installation

```julia
import Pkg;
Pkg.add("StochasticDiffEq");
```

## API

StochasticDiffEq.jl is part of the SciML common interface, but can be used independently of DifferentialEquations.jl. The only requirement is that the user passes a StochasticDiffEq.jl algorithm to `solve`. For example:

```julia
using StochasticDiffEq
f(u, p, t) = 1.01 * u
g(u, p, t) = 0.87 * u
u0 = 0.5
tspan = (0.0, 1.0)
prob = SDEProblem(f, g, u0, tspan)
sol = solve(prob, SOSRI())
```

The in-place syntax (more efficient for systems of equations):

```julia
using StochasticDiffEq
function f!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
function g!(du, u, p, t)
    du[1] = 0.3 * u[1]
    du[2] = 0.3 * u[2]
    du[3] = 0.3 * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = SDEProblem(f!, g!, u0, tspan)
sol = solve(prob, SOSRI())
```

## Subpackages

For more targeted dependencies, you can use individual solver subpackages:

| Package | Solvers |
|---------|---------|
| StochasticDiffEqLowOrder | EM, EulerHeun, LambaEM, LambaEulerHeun, SimplifiedEM, SplitEM, RKMil, RKMilCommute, PCEuler |
| StochasticDiffEqHighOrder | SRI, SRIW1, SRIW2, SOSRI, SOSRI2, SRA, SRA1, SRA2, SRA3, SOSRA, SOSRA2 |
| StochasticDiffEqMilstein | RKMilGeneral, WangLi3SMil_A-F |
| StochasticDiffEqROCK | SROCK1, SROCK2, KomBurSROCK2, SROCKC2, SROCKEM, SKSROCK, TangXiaoSROCK2 |
| StochasticDiffEqImplicit | ImplicitEM, ImplicitEulerHeun, ImplicitRKMil, STrapezoid, SImplicitMidpoint, ISSEM, ISSEulerHeun, SKenCarp |
| StochasticDiffEqWeak | DRI1, RI1, RI3, RI5, RI6, RDI1WM-RDI4WM, W2Ito1, RS1, RS2, PL1WM, NON, COM, IRI1, and more |
| StochasticDiffEqIIF | IIF1M, IIF2M, IIF1Mil |
| StochasticDiffEqLeaping | TauLeaping, CaoTauLeaping, ImplicitTauLeaping, ThetaTrapezoidalTauLeaping |
| StochasticDiffEqRODE | RandomEM, RandomHeun, RandomTamedEM, BAOAB |

## Available Solvers

For the full list of available solvers, refer to the [DifferentialEquations.jl SDE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/), [RODE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/rode_solve/), and [Jump Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/jump_solve/) pages.
