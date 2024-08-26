# OrdinaryDiffEq.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/OrdinaryDiffEq/stable/)

[![codecov](https://codecov.io/gh/SciML/OrdinaryDiffEq.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/SciML/OrdinaryDiffEq.jl)
[![Build Status](https://github.com/SciML/OrdinaryDiffEq.jl/workflows/CI/badge.svg)](https://github.com/SciML/OrdinaryDiffEq.jl/actions?query=workflow%3ACI)
[![Build status](https://badge.buildkite.com/5f39777d009ce94ef1dcf2a4881c68b9fbcaf6f69f1d8b8df2.svg?branch=master)](https://buildkite.com/julialang/ordinarydiffeq-dot-jl)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

OrdinaryDiffEq.jl is a component package in the DifferentialEquations ecosystem. It holds the
ordinary differential equation solvers and utilities. While completely independent
and usable on its own, users interested in using this
functionality should check out [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

## Installation

Assuming that you already have Julia correctly installed, it suffices to import
OrdinaryDiffEq.jl in the standard way:

```julia
import Pkg;
Pkg.add("OrdinaryDiffEq");
```

## API

OrdinaryDiffEq.jl is part of the SciML common interface, but can be used independently of DifferentialEquations.jl. The only requirement is that the user passes an OrdinaryDiffEq.jl algorithm to `solve`. For example, we can solve the [ODE tutorial from the docs](https://docs.sciml.ai/DiffEqDocs/stable/getting_started/#ode_example) using the `Tsit5()` algorithm:

```julia
using OrdinaryDiffEq
f(u, p, t) = 1.01 * u
u0 = 1 / 2
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
using Plots
plot(sol, linewidth = 5, title = "Solution to the linear ODE with a thick line",
    xaxis = "Time (t)", yaxis = "u(t) (in μm)", label = "My Thick Line!") # legend=false
plot!(sol.t, t -> 0.5 * exp(1.01t), lw = 3, ls = :dash, label = "True Solution!")
```

That example uses the out-of-place syntax `f(u,p,t)`, while the inplace syntax (more efficient for systems of equations) is shown in the Lorenz example:

```julia
using OrdinaryDiffEq
function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan)
sol = solve(prob, Tsit5())
using Plots;
plot(sol, idxs = (1, 2, 3))
```

Very fast static array versions can be specifically compiled to the size of your model. For example:

```julia
using OrdinaryDiffEq, StaticArrays
function lorenz(u, p, t)
    SA[10.0(u[2] - u[1]), u[1] * (28.0 - u[3]) - u[2], u[1] * u[2] - (8 / 3) * u[3]]
end
u0 = SA[1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)
sol = solve(prob, Tsit5())
```

For "refined ODEs", like dynamical equations and `SecondOrderODEProblem`s, refer to the [DiffEqDocs](https://diffeq.sciml.ai/dev/types/ode_types/). For example, in [DiffEqTutorials.jl](https://github.com/SciML/DiffEqTutorials.jl) we show how to solve equations of motion using symplectic methods:

```julia
function HH_acceleration!(dv, v, u, p, t)
    x, y = u
    dx, dy = dv
    dv[1] = -x - 2x * y
    dv[2] = y^2 - y - x^2
end
initial_positions = [0.0, 0.1]
initial_velocities = [0.5, 0.0]
prob = SecondOrderODEProblem(HH_acceleration!, initial_velocities, initial_positions, tspan)
sol2 = solve(prob, KahanLi8(), dt = 1 / 10);
```

Other refined forms are IMEX and semi-linear ODEs (for exponential integrators).

## Available Solvers

For the list of available solvers, please refer to the [DifferentialEquations.jl ODE Solvers](https://diffeq.sciml.ai/dev/solvers/ode_solve/), [Dynamical ODE Solvers](http://diffeq.sciml.ai/dev/solvers/dynamical_solve/), and the [Split ODE Solvers](http://diffeq.sciml.ai/dev/solvers/split_ode_solve/) pages.
