# OrdinaryDiffEq.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/OrdinaryDiffEq.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/OrdinaryDiffEq.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/crn27g5aj1r567m5?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/ordinarydiffeq-jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiffEq/OrdinaryDiffEq.jl/badge.svg)](https://coveralls.io/github/JuliaDiffEq/OrdinaryDiffEq.jl)
[![codecov.io](http://codecov.io/github/ChrisRackauckas/OrdinaryDiffEq.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/OrdinaryDiffEq.jl?branch=master)
[![OrdinaryDiffEq](http://pkg.julialang.org/badges/OrdinaryDiffEq_0.6.svg)](http://pkg.julialang.org/?pkg=OrdinaryDiffEq)

OrdinaryDiffEq.jl is a component package in the DifferentialEquations ecosystem. It holds the
ordinary differential equation solvers and utilities. While completely independent
and usable on its own, users interested in using this
functionality should check out [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

## API

OrdinaryDiffEq.jl is part of the JuliaDiffEq common interface, but can be used independently of DifferentialEquations.jl. The only requirement is that the user passes an OrdinaryDiffEq.jl algorithm to `solve`. For example, we can solve the [ODE tutorial from the docs](http://docs.juliadiffeq.org/latest/tutorials/ode_example.html) using the `Tsit5()` algorithm:

```julia
using OrdinaryDiffEq
f(t,u) = 1.01*u
u0=1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
using Plots
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
plot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")
```

That example uses the out-of-place syntax `f(t,u)`, while the inplace syntax (more efficient for systems of equations) is shown in the Lorenz example:

```julia
using OrdinaryDiffEq
function lorenz(t,u,du)
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz,u0,tspan)
sol = solve(prob,Tsit5())
using Plots; plot(sol,vars=(1,2,3))
```

For "refined ODEs" like dynamical equations and `SecondOrderODEProblem`s, refer to the [DiffEqDocs](http://docs.juliadiffeq.org/latest/types/ode_types.html). For example, in [DiffEqTutorials.jl](https://github.com/JuliaDiffEq/DiffEqTutorials.jl) we show how to solve equations of motion using symplectic methods:

```julia
function HH_acceleration(t, u, v, dv)
    x,y  = u
    dx,dy = dv
    dv[1] = -x - 2x*y
    dv[2] = y^2 - y -x^2
end
initial_positions = [0.0,0.1]
initial_velocities = [0.5,0.0]
prob = SecondOrderODEProblem(HH_acceleration,initial_positions,initial_velocities, tspan)
sol2 = solve(prob, KahanLi8(), dt=1/10);
```

## Available Solvers

For the list of available solvers, please refer to the [DifferentialEquations.jl ODE Solvers page](http://docs.juliadiffeq.org/latest/solvers/ode_solve.html#OrdinaryDiffEq.jl-1) and the [Refined ODE Solvers page](http://docs.juliadiffeq.org/latest/solvers/refined_ode_solve.html).
