# OrdinaryDiffEq.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/OrdinaryDiffEq.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/OrdinaryDiffEq.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/crn27g5aj1r567m5?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/ordinarydiffeq-jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiffEq/OrdinaryDiffEq.jl/badge.svg)](https://coveralls.io/github/JuliaDiffEq/OrdinaryDiffEq.jl)
[![codecov.io](http://codecov.io/github/ChrisRackauckas/OrdinaryDiffEq.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/OrdinaryDiffEq.jl?branch=master)
[![OrdinaryDiffEq](http://pkg.julialang.org/badges/OrdinaryDiffEq_0.5.svg)](http://pkg.julialang.org/?pkg=OrdinaryDiffEq)
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

## Available Solvers

For the list of available solvers, please refer to the [DifferentialEquations.jl ODE Solvers page](http://docs.juliadiffeq.org/latest/solvers/ode_solve.html#OrdinaryDiffEq.jl-1) and the [Refined ODE Solvers page](http://docs.juliadiffeq.org/latest/solvers/refined_ode_solve.html).
