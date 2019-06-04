using StaticArrays, Test
using OrdinaryDiffEq

u0 = fill(zero(MVector{2,Float64}), 2)
u0[1] = ones(MVector{2,Float64}) .+ 1
f = (du,u,p,t) -> du .= u
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1e-2)
sol = solve(ode, Tsit5())

u0 = fill(zero(SVector{2,Float64}), 2) .+ 1
u0[1] = ones(SVector{2,Float64}) .+ 1
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1e-2)
sol = solve(ode, Tsit5())

sol = solve(ode, SSPRK22(), dt=1e-2)


u0 = ones(MVector{2,Float64})
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1e-2)
sol = solve(ode, Tsit5(), dt=1e-2)


u0 = ones(SVector{2,Float64})
f = (u,p,t) -> u
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1e-2)
integrator = init(ode, ImplicitEuler())
@test OrdinaryDiffEq.get_W(integrator.cache.nlsolver) isa StaticArrays.LU
sol = solve(ode, ImplicitEuler())
sol = solve(ode, ImplicitEuler(nlsolve=NLAnderson()))
sol = solve(ode, Tsit5(), dt=1e-2)
