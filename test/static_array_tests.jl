using StaticArrays
using DiffEqBase, OrdinaryDiffEq

u0 = zeros(MVector{2,Float64}, 2)
u0[1] = ones(MVector{2,Float64}) + 1
f = (t,u,du) -> du .= u
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1.e-2)
#sol = solve(ode, Tsit5())

u0 = zeros(SVector{2,Float64}, 2) + 1
u0[1] = ones(SVector{2,Float64}) + 1
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1.e-2)
#sol = solve(ode, Tsit5())

sol = solve(ode, SSPRK22(), dt=1.e-2)


u0 = ones(MVector{2,Float64})
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1.e-2)
sol = solve(ode, Tsit5(), dt=1.e-2)


u0 = ones(SVector{2,Float64})
f = (t,u) -> u
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1.e-2)
sol = solve(ode, Tsit5(), dt=1.e-2)
