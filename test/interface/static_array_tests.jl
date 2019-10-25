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

#https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/373
function lorenz_static(u,p,t)
 dx = 10.0*(u[2]-u[1])
 dy = u[1]*(28.0-u[3]) - u[2]
 dz = u[1]*u[2] - (8/3)*u[3]
 @SVector [dx,dy,dz]
end

u0 = @SVector [1.0,0.0,0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz_static,u0,tspan)
solve(prob,dt=0.1,Rosenbrock23(autodiff=false))
