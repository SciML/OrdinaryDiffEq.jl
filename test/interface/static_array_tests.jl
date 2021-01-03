using StaticArrays, Test
using OrdinaryDiffEq
using RecursiveArrayTools

u0 = [fill(2, MVector{2,Float64}), ones(MVector{2,Float64})]
f = (du,u,p,t) -> du .= u
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, Tsit5())
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)

u0 = [fill(2, SVector{2,Float64}), ones(SVector{2,Float64})]
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, Tsit5())
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, SSPRK22(), dt=1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)

u0 = ones(MVector{2,Float64})
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, Tsit5(), dt=1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)

u0 = ones(SVector{2,Float64})
f = (u,p,t) -> u
ode = ODEProblem(f, u0, (0.,1.))
sol = solve(ode, Euler(), dt=1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
integrator = init(ode, ImplicitEuler())
@test OrdinaryDiffEq.get_W(integrator.cache.nlsolver) isa StaticArrays.LU
sol = solve(ode, ImplicitEuler())
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, ImplicitEuler(nlsolve=NLAnderson()))
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)
sol = solve(ode, Tsit5(), dt=1e-2)
@test !any(iszero.(sol(1.0))) && !any(sol(1.0) .== u0)

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

# Check that ArrayPartitions of static vectors work
#https://github.com/SciML/OrdinaryDiffEq.jl/issues/1308
function lorenz_static(u::ArrayPartition,p,t)
 dx = 10.0*(u[2]-u[1])
 dy = u[1]*(28.0-u[3]) - u[2]
 dz = u[1]*u[2] - (8/3)*u[3]
 du1 = @SVector [dx,dy]
 du2 = @SVector [dz]
 ArrayPartition(du1,du2)
end

u01 = @SVector [1.0,0.0]
u02 = @SVector [0.0]
u0ap = ArrayPartition(u01,u02)
probap = ODEProblem(lorenz_static,u0ap,tspan)

sol = solve(prob,dt=1e-2,Heun())
solap = solve(probap,dt=1e-2,Heun())
@test sol(30) == solap(30)

sol = solve(prob,dt=1e-2,Tsit5())
solap = solve(probap,dt=1e-2,Tsit5())
@test sol(30) == solap(30)
