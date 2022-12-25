using OrdinaryDiffEq, Test

const a = Float64[1.0]

function lorenz(u,p,t)
    du1 = 10.0(u[2]-u[1])
    a[1] = u[2]
    du2 = u[1]*(28.0-u[3]) - u[2]
    du3 = u[1]*u[2] - (8/3)*u[3]
    [du1,du2,du3]
   end
u0 = [1.0;0.0;0.0]
tspan = (0.0,1.0)
prob = ODEProblem(lorenz,u0,tspan)
@test_throws OrdinaryDiffEq.FirstAutodiffJacError solve(prob,Rosenbrock23())

function lorenz(u,p,t)
    du1 = 10.0(u[2]-u[1])
    a[1] = t
    du2 = u[1]*(28.0-u[3]) - u[2]
    du3 = u[1]*u[2] - (8/3)*u[3]
    [du1,du2,du3]
end
@test_throws OrdinaryDiffEq.FirstAutodiffTgradError solve(prob,Rosenbrock23())

function lorenz!(du,u,p,t)
    du[1] = 10.0(u[2]-u[1])
    a[1] = u[2]
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
   end
u0 = [1.0;0.0;0.0]
tspan = (0.0,1.0)
prob = ODEProblem(lorenz!,u0,tspan)
@test_throws OrdinaryDiffEq.FirstAutodiffJacError solve(prob,Rosenbrock23())

function lorenz!(du,u,p,t)
    du[1] = 10.0(u[2]-u[1])
    a[1] = t
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end
@test_throws OrdinaryDiffEq.FirstAutodiffTgradError solve(prob,Rosenbrock23())