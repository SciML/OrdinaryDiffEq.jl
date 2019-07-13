using OrdinaryDiffEq, CuArrays, LinearAlgebra, Test
function f(du,u,p,t)
    mul!(du,A,u)
end
function jac(J,u,p,t)
    J .= A
end
ff = ODEFunction(f,jac=jac)
A = cu(-rand(3,3))
u0 = cu([1.0;0.0;0.0])
tspan = (0.0,100.0)
prob = ODEProblem(ff,u0,tspan)

CuArrays.allowscalar(false)
sol = solve(prob,Tsit5())
sol = solve(prob,Rosenbrock23())

prob_nojac = ODEProblem(f,u0,tspan)
@test_broken sol = solve(prob_nojac,Rosenbrock23())

# Test auto-offload
A = -rand(3,3))
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob_nojac2 = ODEProblem(f,u0,tspan)
sol = solve(prob,Rosenbrock23(linsolve=LinSolveGPUFactorize()))
