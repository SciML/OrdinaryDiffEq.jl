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
@test_broken solve(prob,Rosenbrock23()).retcode == :Success
solve(prob,Rosenbrock23(autodiff=false))

prob_nojac = ODEProblem(f,u0,tspan)
@test_broken solve(prob_nojac,Rosenbrock23()).retcode == :Success
@test_broken solve(prob_nojac,Rosenbrock23(autodiff=false)).retcode == :Success

# Test auto-offload
_A = -rand(3,3)
function f2(du,u,p,t)
    mul!(du,_A,u)
end
function jac2(J,u,p,t)
    J .= _A
end
ff2 = ODEFunction(f2,jac=jac2)
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob_num = ODEProblem(ff2,u0,tspan)
sol = solve(prob_num,Rosenbrock23(linsolve=LinSolveGPUFactorize()))
