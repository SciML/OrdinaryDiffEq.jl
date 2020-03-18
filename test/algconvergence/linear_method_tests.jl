using OrdinaryDiffEq, Test, DiffEqDevTools
using LinearAlgebra, Random


# Linear exponential solvers
A = DiffEqArrayOperator([2.0 -1.0; -1.0 2.0])
u0 = ones(2)
prob = ODEProblem(A,u0,(0.0,1.0))
sol1 = solve(prob, LinearExponential(krylov=:off))(1.0)
sol2 = solve(prob, LinearExponential(krylov=:simple))(1.0)
sol3 = solve(prob, LinearExponential(krylov=:adaptive))(1.0)
sol_analytic = exp(1.0 * Matrix(A)) * u0

@test isapprox(sol1, sol_analytic, rtol=1e-10)
@test isapprox(sol2, sol_analytic, rtol=1e-10)
@test isapprox(sol3, sol_analytic, rtol=1e-10)

# u' = A(t)u solvers

function update_func(A,u,p,t)
    A[1,1] = cos(t)
    A[2,1] = sin(t)
    A[1,2] = cos(t)*sin(t)
    A[2,2] = sin(t)^2
end
A = DiffEqArrayOperator(ones(2,2),update_func=update_func)
prob = ODEProblem(A, ones(2), (0.0, 5.))
dts = 1 ./2 .^(10:-1:1)
sol  = solve(prob,OrdinaryDiffEq.MagnusMidpoint(),dt=1/4)

dts = 1 ./2 .^(10:-1:1)
test_setup = Dict(:alg=>Vern9(),:reltol=>1e-14,:abstol=>1e-14)
sim = analyticless_test_convergence(dts,prob,MagnusMidpoint(),test_setup)
@test sim.𝒪est[:l2] ≈ 2 atol=0.2
sim = analyticless_test_convergence(dts,prob,MagnusMidpoint(krylov=true),test_setup)
@test sim.𝒪est[:l2] ≈ 2 atol=0.2

A = DiffEqArrayOperator(ones(2,2),update_func=update_func)
prob = ODEProblem(A, ones(2), (0.0, 5.))
dts = 1 ./2 .^(10:-1:1)
sol  = solve(prob,OrdinaryDiffEq.MagnusLeapfrog(),dt=1/4)

dts = 1 ./2 .^(10:-1:1)
test_setup = Dict(:alg=>Vern9(),:reltol=>1e-14,:abstol=>1e-14)
sim = analyticless_test_convergence(dts,prob,MagnusLeapfrog(),test_setup)
@test sim.𝒪est[:l2] ≈ 2 atol=0.2
sim = analyticless_test_convergence(dts,prob,MagnusLeapfrog(krylov=true),test_setup)
@test sim.𝒪est[:l2] ≈ 2 atol=0.2

function update_func(A,u,p,t)
    A[1,1] = cos(t)*u[1]
    A[2,1] = sin(t)*u[2]
    A[1,2] = -cos(t)*sin(t)
    A[2,2] = -(sin(t)^2) *u[1] * u[2]
end
A = DiffEqArrayOperator(ones(2,2),update_func=update_func)
prob = ODEProblem(A, ones(2), (0.5, 5.))
dts = 1 ./2 .^(10:-1:1)
sol  = solve(prob,OrdinaryDiffEq.LieEuler(),dt=1/4)

dts = 1 ./2 .^(10:-1:1)
test_setup = Dict(:alg=>Vern9(),:reltol=>1e-14,:abstol=>1e-14)
sim = analyticless_test_convergence(dts,prob,LieEuler(),test_setup)
@test sim.𝒪est[:l2] ≈ 1 atol=0.2
sim = analyticless_test_convergence(dts,prob,LieEuler(krylov=true),test_setup)
@test sim.𝒪est[:l2] ≈ 1 atol=0.2
