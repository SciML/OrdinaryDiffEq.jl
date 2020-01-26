using OrdinaryDiffEq, DiffEqDevTools
using Test, Random
Random.seed!(100)

dts = 1 .//2 .^(9:-1:5)
testTol = 0.2

f_dae_linear = (res, du, u, p, t) -> (@. res = du - u)
f_dae_linear_analytic = (du0, u0, p, t) -> @. u0*exp(t)
prob_dae_linear_iip = DAEProblem(
				  DAEFunction(f_dae_linear;analytic=f_dae_linear_analytic),
				  [1.0,1.0], [1.0,1.0], (0.0, 1.0)
				  )

@testset "DAE Solver Convergence Tests (in-place)" begin
	prob = prob_dae_linear_iip

	sim1 = test_convergence(dts,prob,DImplicitEuler())
	@test sim1.𝒪est[:final] ≈ 1 atol=testTol
end

f_dae_linear = (du, u, p, t) -> (@. du - u)
f_dae_linear_analytic = (du0, u0, p, t) -> @. u0*exp(t)
prob_dae_linear_oop = DAEProblem(
				  DAEFunction(f_dae_linear;analytic=f_dae_linear_analytic),
				  1.0, 1.0, (0.0, 1.0)
				  )

@testset "DAE Solver Convergence Tests (out-of-place)" begin
	prob = prob_dae_linear_oop

	sim2 = test_convergence(dts,prob,DImplicitEuler())
	@test sim2.𝒪est[:final] ≈ 1 atol=testTol
end