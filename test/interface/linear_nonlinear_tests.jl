using OrdinaryDiffEq, Test, Random
Random.seed!(123)

using OrdinaryDiffEq, DiffEqOperators, LinearAlgebra
A = 0.01*rand(3, 3)
rn = (du, u, p, t) -> begin
    mul!(du, A, u)
end
u0 = rand(3)
prob = ODEProblem(ODEFunction(rn, jac_prototype=JacVecOperator{Float64}(rn, u0; autodiff=false)), u0, (0, 10.))
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=LinSolveGMRES()));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=LinSolveGMRES(), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=LinSolveGMRES(Pl=lu(A)), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=LinSolveGMRES(Pr=lu(A)), smooth_est=false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob, TRBDF2(autodiff=false, linsolve=LinSolveGMRES(Pl=lu(A), Pr=lu(A)), smooth_est=false));
@test length(sol.t) < 20

sol = @test_nowarn solve(prob, QNDF(autodiff=false, linsolve=LinSolveGMRES()));
@test length(sol.t) < 20
