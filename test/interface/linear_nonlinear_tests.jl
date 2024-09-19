using OrdinaryDiffEq, Test, Random, LinearAlgebra, LinearSolve
Random.seed!(123)

const A = 0.01 * rand(3, 3)
rn = (du, u, p, t) -> begin
    du .= A * u
end
u0 = rand(3)
prob = ODEProblem(rn, u0, (0, 50.0))

function precsl(W, p)
    F = lu(convert(AbstractMatrix, W), check = false)
    return F, IdentityOperator(size(W, 1))
end

function precsr(W, p)
    F = lu(convert(AbstractMatrix, W), check = false)
    IdentityOperator(size(W, 1)), F
end

function precslr(W, p)
    F = lu(convert(AbstractMatrix, W), check = false)
    F, F
end

@testset "precs" begin
    @testset "$linsolve" for linsolve in (KrylovJL_GMRES(),
                              KrylovJL_GMRES(precs = precsl),
                              KrylovJL_GMRES(precs = precsr),
                              KrylovJL_GMRES(precs = precslr))
        sol = @test_nowarn solve(prob, TRBDF2(;linsolve,
            smooth_est = false, concrete_jac = true), maxiters=20)
        sol = @test_nowarn solve(prob, Rodas5P(autodiff = false;
            linsolve, concrete_jac = true), maxiters=30)
    end
    @testset "$solver" for solver in (Rosenbrock23, FBDF, QNDF)
        sol = @test_nowarn solve(prob, solver(
            linsolve=KrylovJL_GMRES(precs = precslr), concrete_jac = true),
            maxiters=30)
    end
end

sol = @test_nowarn solve(prob, TRBDF2(autodiff = false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob,
    TRBDF2(autodiff = false, linsolve = KrylovJL_GMRES()));
@test length(sol.t) < 20
solref = @test_nowarn solve(prob,
    TRBDF2(autodiff = false, linsolve = KrylovJL_GMRES(),
        smooth_est = false));
@test length(sol.t) < 20
sol = @test_nowarn solve(prob,
    QNDF(autodiff = false, linsolve = KrylovJL_GMRES(),
        concrete_jac = true));
@test length(sol.t) < 25
