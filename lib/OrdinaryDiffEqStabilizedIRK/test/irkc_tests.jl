
@testset "Power Iteration of Runge-Kutta-Chebyshev Tests" begin
    Random.seed!(123)
    eigen_est = (integrator) -> integrator.eigen_est = 1.5e2
    for iip in [true, false], Alg in [IRKC]
        alg = Alg()
        println(typeof(alg))
        A = randn(20, 20)
        B = randn(20, 20)
        test_f1 = !iip ? (u, p, t) -> A * u : (du, u, p, t) -> mul!(du, A, u)
        test_f2 = !iip ? (u, p, t) -> B * u : (du, u, p, t) -> mul!(du, B, u)
        ff_split = SplitFunction{iip}(test_f1, test_f2)
        prob = SplitODEProblem{iip}(ff_split, randn(20, 1), (0.0, 1.0))
        integrator = init(prob, alg)
        eigm = maximum(abs.(eigvals(A)))
        maxeig!(integrator, integrator.cache)
        eigest = integrator.eigen_est
        @test eigest≈eigm rtol=0.1eigm

        A = A - 1e2I
        test_f1 = !iip ? (u, p, t) -> A * u : (du, u, p, t) -> mul!(du, A, u)
        prob = SplitODEProblem{iip}(SplitFunction{iip}(test_f1, test_f2), ones(20),
            (0.0, 1.0))
        @test_nowarn solve(prob, alg)
    end
end