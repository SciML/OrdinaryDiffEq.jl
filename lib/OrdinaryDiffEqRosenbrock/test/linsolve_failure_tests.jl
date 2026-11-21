using OrdinaryDiffEqRosenbrock, LinearSolve, Test

# ∂w'/∂c = -Inf at c = 0, so W is non-finite and LU reports a failed factorization
function f!(du, u, p, t)
    d, c, w = u
    du[1] = -d
    du[2] = d - c
    du[3] = -(c < 0 ? NaN : c^0.1) - w
    return nothing
end
prob = ODEProblem(f!, [1.0, 0.0, 0.0], (0.0, 10.0))

# W = I/(γ dt) - λ I is finite but exactly singular at dt = 0.01
γ = OrdinaryDiffEqRosenbrock.Rodas5PTableau(Float64, Float64).gamma
λ = 1 / (γ * 0.01)
singular_prob = ODEProblem((du, u, p, t) -> (du .= λ .* u), [1.0], (0.0, 0.05))
singular_prob_oop = ODEProblem((u, p, t) -> λ .* u, [1.0], (0.0, 0.05))

@testset "Failed linear solve rejects the step: $(nameof(typeof(alg)))" for alg in (
        Rodas5P(; linsolve = LUFactorization()),
        Rodas4(; linsolve = LUFactorization()),
        Rosenbrock23(; linsolve = LUFactorization()),
        Rosenbrock32(; linsolve = LUFactorization()),
    )
    sol = solve(prob, alg)
    @test !SciMLBase.successful_retcode(sol)
    @test sol.t[end] == 0.0

    sol = solve(prob, alg; adaptive = false, dt = 1.0e-3)
    @test sol.retcode == ReturnCode.ConvergenceFailure
    @test sol.t[end] == 0.0
end

@testset "Failed linear solve is not accepted under force_dtmin" begin
    sol = solve(prob, Rodas5P(; linsolve = LUFactorization()); dt = 1.0e-3, dtmin = 1.0e-3, force_dtmin = true)
    @test !SciMLBase.successful_retcode(sol)
    @test sol.t[end] == 0.0
end

@testset "Singular W" begin
    sol = solve(singular_prob, Rodas5P(; linsolve = LUFactorization()); dt = 0.01)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end][1] ≈ exp(λ * 0.05) rtol = 1.0e-2

    sol = solve(singular_prob, Rodas5P(; linsolve = LUFactorization()); dt = 0.01, adaptive = false)
    @test sol.retcode == ReturnCode.ConvergenceFailure

    sol = solve(singular_prob_oop, Rodas5P(); dt = 0.01)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end][1] ≈ exp(λ * 0.05) rtol = 1.0e-2

    sol = solve(singular_prob_oop, Rodas5P(); dt = 0.01, adaptive = false)
    @test sol.retcode == ReturnCode.ConvergenceFailure
end
