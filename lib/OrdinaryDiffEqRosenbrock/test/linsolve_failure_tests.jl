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

# W = λ I - I/(γ dt) is finite but exactly singular at dt = 0.01
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
    for alg in (
                Rodas5P(; linsolve = LUFactorization()),
                Rodas4(; linsolve = LUFactorization()),
                Rosenbrock23(; linsolve = LUFactorization()),
                Rosenbrock32(; linsolve = LUFactorization()),
            ), dt in (1.0e-3, 4.0e-3)
        sol = solve(prob, alg; dt, dtmin = 1.0e-3, force_dtmin = true, maxiters = 25)
        @test sol.retcode == ReturnCode.ConvergenceFailure
        @test sol.t[end] == 0.0
        @test sol.u[end] == prob.u0
        @test sol.stats.naccept == 0
        @test sol.stats.nw <= 3
    end

    # Exercise each out-of-place W-factorization guard with an exactly singular
    # finite matrix, including integration in the negative time direction.
    for (alg, gamma) in (
                (Rodas5P(), γ),
                (Rosenbrock23(), OrdinaryDiffEqRosenbrock.Rosenbrock23Tableau(Float64).d),
                (Rosenbrock32(), OrdinaryDiffEqRosenbrock.Rosenbrock32Tableau(Float64).d),
            ), direction in (-1, 1)
        dt = direction * 0.01
        lambda = 1 / (gamma * dt)
        singular = ODEProblem((u, p, t) -> lambda .* u, [1.0], (0.0, direction * 0.05))
        sol = solve(singular, alg; dt, dtmin = 0.01, force_dtmin = true, maxiters = 25)
        @test sol.retcode == ReturnCode.ConvergenceFailure
        @test sol.t[end] == 0.0
        @test sol.u[end] == singular.u0
        @test sol.stats.naccept == 0
        @test sol.stats.nw <= 1
    end
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

# Inject a failure after successful stages, without depending on a particular
# factorization backend or a matrix becoming singular at a later stage.
struct StageFailureLinearSolve <: LinearSolve.AbstractSolveFunction
    calls::Base.RefValue{Int}
    fail_at::Int
    retcode::SciMLBase.ReturnCode.T
end
LinearSolve.needs_concrete_A(::StageFailureLinearSolve) = true
function SciMLBase.solve!(cache::LinearSolve.LinearCache, alg::StageFailureLinearSolve; kwargs...)
    alg.calls[] += 1
    u = cache.A \ cache.b
    retcode = alg.calls[] == alg.fail_at ? alg.retcode : ReturnCode.Success
    return SciMLBase.build_linear_solution(alg, u, nothing, cache; retcode)
end

decay_prob = ODEProblem((du, u, p, t) -> (du .= -u), [1.0], (0.0, 0.3))

@testset "Failed stages preserve the accepted state: $ctor" for ctor in (
        Rodas5P, Rodas4, Rosenbrock23, Rosenbrock32,
    )
    calls = Ref(0)
    counting_solver = StageFailureLinearSolve(calls, typemax(Int), ReturnCode.Failure)
    reference = solve(remake(decay_prob; tspan = (0.0, 0.1)), ctor(; linsolve = counting_solver); dt = 0.1, adaptive = false)
    stages = calls[]

    for failed_stage in 2:stages, save_start in (false, true)
        linsolve = StageFailureLinearSolve(Ref(0), stages + failed_stage, ReturnCode.Failure)
        integrator = init(
            decay_prob, ctor(; linsolve); dt = 0.1, adaptive = false,
            save_start, save_everystep = false, save_end = true
        )
        sol = solve!(integrator)
        @test sol.retcode == ReturnCode.ConvergenceFailure
        @test sol.t[end] == 0.1
        @test sol.u[end] ≈ reference.u[end]
        @test integrator.u == integrator.uprev

        linsolve = StageFailureLinearSolve(Ref(0), failed_stage, ReturnCode.Failure)
        sol = solve(
            decay_prob, ctor(; linsolve); dt = 0.1, adaptive = false,
            save_start, save_everystep = false, save_end = true
        )
        @test sol.retcode == ReturnCode.ConvergenceFailure
        @test sol.t[end] == 0.0
        @test sol.u[end] == decay_prob.u0
    end
end

@testset "Transient failure recovers and Default is not rejected: $ctor" for ctor in (
        Rodas5P, Rodas4, Rosenbrock23, Rosenbrock32,
    )
    for force_dtmin in (false, true)
        linsolve = StageFailureLinearSolve(Ref(0), 2, ReturnCode.Failure)
        sol = solve(decay_prob, ctor(; linsolve); dt = 0.1, dtmin = 0.025, force_dtmin, abstol = 1.0e3, reltol = 1.0e3)
        reference = solve(decay_prob, ctor(); dt = 0.05, dtmin = 0.025, force_dtmin, abstol = 1.0e3, reltol = 1.0e3)
        @test SciMLBase.successful_retcode(sol)
        @test sol.t == reference.t
        @test sol.u ≈ reference.u
    end

    linsolve = StageFailureLinearSolve(Ref(0), 2, ReturnCode.Default)
    sol = solve(decay_prob, ctor(; linsolve); dt = 0.1, adaptive = false)
    reference = solve(decay_prob, ctor(); dt = 0.1, adaptive = false)
    @test SciMLBase.successful_retcode(sol)
    @test sol.u ≈ reference.u
end

@testset "Failure at dtmin preserves a previously accepted state: $ctor" for ctor in (
        Rodas5P, Rodas4, Rosenbrock23, Rosenbrock32,
    )
    for direction in (-1, 1)
        dt = direction * 0.1
        calls = Ref(0)
        linsolve = StageFailureLinearSolve(calls, typemax(Int), ReturnCode.Failure)
        reference = solve(remake(decay_prob; tspan = (0.0, dt)), ctor(; linsolve); dt, dtmin = 0.1, dtmax = 0.1, force_dtmin = true, abstol = 1.0e3, reltol = 1.0e3)
        stages = calls[]
        linsolve = StageFailureLinearSolve(Ref(0), stages + 2, ReturnCode.Failure)
        sol = solve(
            remake(decay_prob; tspan = (0.0, direction * 0.3)), ctor(; linsolve);
            dt, dtmin = 0.1, dtmax = 0.1, force_dtmin = true, abstol = 1.0e3, reltol = 1.0e3,
            save_everystep = false, save_end = true, maxiters = 25
        )
        @test sol.retcode == ReturnCode.ConvergenceFailure
        @test sol.t[end] == dt
        @test sol.u[end] ≈ reference.u[end]
        @test sol.stats.naccept == 1
    end
end
