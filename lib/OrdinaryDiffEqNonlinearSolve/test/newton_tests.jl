using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using DiffEqDevTools
using DiffEqBase
using LineSearches
using Test

using ODEProblemLibrary: prob_ode_lorenz, prob_ode_orego

for prob in (prob_ode_lorenz, prob_ode_orego)
    sol1 = solve(prob, Trapezoid(), reltol = 1.0e-12, abstol = 1.0e-12)
    @test sol1.retcode == SciMLBase.ReturnCode.Success
    sol2 = solve(
        prob, Trapezoid(nlsolve = NLNewton(relax = BackTracking())),
        reltol = 1.0e-12, abstol = 1.0e-12
    )
    @test sol2.retcode == SciMLBase.ReturnCode.Success
    @test sol2.stats.nf <= sol1.stats.nf + 20
end

@testset "non-descent is not nonlinear convergence" begin
    for T in (Float32, Float64)
        ϕ(α) = one(α) + α
        dϕ(α) = one(α)
        ϕdϕ(α) = (ϕ(α), dϕ(α))
        @test_throws ErrorException OrdinaryDiffEqNonlinearSolve.linesearch_step(
            BackTracking(), ϕ, dϕ, ϕdϕ, one(T)
        )
    end

    function rhs!(du, u, p, t)
        du[1] = (30 - p[1] * t) * u[1] + 1
    end
    prob = ODEProblem(rhs!, [0.0], (0.0, 0.1), [800.0])
    alg = ImplicitEuler(nlsolve = NLNewton(relax = BackTracking()))
    integrator = init(
        prob, alg; dt = 0.1, adaptive = false, abstol = 1.0e-12, reltol = 1.0e-12
    )
    @test_throws ErrorException solve!(integrator)
    @test integrator.t == 0.0

    for relaxation in (nothing, BackTracking())
        fresh = solve(
            prob, ImplicitEuler(nlsolve = NLNewton(; relax = relaxation, always_new = true));
            dt = 0.1, adaptive = false, abstol = 1.0e-12, reltol = 1.0e-12
        )
        @test SciMLBase.successful_retcode(fresh)
        @test fresh.u[end][1] ≈ 1 / 60
        @test abs(6 * fresh.u[end][1] - 0.1) < 1.0e-12
    end

    integrator.p[1] = 0.0
    reinit!(integrator, [0.0])
    reused = solve!(integrator)
    @test SciMLBase.successful_retcode(reused)
    @test reused.u[end][1] ≈ -0.05
    @test abs(-2 * reused.u[end][1] - 0.1) < 1.0e-12
    @test prob.u0 == [0.0]
end

@testset "non-negative FD slope is salvaged for BackTracking" begin
    for T in (Float32, Float64), slope in (T(1), T(NaN), T(Inf), T(-Inf))
        ϕ = α -> abs2(one(α) - oftype(α, 4) * α)
        dϕ = α -> oftype(α, -8) * (one(α) - oftype(α, 4) * α)
        ϕdϕ = α -> iszero(α) ? (ϕ(α), slope) : (ϕ(α), dϕ(α))

        ϕ0 = ϕ(zero(T))
        @test ϕ(one(ϕ0)) > ϕ0

        α = OrdinaryDiffEqNonlinearSolve.linesearch_step(
            BackTracking(), ϕ, dϕ, ϕdϕ, one(ϕ0)
        )
        @test α < one(α)
        @test ϕ(α) < ϕ0

        stalled(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0) = (eps(α0) / 2, ϕ(eps(α0) / 2))
        @test_throws ErrorException OrdinaryDiffEqNonlinearSolve.linesearch_step(
            stalled, ϕ, dϕ, ϕdϕ, one(ϕ0)
        )
    end

    for slope in (-2.0e-5, 0.0, 1.0)
        ϕ = α -> abs(1 - 2.0e-5 * α)
        dϕ = α -> -2.0e-5
        ϕdϕ = α -> (ϕ(α), slope)
        @test OrdinaryDiffEqNonlinearSolve.linesearch_step(
            BackTracking(), ϕ, dϕ, ϕdϕ, 1.0
        ) == 1.0
    end

    ϕ = α -> abs(α)
    dϕ = α -> one(α)
    ϕdϕ = α -> (ϕ(α), dϕ(α))
    @test iszero(
        OrdinaryDiffEqNonlinearSolve.linesearch_step(
            BackTracking(), ϕ, dϕ, ϕdϕ, 1.0
        )
    )
end
