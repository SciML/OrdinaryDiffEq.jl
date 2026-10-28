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
            BackTracking(), ϕ, dϕ, ϕdϕ, one(T), [one(T)], [one(T)]
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
            BackTracking(), ϕ, dϕ, ϕdϕ, one(ϕ0), [one(ϕ0)], [one(ϕ0)]
        )
        @test α < one(α)
        @test ϕ(α) < ϕ0

        stalled(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0) = (eps(α0) / 2, ϕ(eps(α0) / 2))
        @test_throws ErrorException OrdinaryDiffEqNonlinearSolve.linesearch_step(
            stalled, ϕ, dϕ, ϕdϕ, one(ϕ0), [one(ϕ0)], [one(ϕ0) / 2]
        )
    end

    for slope in (-2.0e-5, 0.0, 1.0)
        ϕ = α -> abs(1 - 2.0e-5 * α)
        dϕ = α -> -2.0e-5
        ϕdϕ = α -> (ϕ(α), slope)
        @test OrdinaryDiffEqNonlinearSolve.linesearch_step(
            BackTracking(), ϕ, dϕ, ϕdϕ, 1.0, [1.0], [1.0]
        ) == 1.0
    end

    ϕ = α -> abs(α)
    dϕ = α -> one(α)
    ϕdϕ = α -> (ϕ(α), dϕ(α))
    @test iszero(
        OrdinaryDiffEqNonlinearSolve.linesearch_step(
            BackTracking(), ϕ, dϕ, ϕdϕ, 1.0, [1.0], [1.0]
        )
    )
end

@testset "nearby line-search minimum" begin
    function halving_search(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
        iszero(ϕ0) && return zero(α0), ϕ0
        α = α0
        while α > eps(α0)
            ϕα = ϕ(α)
            ϕα < ϕ0 && return α, ϕα
            α /= 2
        end
        error("No decreasing step found")
    end

    @testset "helper $T" for (T, root) in (
            (Float32, sqrt(eps(Float32)) / 16), (Float64, 1.0e-9),
        )
        ϕ(α) = abs(one(T) - α / root)
        dϕ(α) = (ϕ(α + sqrt(eps(T))) - ϕ(α)) / sqrt(eps(T))
        ϕdϕ(α) = (ϕ(α), dϕ(α))
        @test dϕ(zero(T)) > 0
        for search in (halving_search, BackTracking())
            α = OrdinaryDiffEqNonlinearSolve.linesearch_step(
                search, ϕ, dϕ, ϕdϕ, one(T), [one(T)], [one(T)]
            )
            @test α > 0
            @test ϕ(α) < ϕ(zero(T))
        end
    end

    @testset "analytic implicit stage" begin
        function rhs!(du, u, p, t)
            du[1] = (10 - 1.0e-8 - 80t) * u[1] + 1
        end
        prob = ODEProblem(rhs!, [0.0], (0.0, 0.1))
        sol = solve(
            prob, ImplicitEuler(nlsolve = NLNewton(relax = halving_search, max_iter = 100));
            dt = 0.1, adaptive = false, abstol = 1.0e-12, reltol = 0.0
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end][1] ≈ 0.1 / (0.8 + 1.0e-9)
        @test abs((0.8 + 1.0e-9) * sol.u[end][1] - 0.1) < 1.0e-12
    end
end

@testset "sub-epsilon line-search minimum" begin
    # A Jacobian stale by ~eps makes the Newton direction overshoot so far that
    # the line minimum sits below eps(1); only sub-epsilon steps make progress.
    function deep_halving(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
        iszero(ϕ0) && return zero(α0), ϕ0
        α = α0
        while α > eps(α0)^2
            ϕα = ϕ(α)
            ϕα < ϕ0 && return α, ϕα
            α /= 2
        end
        error("No decreasing step found")
    end

    for (T, detune, stiffen, abstol, exact) in (
            (Float32, 1.0e-6, 8000, 1.0e-6, 0.00125),
            (Float64, 1.0e-15, 80, 1.0e-12, 0.125),
            (Float64, 1.0e-15, 8000, 1.0e-12, 0.00125),
        )
        function rhs!(du, u, p, t)
            du[1] = (10 - p[1] - stiffen * t) * u[1] + 1
        end
        prob = ODEProblem(rhs!, T[0], (T(0), T(0.1)), T[detune])
        for search in (deep_halving, BackTracking())
            sol = solve(
                prob, ImplicitEuler(nlsolve = NLNewton(relax = search, max_iter = 100));
                dt = T(0.1), adaptive = false, abstol = abstol, reltol = 0.0
            )
            @test SciMLBase.successful_retcode(sol)
            @test sol.u[end][1] ≈ exact
        end
    end
end
