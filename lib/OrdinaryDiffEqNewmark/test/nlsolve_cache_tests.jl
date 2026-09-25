using OrdinaryDiffEqNewmark, Test
using NonlinearSolve: NewtonRaphson, RobustMultiNewton, FastShortcutNonlinearPolyalg,
    SimpleNewtonRaphson, SimpleTrustRegion
using SciMLBase

# u'' = -u with u(0) = 1, u'(0) = 0, so u(1) = cos(1) and u'(1) = -sin(1).
function accel!(dv, v, u, p, t)
    @. dv = -u
    return nothing
end
prob = SecondOrderODEProblem(accel!, [0.0], [1.0], (0.0, 1.0))

newmark_beta(nlsolve) = NewmarkBeta(; nlsolve)
generalized_alpha(nlsolve) = GeneralizedAlpha(rho_inf = 0.8; nlsolve)

schemes = ["NewmarkBeta" => newmark_beta, "GeneralizedAlpha" => generalized_alpha]

function discretization_params(alg::NewmarkBeta)
    return (; alg.β, alg.γ, αm = zero(alg.β), αf = zero(alg.β))
end
discretization_params(alg::GeneralizedAlpha) = (; alg.β, alg.γ, alg.αm, alg.αf)

# Inner solvers whose caches expose neither the parameter object nor the solved state at
# top level: `NonlinearSolvePolyAlgorithmCache` for the polyalgorithms, and
# `NonlinearSolveNoInitCache` for the SimpleNonlinearSolve algorithms, which have no
# `__init` of their own.
polyalgorithms = [
    "RobustMultiNewton" => RobustMultiNewton(),
    "FastShortcutNonlinearPolyalg" => FastShortcutNonlinearPolyalg(),
]
simple_solvers = [
    "SimpleNewtonRaphson" => SimpleNewtonRaphson(),
    "SimpleTrustRegion" => SimpleTrustRegion(),
]
inner_solvers = vcat(polyalgorithms, simple_solvers)

# The residual the schemes hand to the inner solver is linear in aₙ₊₁ for u'' = -u, so
# each step inverts in closed form. Consuming an explicit step sequence makes the
# reference sensitive to the dt used per step, which the solver reads out of its
# discretization cache.
function reference_endpoint(dts; β, γ, αm, αf)
    v, u = 0.0, 1.0
    for dt in dts
        aₙ = -u
        aₙ₊₁ = (
            -αm * aₙ - αf * u - (1 - αf) * (u + dt * v + dt^2 * (0.5 - β) * aₙ)
        ) / ((1 - αm) + (1 - αf) * β * dt^2)
        v, u = v + dt * ((1 - γ) * aₙ + γ * aₙ₊₁),
            u + dt * v + dt^2 / 2 * ((1 - 2β) * aₙ + 2β * aₙ₊₁)
    end
    return v, u
end

# The step-failure branch is unreachable from a problem every inner solver converges on,
# so exercising it needs an inner solve that reports failure on demand.
struct FailingNonlinearSolver <: SciMLBase.AbstractNonlinearAlgorithm end

mutable struct FailingNonlinearCache{P}
    prob::P
end

function SciMLBase.__init(
        prob::SciMLBase.AbstractNonlinearProblem, ::FailingNonlinearSolver, args...;
        kwargs...
    )
    return FailingNonlinearCache(prob)
end

function SciMLBase.reinit!(cache::FailingNonlinearCache, u0; p = cache.prob.p, kwargs...)
    cache.prob = SciMLBase.remake(cache.prob; u0, p)
    return cache
end

# Reports the failure alongside a finite, plausible-looking iterate, so a caller that
# ignores the retcode keeps integrating with a wrong acceleration instead of blowing up.
function SciMLBase.solve!(cache::FailingNonlinearCache)
    return SciMLBase.build_solution(
        cache.prob, FailingNonlinearSolver(), cache.prob.u0, copy(cache.prob.u0);
        retcode = ReturnCode.ConvergenceFailure
    )
end

@testset "Newmark inner nonlinear solvers" begin
    @testset "$scheme" for (scheme, make_alg) in schemes
        alg = make_alg(NewtonRaphson())
        params = discretization_params(alg)

        ref = solve(prob, alg, dt = 0.1, adaptive = false)
        ref_tstop = solve(prob, alg, dt = 0.1, adaptive = false, tstops = [0.25])
        @test ref.retcode == ReturnCode.Success
        @test ref_tstop.retcode == ReturnCode.Success

        # Absolute anchors on both components, so an error shared by every inner solver
        # still shows up.
        @test ref.u[end].x[1][1] ≈ -sin(1.0) atol = 2.0e-3
        @test ref.u[end].x[2][1] ≈ cos(1.0) atol = 2.0e-3

        # The tstop clips one step to dt = 0.05; pin both step sequences against the
        # closed-form map so a stale dt in the discretization cache cannot pass.
        @test sort(unique(round.(diff(ref_tstop.t), digits = 12))) == [0.05, 0.1]
        for sol in (ref, ref_tstop)
            v, u = reference_endpoint(diff(sol.t); params...)
            @test sol.u[end].x[1][1] ≈ v atol = 1.0e-12
            @test sol.u[end].x[2][1] ≈ u atol = 1.0e-12
        end

        @testset "$name" for (name, nlsolve) in inner_solvers
            sol = solve(prob, make_alg(nlsolve), dt = 0.1, adaptive = false)
            @test sol.retcode == ReturnCode.Success
            @test sol.u[end].x[1][1] ≈ ref.u[end].x[1][1] rtol = 1.0e-8
            @test sol.u[end].x[2][1] ≈ ref.u[end].x[2][1] rtol = 1.0e-8
        end

        @testset "$name across a clipped step" for (name, nlsolve) in polyalgorithms
            sol = solve(
                prob, make_alg(nlsolve), dt = 0.1, adaptive = false, tstops = [0.25]
            )
            @test sol.retcode == ReturnCode.Success
            @test sol.u[end].x[1][1] ≈ ref_tstop.u[end].x[1][1] rtol = 1.0e-8
            @test sol.u[end].x[2][1] ≈ ref_tstop.u[end].x[2][1] rtol = 1.0e-8
        end

        # The discretization state is one long-lived object that every step refreshes,
        # not a fresh one recycling buffers out of the inner cache. Stepping across a
        # tstop varies dt, which a cache refreshed only at init would miss.
        @testset "discretization cache is refreshed per step" begin
            integrator = init(prob, alg, dt = 0.1, adaptive = false, tstops = [0.25])
            evalcache = integrator.cache.evalcache
            step!(integrator)                     # 0.0 -> 0.1
            step!(integrator)                     # 0.1 -> 0.2
            step!(integrator)                     # 0.2 -> 0.25, dt clipped to 0.05
            @test integrator.cache.evalcache === evalcache
            @test evalcache.t ≈ 0.2 atol = 1.0e-12
            @test evalcache.dt ≈ 0.05 atol = 1.0e-12
            step!(integrator)                     # 0.25 -> 0.35, dt back to 0.1
            @test evalcache.t ≈ 0.25 atol = 1.0e-12
            @test evalcache.dt ≈ 0.1 atol = 1.0e-12
        end

        @testset "failing inner solve fails the step" begin
            sol = solve(
                prob, make_alg(FailingNonlinearSolver()), dt = 0.1,
                adaptive = false, maxiters = 20
            )
            @test !SciMLBase.successful_retcode(sol.retcode)
            @test sol.t[end] == 0.0
        end

        # `evalcache.p` has to track `integrator.p`: a callback may rebind it, and the
        # residual reads the parameters out of the discretization cache. The reference
        # is two chained constant-parameter solves, which a correct solver reproduces
        # exactly and a cache holding its construction-time parameters does not.
        @testset "rebound parameters reach the inner solve" begin
            function stiff_accel!(dv, v, u, p, t)
                @. dv = -p[1] * u
                return nothing
            end
            pprob = SecondOrderODEProblem(
                stiff_accel!, [0.0], [1.0], (0.0, 1.0), [1.0]
            )
            rebind = SciMLBase.DiscreteCallback(
                (u, t, integ) -> t ≥ 0.5 - 1.0e-12,
                integ -> (integ.p = [4.0]; nothing);
                save_positions = (false, false)
            )

            first_half = solve(
                SciMLBase.remake(pprob; tspan = (0.0, 0.5)), alg,
                dt = 0.1, adaptive = false
            )
            second_half = solve(
                SciMLBase.remake(
                    pprob; tspan = (0.5, 1.0), u0 = copy(first_half.u[end]), p = [4.0]
                ), alg, dt = 0.1, adaptive = false
            )
            vref = second_half.u[end].x[1][1]
            uref = second_half.u[end].x[2][1]

            solvers = vcat(["NewtonRaphson" => NewtonRaphson()], inner_solvers)
            @testset "$name" for (name, nlsolve) in solvers
                sol = solve(
                    pprob, make_alg(nlsolve), dt = 0.1,
                    adaptive = false, callback = rebind
                )
                @test sol.retcode == ReturnCode.Success
                @test sol.u[end].x[1][1] ≈ vref atol = 1.0e-12
                @test sol.u[end].x[2][1] ≈ uref atol = 1.0e-12
            end
        end
    end
end
