module ExternalNonlinearSolverClient

    using CommonSolve: init, step!
    using DiffEqBase: DEVerbosity
    using OrdinaryDiffEqCore:
        AbstractNLSolver, Convergence, Divergence, FastConvergence, SlowConvergence,
        VerySlowConvergence, isfirststage
    using OrdinaryDiffEqNonlinearSolve:
        NLAnderson, NLFunctional, NLNewton, anderson, anderson!, build_nlsolver,
        can_smooth_est, compute_step!, du_alias_or_new, initial_η, markfirststage!,
        nlsolve!, nlsolvefail
    using OrdinaryDiffEqSDIRK: ImplicitEuler
    using SciMLBase: ODEProblem

    function build_solver(prob, nlalg, iip; dt = 0.1)
        alg = ImplicitEuler(nlsolve = nlalg)
        u = copy(prob.u0)
        return build_nlsolver(
            alg, u, copy(u), prob.p, first(prob.tspan), dt, prob.f, zero(u),
            Float64, Float64, Float64, 1.0, 1.0, iip, DEVerbosity()
        )
    end

    function solve_stage(prob, nlalg, iip; dt = 0.1)
        alg = ImplicitEuler(nlsolve = nlalg)
        integrator = init(prob, alg; dt, adaptive = false)
        nlsolver = build_solver(prob, nlalg, iip; dt)
        initial_rate = initial_η(nlsolver, integrator)
        z = nlsolve!(nlsolver, integrator)
        return (; nlsolver, z, initial_rate)
    end

    function solve_integrator_step(prob, nlalg; dt = 0.1)
        integrator = init(prob, ImplicitEuler(nlsolve = nlalg); dt, adaptive = false)
        step!(integrator)
        return integrator.u
    end

end

using .ExternalNonlinearSolverClient
using Test

const ExternalClient = ExternalNonlinearSolverClient

@testset "Developer nonlinear-solver API" begin
    oop_prob = ExternalClient.ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    iip_prob = ExternalClient.ODEProblem(
        (du, u, p, t) -> (du .= -u), [1.0], (0.0, 1.0)
    )

    functional = ExternalClient.build_solver(
        iip_prob, ExternalClient.NLFunctional(), Val(true)
    )
    @test functional isa ExternalClient.AbstractNLSolver
    @test !ExternalClient.can_smooth_est(functional)
    @test ExternalClient.markfirststage!(functional) === nothing

    du = ExternalClient.du_alias_or_new(functional, zeros(1))
    du_again = ExternalClient.du_alias_or_new(functional, zeros(1))
    @test du === du_again
    du[1] = 2.0
    @test ExternalClient.du_alias_or_new(functional, zeros(1)) == [2.0]

    newton = ExternalClient.build_solver(iip_prob, ExternalClient.NLNewton(), Val(true))
    @test ExternalClient.can_smooth_est(newton)
    @test ExternalClient.markfirststage!(newton) === nothing
    @test ExternalClient.isfirststage(newton)

    oop_newton_integrator = ExternalClient.init(
        oop_prob, ExternalClient.ImplicitEuler(); dt = 0.1, adaptive = false
    )
    oop_newton = ExternalClient.build_solver(
        oop_prob, ExternalClient.NLNewton(), Val(false)
    )
    for invalid_cache in (0.5, nothing)
        @test_throws ArgumentError ExternalClient.nlsolve!(
            oop_newton, oop_newton_integrator, invalid_cache, false
        )
    end

    @test !ExternalClient.nlsolvefail(ExternalClient.FastConvergence)
    @test !ExternalClient.nlsolvefail(ExternalClient.Convergence)
    @test ExternalClient.nlsolvefail(ExternalClient.SlowConvergence)
    @test ExternalClient.nlsolvefail(ExternalClient.VerySlowConvergence)
    @test ExternalClient.nlsolvefail(ExternalClient.Divergence)

    alg = ExternalClient.NLFunctional(κ = 1.0e-10, max_iter = 20)
    integrator = ExternalClient.init(
        oop_prob, ExternalClient.ImplicitEuler(nlsolve = alg); dt = 0.1,
        adaptive = false
    )
    single_step_solver = ExternalClient.build_solver(oop_prob, alg, Val(false))
    @test ExternalClient.initial_η(single_step_solver, integrator) >= 0
    @test isfinite(ExternalClient.compute_step!(single_step_solver, integrator))

    oop_result = ExternalClient.solve_stage(oop_prob, alg, Val(false))
    @test !ExternalClient.nlsolvefail(oop_result.nlsolver)
    @test oop_result.initial_rate >= 0
    @test oop_result.z ≈ -1 / 11 atol = 1.0e-10

    anderson_alg = ExternalClient.NLAnderson(
        κ = 1.0e-10, max_iter = 20, max_history = 3, aa_start = 1
    )
    oop_anderson = ExternalClient.solve_stage(oop_prob, anderson_alg, Val(false))
    @test !ExternalClient.nlsolvefail(oop_anderson.nlsolver)
    @test oop_anderson.z ≈ -1 / 11 atol = 1.0e-10

    @test ExternalClient.solve_integrator_step(iip_prob, anderson_alg) ≈ [10 / 11] atol =
        1.0e-10
end
