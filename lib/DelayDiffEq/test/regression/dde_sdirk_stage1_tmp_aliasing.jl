using Test, DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqCore: IController, perform_step!
using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqSDIRK
using SciMLBase: ReturnCode

# Regression test for the nlsolver.tmp/uprev aliasing bug found while working
# on SciML/OrdinaryDiffEq.jl#3648.
#
# `get_tmp_cache` for Newton-based algorithms returns `cache.nlsolver.tmp`
# (OrdinaryDiffEqCore/src/integrators/integrator_interface.jl). DelayDiffEq's
# discontinuity tracking for state-dependent lags (`track.jl`,
# `discontinuity_function`) grabs that buffer via `get_tmp_cache(integrator)`
# and interpolates the delayed state into it in place. The implicit-first-stage
# branch of the ESDIRK/IMEX `perform_step!` used to assign
# `nlsolver.tmp = uprev` directly, aliasing the two arrays: any later write
# through `get_tmp_cache` silently corrupted `uprev`. With `dependent_lags`
# set, discontinuity tracking runs on every rejected adaptive step, so
# `MethodOfSteps(ImplicitEuler())` drives dt to zero and aborts with
# `ReturnCode.Unstable` instead of integrating normally.
@testset "implicit-first-stage nlsolver.tmp does not alias uprev" begin
    prob = remake(
        DDEProblemLibrary.prob_dde_constant_1delay_ip;
        constant_lags = nothing, dependent_lags = ((u, p, t) -> 1,)
    )

    integrator = init(prob, MethodOfSteps(ImplicitEuler()); dt = 0.1)
    step!(integrator)
    @test integrator.cache.nlsolver.tmp !== integrator.uprev

    sol = solve(prob, MethodOfSteps(ImplicitEuler()))
    @test sol.retcode == ReturnCode.Success
end

@testset "failed implicit first stage preserves uprev" begin
    f!(du, u, h, p, t) = (du[1] = -u[1]^2; nothing)
    history(p, t; idxs = nothing) = idxs === nothing ? [1.0] : 1.0
    prob = DDEProblem(
        f!, history, (0.0, 1.0); dependent_lags = ((u, p, t) -> 0.1,)
    )
    for method in (ImplicitEuler, SDIRK2)
        alg = method(nlsolve = NLNewton(max_iter = 1))
        integrator = init(prob, MethodOfSteps(alg); dt = 0.1, adaptive = false)
        uprev = copy(integrator.uprev)
        perform_step!(integrator, integrator.cache)
        @test integrator.force_stepfail
        @test integrator.cache.nlsolver.c == 1
        @test integrator.uprev == uprev
        @test integrator.cache.nlsolver.tmp !== integrator.uprev
        fill!(integrator.cache.nlsolver.tmp, 123)
        @test integrator.uprev == uprev

        integrator.cache.nlsolver.maxiters = 10
        integrator.force_stepfail = false
        perform_step!(integrator, integrator.cache)
        @test !integrator.force_stepfail
        @test integrator.uprev == uprev
        @test integrator.cache.nlsolver.tmp !== integrator.uprev
    end
end

@testset "multi-stage implicit-first DDE fixed point" begin
    τ = 0.01
    calls = Ref(0)
    function robertson!(du, u, h, p, t)
        calls[] += 1
        calls[] <= 100_000 || error("SDIRK2 exceeded the RHS evaluation budget")
        delayed_u2 = h(p, t - τ; idxs = 2)
        reaction = 0.04 * u[1] - 10_000 * delayed_u2 * u[3]
        decay = 3.0e7 * u[2]^2
        du[1] = -reaction
        du[2] = reaction - decay
        du[3] = decay
        return nothing
    end
    function robertson(u, h, p, t)
        du = similar(u)
        robertson!(du, u, h, p, t)
        return du
    end
    history(p, t; idxs = nothing) = idxs === nothing ? [1.0, 0.0, 0.0] : 0.0

    for f in (robertson!, robertson), growth in (nothing, 10, 10_000)
        calls[] = 0
        prob = DDEProblem(f, history, (0.0, 1.0); constant_lags = [τ])
        method = SDIRK2()
        controller = growth === nothing ? nothing :
            IController(method; qmax = 10, qmax_first_step = growth)
        sol = solve(
            prob, MethodOfSteps(method); controller, reltol = 0.1, abstol = 1.0e-4,
            dt = 1.0e-6, maxiters = 1000, save_everystep = false
        )

        @test sol.retcode == ReturnCode.Success
        @test sol.t[end] == 1.0
        @test all(isfinite, sol.u[end])
        @test sum(sol.u[end]) ≈ 1.0 atol = 1.0e-12
    end
end
