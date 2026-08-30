using Test, DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqCore: IController
using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqSDIRK
using SciMLBase: ReturnCode, get_tmp_cache

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
    tmpcache = get_tmp_cache(integrator)
    @test tmpcache !== nothing
    @test first(tmpcache) !== integrator.uprev

    sol = solve(prob, MethodOfSteps(ImplicitEuler()))
    @test sol.retcode == ReturnCode.Success
end

@testset "multi-stage implicit-first tmp is not uprev" begin
    function f!(du, u, h, p, t)
        du[1] = -u[1]
        return nothing
    end
    history(p, t; idxs = nothing) = idxs === nothing ? [1.0] : 1.0
    prob = DDEProblem(f!, history, (0.0, 0.2); constant_lags = [0.1])
    # Same `nlsolve` as DelayDiffEq's working SDIRK MethodOfSteps tests.
    alg = SDIRK2(nlsolve = NLNewton(fast_convergence_cutoff = 0))
    integrator = init(prob, MethodOfSteps(alg); dt = 0.05, adaptive = false, maxiters = 20)
    step!(integrator)
    @test integrator.cache.nlsolver.tmp !== integrator.uprev
    tmpcache = get_tmp_cache(integrator)
    @test tmpcache !== nothing
    @test first(tmpcache) !== integrator.uprev
end

# SciML/OrdinaryDiffEq.jl#4389: MethodOfSteps(SDIRK2()) on a delayed Robertson
# problem spun inside LinearSolve QR when the first-step qmax was the ordinary
# `qmax` rather than 10_000. Recomputing the inner SDIRK step with a frozen `W`
# (ODE `repeat_step=true`) is wrong once MethodOfSteps has updated the history.
@testset "multi-stage implicit-first DDE fixed point" begin
    τ = 0.01
    function robertson!(du, u, h, p, t)
        delayed_u2 = h(p, t - τ; idxs = 2)
        reaction = 0.04 * u[1] - 10_000 * delayed_u2 * u[3]
        decay = 3.0e7 * u[2]^2
        du[1] = -reaction
        du[2] = reaction - decay
        du[3] = decay
        return nothing
    end
    history(p, t; idxs = nothing) = idxs === nothing ? [1.0, 0.0, 0.0] : 0.0

    prob = DDEProblem(robertson!, history, (0.0, 0.5); constant_lags = [τ])
    method = SDIRK2()
    controller = IController(method; qmax = 10, qmax_first_step = 10)
    sol = solve(
        prob, MethodOfSteps(method); controller, reltol = 0.1, abstol = 1.0e-4,
        dt = 1.0e-6, maxiters = 1000, save_everystep = false
    )

    @test sol.retcode == ReturnCode.Success
end
