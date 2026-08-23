using Test, DelayDiffEq, DDEProblemLibrary
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
