using OrdinaryDiffEqNonlinearSolve: NLNewton, relax!
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

# LineSearches throws when the Newton step is not a descent direction;
# relax! must keep the undamped step there instead of letting the
# exception escape solve.
@testset "relax! keeps the Newton step along a non-descent direction" begin
    integrator = init(
        prob_ode_orego, Trapezoid(nlsolve = NLNewton(relax = BackTracking())),
        reltol = 1.0e-12, abstol = 1.0e-12
    )
    SciMLBase.step!(integrator)
    nlsolver = integrator.cache.nlsolver
    f = OrdinaryDiffEqCore.nlsolve_f(integrator)
    dz0 = ones(length(integrator.u))
    dz = copy(dz0)
    @test relax!(dz, nlsolver, integrator, f) === dz
    @test dz == dz0
end

# Backward Euler for u' = -2 - u^3 from u = 0 with dt = 1 solves u^3 + u + 2 = 0,
# whose unique real root is -1.
# The Newton solve stops at abstol = 1e-6, so the endpoint is only determined to that tolerance.
@testset "Float32 BackTracking damps Newton steps" begin
    f!(du, u, p, t) = (du[1] = -2 - u[1]^3)
    prob = ODEProblem(f!, [0.0f0], (0.0f0, 1.0f0))
    sol = solve(
        prob, ImplicitEuler(nlsolve = NLNewton(relax = BackTracking(), max_iter = 50));
        adaptive = false, dt = 1.0f0, abstol = 1.0f-6, reltol = 1.0f-6
    )
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test sol.t[end] == 1.0f0
    @test sol.u[end][1] ≈ -1.0f0 atol = 1.0f-6
end
