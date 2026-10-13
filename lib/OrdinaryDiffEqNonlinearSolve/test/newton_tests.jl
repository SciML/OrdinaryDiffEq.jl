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
