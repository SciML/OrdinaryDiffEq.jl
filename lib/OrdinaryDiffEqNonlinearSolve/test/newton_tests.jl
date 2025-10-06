using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using DiffEqDevTools
using DiffEqBase
using LineSearches
using Test

using ODEProblemLibrary: prob_ode_lorenz, prob_ode_orego

for prob in (prob_ode_lorenz, prob_ode_orego)
    sol1 = solve(prob, Trapezoid(), reltol = 1e-12, abstol = 1e-12)
    @test sol1.retcode == DiffEqBase.ReturnCode.Success
    sol2 = solve(prob, Trapezoid(nlsolve = NLNewton(relax = BackTracking())),
        reltol = 1e-12, abstol = 1e-12)
    @test sol2.retcode == DiffEqBase.ReturnCode.Success
    @test sol2.stats.nf <= sol1.stats.nf + 20
end
