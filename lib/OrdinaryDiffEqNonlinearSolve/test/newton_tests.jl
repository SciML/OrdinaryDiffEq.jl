using OrdinaryDiffEqNonlinearSolve: NLNewton
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


@testset "fixed-step ImplicitEuler respects solve-level adaptivity (#4905)" begin
    # ImplicitEuler supports adaptive stepping by default, but this solve explicitly
    # disables adaptivity. The nonlinear iteration must use the resolved solve option,
    # not the algorithm's default adaptivity trait.
    rhs(u, p, t) = -10u^2
    prob = ODEProblem(rhs, 1.0, (0.0, 0.5))
    dt = 0.05
    sol = solve(
        prob, ImplicitEuler(); adaptive = false, dt = dt,
        abstol = 1.0e-12, reltol = 1.0e-12
    )

    # Backward Euler's scalar update solves x + 10dt*x^2 = u_n.
    expected = 1.0
    for _ in 1:10
        expected = 2expected / (1 + sqrt(1 + 40dt * expected))
    end
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test length(sol.t) == 11
    @test sol.u[end] ≈ expected atol = 1.0e-9
end
