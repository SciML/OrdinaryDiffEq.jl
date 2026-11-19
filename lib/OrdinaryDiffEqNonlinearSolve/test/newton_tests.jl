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
    # Use a sufficiently small fixed step that the nonlinear iteration reaches the
    # branch governed by solve-level adaptivity after the first accepted step.
    rhs(u, p, t) = -10u^2
    prob = ODEProblem(rhs, 1.0, (0.0, 0.5))
    dt = 0.01
    integ = init(
        prob, ImplicitEuler(); adaptive = false, dt = dt,
        abstol = 1.0e-3, reltol = 1.0e-3
    )

    # Backward Euler's scalar update solves x + 10dt*x^2 = u_n.
    expected = 1.0
    nsteps = 0
    while integ.t < 0.5
        previous_u = integ.u
        previous_iterations = integ.stats.nnonliniter
        step!(integ)
        iterations = integ.stats.nnonliniter - previous_iterations
        nsteps += 1

        # Every accepted step has a meaningful state correction in this problem.
        if integ.accept_step && abs(integ.u - previous_u) > 1.0e-10
            @test iterations >= 2
        end
        expected = 2expected / (1 + sqrt(1 + 40dt * expected))
    end
    @test SciMLBase.successful_retcode(integ.sol)
    @test nsteps == 50
    @test integ.u ≈ expected atol = 1.0e-6
end
