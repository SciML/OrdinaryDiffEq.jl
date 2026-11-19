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
    rhs(u, p, t) = -8u
    prob = ODEProblem(rhs, 1.0, (0.0, 1.5))
    dt = 1 / 8
    integ = init(
        prob, ImplicitEuler(); adaptive = false, dt = dt,
        abstol = 1.0e-3, reltol = 1.0e-3
    )

    # Backward Euler gives u[n+1] = u[n] / (1 + 8dt) = u[n] / 2.
    expected = 1.0
    for n in 1:12
        previous_t = integ.t
        previous_iterations = integ.stats.nnonliniter
        step!(integ)

        # Bound the loop so a failed step cannot leave the test stuck at one time.
        @test integ.accept_step
        @test integ.t > previous_t
        @test integ.t ≈ n * dt

        iterations = integ.stats.nnonliniter - previous_iterations
        @test iterations >= 2

        expected /= 2
        @test integ.u ≈ expected atol = 1.0e-12
    end
    @test integ.stats.nreject == 0
end
