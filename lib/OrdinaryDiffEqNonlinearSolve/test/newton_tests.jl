using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqBDF: DFBDF
using DiffEqDevTools
using DiffEqBase
using LineSearches
using Test

using ODEProblemLibrary: prob_ode_lorenz, prob_ode_orego

@testset "a rebuilt W keeps the convergence rate by default" begin
    # u' = -u^2 with u(0) = 1, so u(1) = 1/2
    decay!(res, du, u, p, t) = (res[1] = du[1] + u[1]^2; nothing)
    prob = DAEProblem(decay!, [-1.0], [1.0], (0.0, 1.0); differential_vars = [true])
    integ = init(prob, DFBDF(); abstol = 1.0e-6, reltol = 1.0e-6)
    kept = 0
    for _ in 1:200
        integ.t < 1 || break
        nw, iters = integ.stats.nw, integ.stats.nnonliniter
        step!(integ)
        # Past the first steps, a solve with a new W can still stop after one iteration.
        kept += integ.iter > 10 && integ.stats.nw > nw &&
            integ.stats.nnonliniter == iters + 1
    end
    @test integ.t == 1
    @test integ.u[1] ≈ 0.5 atol = 1.0e-6
    @test kept > 0
end

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

        # Stop after the first rejected step so a failed solve cannot leave the
        # test repeating at the same time.
        if !integ.accept_step
            @test integ.accept_step
            break
        end

        @test integ.t > previous_t
        @test integ.t ≈ n * dt

        iterations = integ.stats.nnonliniter - previous_iterations
        @test iterations >= 2

        expected /= 2
        @test integ.u ≈ expected atol = 1.0e-12
    end
    @test integ.stats.nreject == 0
end
