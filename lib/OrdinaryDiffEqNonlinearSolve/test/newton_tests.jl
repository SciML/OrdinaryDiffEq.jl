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

# A predictor that overshoots the stage far enough that the exact Newton step's
# finite-difference slope at α = 0 reads as ascent, even though the full step solves
# the stage exactly. relax! must not report this as divergence (#4894 round 1).
let
    decay!(du, u, p, t) = (du .= -100 .* u)
    overshoot_predictor!(u, uprev, p, t, dt) = (u .= -2 .* uprev)
    exact = 1 / 11
    for adaptive in (false, true)
        alg = ImplicitEuler(nlsolve = NLNewton(relax = BackTracking(), predictor = overshoot_predictor!))
        sol = solve(
            ODEProblem(decay!, [1.0], (0.0, 0.1)), alg;
            dt = 0.1, adaptive, abstol = 1.0e-9, reltol = 1.0e-6
        )
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.u[end][1] ≈ exact rtol = 1.0e-6
    end
end
