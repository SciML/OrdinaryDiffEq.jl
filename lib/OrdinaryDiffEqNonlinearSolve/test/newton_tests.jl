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

# A predictor that overshoots the stage far enough that the finite-difference slope
# at α = 0 reads as ascent, even though the full step solves the stage exactly.
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

# A predictor far enough from the stage's unique root that the full Newton step
# overshoots badly (merit function increases), but a damped step reaches it. Checked
# against both a line search with a descent precondition (BackTracking) and one
# without (Static), which must not be overridden by non-descent handling.
let
    F(u) = u * ((u + 2.5)^2 + 71 / 36)
    f!(du, u, p, t) = (du[1] = u[1] - 1 - F(u[1]))
    overshoot_predictor!(u, uprev, p, t, dt) = (u[1] = -2)
    for relax in (BackTracking(), LineSearches.Static())
        alg = ImplicitEuler(
            nlsolve = NLNewton(; relax, always_new = true, max_iter = 50, predictor = overshoot_predictor!)
        )
        sol = solve(
            ODEProblem(f!, [1.0], (0.0, 1.0)), alg;
            dt = 1.0, adaptive = false, abstol = 1.0e-9, reltol = 1.0e-6
        )
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.u[end][1] ≈ 0 atol = 1.0e-6
    end
end
