using OrdinaryDiffEqNonlinearSolve: NLNewton
using OrdinaryDiffEqCore
using OrdinaryDiffEqSDIRK
using DiffEqDevTools
using DiffEqBase
using LineSearches
using LinearAlgebra
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
    # Orego uses 967/992192 = 0.0975% extra RHS calls (966 nonlinear iterations
    # and one accepted step); a 0.1% work budget covers this measured overhead.
    nf_budget = prob === prob_ode_orego ? 0.001 * sol1.stats.nf : 20
    @test sol2.stats.nf <= sol1.stats.nf + nf_budget
    @test sol2.stats.nw <= sol1.stats.nw + 100
end

# The full Newton step solves this linear stage despite the overshooting predictor.
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

# Per-solve line-search overhead stays constant as the state size grows.
let
    function alloc_for(n, relax)
        A = Diagonal(fill(-1000.0, n))
        f!(du, u, p, t) = mul!(du, A, u)
        jac!(J, u, p, t) = (J .= A)
        fn = ODEFunction(f!, jac = jac!, jac_prototype = Matrix(A))
        prob = ODEProblem(fn, ones(n), (0.0, 0.1))
        alg = ImplicitEuler(nlsolve = NLNewton(; relax))
        solve(prob, alg, dt = 0.001, adaptive = false, save_everystep = false)
        return @allocated solve(prob, alg, dt = 0.001, adaptive = false, save_everystep = false)
    end
    extra(n) = alloc_for(n, BackTracking()) - alloc_for(n, nothing)
    small, large = extra(16), extra(4096)
    @test large < 10small
end

# Unrelated exceptions, invalid configurations and exhausted searches propagate.
module ForeignLineSearch
    struct LineSearchException <: Exception
        message::String
    end
end
struct AlwaysThrows end
(::AlwaysThrows)(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0) = throw(ForeignLineSearch.LineSearchException("unrelated"))

let
    decay!(du, u, p, t) = (du .= -100 .* u)
    alg(relax) = ImplicitEuler(nlsolve = NLNewton(; relax))
    prob = ODEProblem(decay!, [1.0], (0.0, 0.1))

    @test_throws ForeignLineSearch.LineSearchException solve(prob, alg(AlwaysThrows()); dt = 0.1, adaptive = false)

    @test_throws LineSearches.LineSearchException solve(
        prob, alg(MoreThuente(maxfev = 0)); dt = 0.1, adaptive = false
    )

    f!(du, u, p, t) = (du[1] = u[1] + 1 - exp(u[1]))
    overshoot_predictor!(u, uprev, p, t, dt) = (u[1] = -10)
    exhausted_alg = ImplicitEuler(
        nlsolve = NLNewton(
            relax = BackTracking(iterations = 1), always_new = true,
            max_iter = 50, predictor = overshoot_predictor!
        )
    )
    @test_throws LineSearches.LineSearchException solve(
        ODEProblem(f!, [1.0], (0.0, 1.0)), exhausted_alg; dt = 1.0, adaptive = false
    )
    # Control: the same stage with no iteration budget limit converges.
    full_alg = ImplicitEuler(
        nlsolve = NLNewton(
            relax = BackTracking(), always_new = true,
            max_iter = 50, predictor = overshoot_predictor!
        )
    )
    sol = solve(ODEProblem(f!, [1.0], (0.0, 1.0)), full_alg; dt = 1.0, adaptive = false)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test sol.u[end][1] ≈ log(2) rtol = 1.0e-6
end

include("stalled_linesearch_tests.jl")

# Exact linear corrections include zero steps and steps below the predictor's ULP.
let
    f!(du, u, p, t) = (du[1] = 2 * (u[1] - 1) + p)
    jac!(J, u, p, t) = (J[1, 1] = 2.0)
    unchanged!(u, uprev, p, t, dt) = (u .= uprev)
    prevfloat_pred!(u, uprev, p, t, dt) = (u .= prevfloat.(uprev))
    for c in (1.0e-8, 1.0e-9, 1.0e-10), predictor in (unchanged!, prevfloat_pred!),
            relax in (BackTracking(), LineSearches.Static(), nothing)

        alg = ImplicitEuler(nlsolve = NLNewton(; relax, predictor, always_new = true))
        sol = solve(
            ODEProblem(ODEFunction(f!; jac = jac!), [1.0], (0.0, 1.0), c), alg;
            dt = 1.0, adaptive = false, abstol = 1.0e-12, reltol = 0.0
        )
        @test sol.retcode == SciMLBase.ReturnCode.Success
        @test sol.u[end][1] ≈ 1 - c atol = 1.0e-12 rtol = 0.0
    end
end
