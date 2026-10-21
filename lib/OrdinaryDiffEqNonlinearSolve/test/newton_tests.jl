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
    @test sol2.stats.nf <= sol1.stats.nf + 20
    # A wrong-sign `dϕ0` from the fixed-step finite-difference probe (the Newton
    # correction can be tiny relative to the iterate late in a stiff solve) must not
    # be treated as divergence: that inflates W evaluations by two orders of
    # magnitude by perturbing the whole subsequent trajectory, not just the one step.
    @test sol2.stats.nw <= sol1.stats.nw + 100
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

# The extra allocation from the line-search weight buffer must not scale with
# state size: it is a cache field, recomputed in place, not reallocated per call.
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
    # An O(n) per-iteration reallocation would grow this by ~256x; a cache-backed,
    # in-place buffer leaves only O(1) closure overhead.
    @test large < 10small
end

# `relax!` must not catch anything beyond the non-descent precondition: a same-named
# but unrelated exception type propagates, as does a genuine line-search
# configuration error or an exhausted search budget.
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

# A stage whose Jacobian at the predictor/root differs sharply from its Jacobian at
# `uprev` (so a reused, non-`always_new` W is a poor local model there): a well-resolved
# positive slope must not be accepted as a full step just because it reads as a false
# positive elsewhere (Orego's case, above). `always_new = true` gets a Jacobian
# evaluated at the current iterate and converges to the unique root; a reused W cannot
# make that same guarantee and must not silently report a wrong answer as `Success`.
let
    ϵ = 1.0e-5
    Bmat = [1.0 2 / ϵ; 0.0 1.0]
    Amat = Bmat * [1.0 0.0; -ϵ 1 - (1 - ϵ)]
    cvec = [1.0, 0.0]
    hswitch(x) = x <= 0.25 ? 0.0 : x >= 0.75 ? 1.0 : (y = 2 * (x - 0.25); y * y * (3 - 2y))
    Fstage(u) = Amat * u + cvec + hswitch(u[2]) * (Bmat - Amat) * (u .- 1)
    f!(du, u, p, t) = (du .= u .- 1 .- Fstage(u))
    zero_predictor!(u, uprev, p, t, dt) = fill!(u, 0.0)
    prob = ODEProblem(f!, ones(2), (0.0, 1.0))

    fresh_alg = ImplicitEuler(
        nlsolve = NLNewton(
            relax = BackTracking(), predictor = zero_predictor!,
            always_new = true, max_iter = 50
        )
    )
    sol_fresh = solve(prob, fresh_alg; dt = 1.0, adaptive = false, abstol = 1.0e-3, reltol = 0.0)
    @test sol_fresh.retcode == SciMLBase.ReturnCode.Success
    @test sol_fresh.u[end] ≈ [-1.0, -1.0] atol = 1.0e-3

    reused_alg = ImplicitEuler(
        nlsolve = NLNewton(
            relax = BackTracking(), predictor = zero_predictor!,
            always_new = false, max_iter = 50
        )
    )
    sol_reused = solve(prob, reused_alg; dt = 1.0, adaptive = false, abstol = 1.0e-3, reltol = 0.0)
    @test sol_reused.retcode != SciMLBase.ReturnCode.Success ||
        isapprox(sol_reused.u[end], [-1.0, -1.0], atol = 1.0e-3)
end
