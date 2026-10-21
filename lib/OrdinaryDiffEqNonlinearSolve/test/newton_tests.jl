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
    # Orego's exact trajectory shifts slightly (nf 992192 -> 993159, +967 of ~1e6,
    # Lorenz is unaffected) now that relax!'s slope at α=0 is the analytic one
    # (dϕ0 = -ϕ0) rather than a finite difference: BackTracking occasionally
    # backtracks a handful of extra half-steps on an otherwise-identical trajectory
    # (nw/nreject/njacs below are unchanged). The +20 bound predates that; +2000
    # keeps the same intent (catch a real blowup, not a few extra evaluations).
    @test sol2.stats.nf <= sol1.stats.nf + 2000
    # A finite-difference `dϕ0` can be the wrong sign when the Newton correction is
    # tiny relative to the iterate late in a stiff solve; mistaking that for
    # divergence inflates W evaluations by two orders of magnitude by perturbing the
    # whole subsequent trajectory, not just the one step.
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

# `relax!`'s slope at α=0 is now the analytic one for the frozen-weight merit
# (dϕ0 = -ϕ0, since dz solves W*dz = R(z): see the comment in newton.jl), not a
# finite difference, so it is always a true descent direction for an exact,
# current-iterate W: `always_new = true` gets that and converges to the root, with
# or without an unrelated, large, unchanging padding component (zero correction)
# elsewhere in the state.
#
# For a non-`always_new` solve, W is linearized once per stage (typically at
# `uprev`), not at the current iterate, so the guarantee does not hold on a stage
# whose Jacobian there differs sharply from its Jacobian near the predictor/root:
# BackTracking can fail its own acceptance test and give up at a vacuous, near-zero
# step with the residual untouched, and nothing downstream of `relax!` recognizes
# that as non-convergence. This is a known, unresolved limitation, not covered by
# this fix: it reports a wrong `Success` rather than failing, and is NOT asserted
# against below. Fixing it needs `calc_J!`/`calc_W!` to support relinearizing at an
# arbitrary point independently of the static `always_new` flag, which does not
# currently exist and is out of scope here (see the PR discussion).
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

    alg = ImplicitEuler(
        nlsolve = NLNewton(
            relax = BackTracking(), predictor = zero_predictor!,
            always_new = true, max_iter = 50
        )
    )
    sol = solve(prob, alg; dt = 1.0, adaptive = false, abstol = 1.0e-3, reltol = 0.0)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test norm(Fstage(sol.u[end]), Inf) <= 1.0e-3

    function padded!(du, u, p, t)
        du[1:2] .= u[1:2] .- 1 .- Fstage(u[1:2])
        return du[3] = p
    end
    function padded_pred!(u, uprev, p, t, dt)
        u[1] = 0.0
        u[2] = 0.0
        return u[3] = uprev[3] + p * dt
    end
    for magnitude in (1.0e7, 1.0e9, 1.0e12)
        padded_prob = ODEProblem(padded!, [1.0, 1.0, 0.0], (0.0, 1.0), magnitude)
        padded_alg = ImplicitEuler(
            nlsolve = NLNewton(
                relax = BackTracking(), predictor = padded_pred!,
                always_new = true, max_iter = 50
            )
        )
        padded_sol = solve(padded_prob, padded_alg; dt = 1.0, adaptive = false, abstol = 1.0e-3, reltol = 0.0)
        @test padded_sol.retcode == SciMLBase.ReturnCode.Success
        @test norm(Fstage(padded_sol.u[end][1:2]), Inf) <= 1.0e-3
    end
end

# A fresh, exact Jacobian whose Newton correction is exactly representable (the
# linear stage solves in one step) but tiny or zero relative to the iterate: the
# analytic slope must not be confused by a vanishing or stationary probe the way a
# finite-difference one was (round 6).
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
