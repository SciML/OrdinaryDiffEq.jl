"""
    issuccess_W(W) -> Bool

Return whether the factorized system matrix `W` is nonsingular / the factorization
succeeded. For a `Factorization` it forwards to `LinearAlgebra.issuccess`; for a
scalar `W` it checks `!iszero(W)`; otherwise it returns `true`.
"""
issuccess_W(W::LinearAlgebra.Factorization) = LinearAlgebra.issuccess(W)
issuccess_W(W::Number) = !iszero(W)
issuccess_W(::Any) = true

"""
    set_linear_reltol!(linsolve, reltol; default_algorithm = false) -> linsolve

Point the LinearSolve.jl cache `linsolve` at relative tolerance `reltol`.

`SciMLBase.solve!(::LinearCache, alg; kwargs...)` takes its tolerances from the cache and
drops any handed to it as keywords, so `solve!(linsolve; reltol)` sets nothing and an
iterative solve keeps running at LinearSolve's `sqrt(eps)` default however tight or loose
the integrator is. `LinearSolve.update_tolerances!` is the entry point that does take
effect.

Concrete factorization algorithms are left alone because they have no tolerance to set.
When `default_algorithm` is true, the cache may select an iterative algorithm whose
tolerance must be updated; factorization paths safely ignore it. A non-scalar `reltol`
(per-component tolerances) is left alone too — `LinearCache` holds a single scalar
tolerance and the iterative solvers compare against a scalar bound.
"""
function set_linear_reltol!(linsolve, reltol; default_algorithm = false)
    if reltol isa Number &&
            (default_algorithm || !LinearSolve.needs_concrete_A(linsolve.alg)) &&
            linsolve.reltol != reltol
        LinearSolve.update_tolerances!(linsolve; reltol)
    end
    return linsolve
end

"""
    ScaledResidualTolerance(b₁, work, reltol, recompute)

Absolute tolerance `reltol * ‖Pl \\ b₁‖` for an iterative solve, measured in the left
preconditioner `Pl` the solve will use. `work` is scratch of `b₁`'s size. The norm is taken
only once `Pl` is installed, since a `precs` callback may return a different `Pl` for the
system being solved, and is recomputed when `recompute` is true or the preconditioners were
refreshed for this solve; otherwise the cache keeps the tolerance of the earlier solve.
"""
struct ScaledResidualTolerance{B, W, T}
    b₁::B
    work::W
    reltol::T
    recompute::Bool
end

"""
    refresh_precs!(linsolve) -> Bool

Install the preconditioners `linsolve.alg.precs` returns for the cache's current `A` and `p`
if they are stale, as `solve!` would, and return whether it did. `solve!` then finds them
fresh, so a `precs` callback runs once per refresh and the tolerances of the solve can be
measured in the preconditioner it uses.
"""
function refresh_precs!(linsolve)
    (hasproperty(linsolve.alg, :precs) && !isnothing(linsolve.alg.precs)) || return false
    linsolve.precsisfresh || return false
    linsolve.Pl, linsolve.Pr = linsolve.alg.precs(linsolve.A, linsolve.p)
    linsolve.precsisfresh = false
    return true
end

set_linear_abstol!(linsolve, abstol::Number, refreshed) =
    LinearSolve.update_tolerances!(linsolve; abstol)
function set_linear_abstol!(linsolve, tol::ScaledResidualTolerance, refreshed)
    (tol.recompute || refreshed) || return nothing
    LinearAlgebra.ldiv!(tol.work, linsolve.Pl, tol.b₁)
    return LinearSolve.update_tolerances!(linsolve; abstol = tol.reltol * LinearAlgebra.norm(tol.work))
end

"""
    linear_residual!(integrator, linsolve, rlin, work, b₁) -> (; residual, bnorm, roundoff, rtol)

Leave the true residual `b - A u` of `linsolve` in `rlin` and return its norm `residual`
in the left preconditioner `Pl`, with `bnorm = ‖Pl \\ b‖`, the `roundoff`
`100 eps (‖Pl \\ b‖ + ‖Pl \\ A u‖ + ‖Pl \\ b₁‖)` of forming it and of evaluating `b`
itself, which is at the scale of the first residual `b₁`, and the cache's relative
tolerance `rtol`. `work` is scratch of `b`'s size. Charges `stats.nf` for the product with `A`.

`A` is applied to `u / ‖u‖` and the product scaled by `‖u‖`, so a finite-difference `A`
perturbs by the same size as in the solve. A finite-difference `A` is linear only to about
`sqrt(eps)`, so `rtol` is floored at `cbrt(eps)` for it; an exact `A` keeps `reltol`.
"""
function linear_residual!(integrator, linsolve, rlin, work, b₁)
    (; A, b, u, Pl, reltol) = linsolve
    T = real(eltype(b))
    unorm = LinearAlgebra.norm(u)
    Auscaled = zero(T)
    if iszero(unorm)
        copyto!(rlin, b)
    else
        # `A` is linear, so apply it to the unit vector `u / ‖u‖`, as the solve does to its
        # basis vectors: a finite-difference step sized for a tiny `u` can round away.
        @. work = u / unorm
        LinearAlgebra.mul!(rlin, A, work)
        drain_jvp_count!(integrator, unwrap_alg(integrator, true), A)
        @. rlin = unorm * rlin
        LinearAlgebra.ldiv!(work, Pl, rlin)
        Auscaled = LinearAlgebra.norm(work)
        @. rlin = b - rlin
    end
    LinearAlgebra.ldiv!(work, Pl, rlin)
    residual = LinearAlgebra.norm(work)
    LinearAlgebra.ldiv!(work, Pl, b)
    bnorm = LinearAlgebra.norm(work)
    LinearAlgebra.ldiv!(work, Pl, b₁)
    roundoff = 100 * eps(T) * (bnorm + Auscaled + LinearAlgebra.norm(work))
    rtol = uses_finite_differences(A) ? max(reltol, cbrt(eps(T))) : reltol
    return (; residual, bnorm, roundoff, rtol)
end

"""
    linear_residual_converged!(integrator, linsolve, rlin, work, b₁; evidence = false) -> Bool

Whether `linsolve.u` solves `linsolve.A * u = linsolve.b` to the cache's relative tolerance on
its true residual (see `linear_residual!`): `residual ≤ rtol * ‖Pl \\ b‖ + roundoff`.
`abstol` is deliberately not part of the bound: a Krylov solve may stop on it, but an
allowance taken from another residual can swallow the whole of this one, so it must not
decide accuracy. A solver can also report convergence on its recurrence residual while the
iterate is far from the solution: a badly scaled preconditioner can end it at `u = 0` after
one iteration. With `evidence = true` the tolerance is `sqrt(rtol)`, a guard against gross
failure that accepts an iterate that clearly reduced the residual.
"""
function linear_residual_converged!(integrator, linsolve, rlin, work, b₁; evidence = false)
    (; residual, bnorm, roundoff, rtol) = linear_residual!(integrator, linsolve, rlin, work, b₁)
    return residual <= (evidence ? sqrt(rtol) : rtol) * bnorm + roundoff
end

"""
    uses_finite_differences(W) -> Bool

Whether the Jacobian-vector products of `W` are finite differences, which are linear only to
about `sqrt(eps)`.
"""
uses_finite_differences(W) = jvp_counter(W) isa JVPCache && jvp_counter(W).jvp_op isa DAESplitPushforward

"""
    dolinsolve(integrator, linsolve; alg = unwrap_alg(integrator, true), A = nothing, linu = nothing, b = nothing, reltol = …, abstol = nothing) -> linres

Solve the linear system with the LinearSolve.jl cache `linsolve`, optionally
resetting its matrix `A`, unknown `linu`, right-hand side `b`, and tolerances
`reltol` (see `set_linear_reltol!`) and `abstol` (left unchanged when `nothing`; set only
for an iterative solver; a `ScaledResidualTolerance` is evaluated once the preconditioners
the solve uses are installed). `alg` is the integrator's concrete member
algorithm; callers holding a cache that stores it (e.g. `cache.alg`) should
pass it so the call stays a static dispatch. Charges `stats.nf` for the
Jacobian-vector products the solve applied (see `drain_jvp_count!`) and
returns the LinearSolve result.
"""
function dolinsolve(
        integrator, linsolve;
        alg = integrator === nothing ? nothing : unwrap_alg(integrator, true),
        A = nothing, linu = nothing, b = nothing,
        reltol = integrator === nothing ? nothing : integrator.opts.reltol,
        abstol = nothing
    )
    b !== nothing && (linsolve.b = b)
    linu !== nothing && (linsolve.u = linu)
    return _dolinsolve(alg, integrator, linsolve, A, reltol, abstol)
end

function _dolinsolve(_alg, integrator, linsolve, A, reltol, abstol)
    if !isnothing(A)
        if integrator isa DEIntegrator
            (; u, p, t) = integrator
            du = hasproperty(integrator, :du) ? integrator.du : nothing
            linsolve.p = (du, u, p, t)
        end
        LinearSolve.reinit!(linsolve; A)
    end

    set_linear_reltol!(
        linsolve, reltol;
        default_algorithm = hasproperty(_alg, :linsolve) && isnothing(_alg.linsolve)
    )
    if abstol !== nothing && !LinearSolve.needs_concrete_A(linsolve.alg)
        set_linear_abstol!(linsolve, abstol, refresh_precs!(linsolve))
    end
    linres = solve!(linsolve)

    # TODO: this ignores the add of the `f` count for add_steps!
    drain_jvp_count!(integrator, _alg, linsolve.A)

    return linres
end

"""
    default_krylov_warm_start(linsolver) -> linsolver

Resolve a Krylov linear solver left at `LinearSolve.WarmStart.Auto` to the mode
appropriate for a Newton-based integrator, `LinearSolve.WarmStart.Hegedus`
(Hegedüs-scaled reuse of the previous solution — the recommended mode for the
sequence of correlated preconditioned solves inside a Newton iteration).

Only the `Auto` default is resolved: an explicit `WarmStart.None`/`Previous`/
`Hegedus` or a non-Krylov solver is returned unchanged. This is called on the
Newton nonlinear-solver path only; Rosenbrock/W-method integrators never call
it, so their `Auto` solver stays a cold start (warm starting is unsafe there —
no outer Newton iteration absorbs the within-tolerance stage-solve
perturbation).

Warm starting the Newton path is only sound because LinearSolve discards a
Hegedüs guess that cannot reduce the residual (SciML/LinearSolve.jl#1123). The
guess is a previous Newton increment, and those converge to round-off; without
that check the rescaling amplified last-ulp input differences by ~1e13
(SciML/OrdinaryDiffEq.jl#4034). The `[compat]` floor on LinearSolve is what
guarantees the check is present.
"""
function default_krylov_warm_start(linsolver)
    (
        linsolver isa LinearSolve.KrylovJL &&
            linsolver.warm_start === LinearSolve.WarmStart.Auto
    ) || return linsolver
    return SciMLBase.remake(linsolver; warm_start = LinearSolve.WarmStart.Hegedus)
end

"""
    wrapprecs(linsolver, W, weight) -> linsolver

Attach a diagonal (`weight`-based) left/right preconditioner to `linsolver` when it
supports `precs` and none was supplied, returning the reconfigured solver;
otherwise return `linsolver` unchanged.
"""
function wrapprecs(linsolver, W, weight)
    if hasproperty(linsolver, :precs) && isnothing(linsolver.precs)
        Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight)))
        Pr = Diagonal(_vec(weight))
        precs = Returns((Pl, Pr))
        return remake(linsolver; precs)
    else
        return linsolver
    end
end

Base.resize!(p::LinearSolve.LinearCache, i) = p
