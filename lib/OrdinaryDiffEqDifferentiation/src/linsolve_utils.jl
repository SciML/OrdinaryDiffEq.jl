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
    dolinsolve(integrator, linsolve; alg = unwrap_alg(integrator, true), A = nothing, linu = nothing, b = nothing, reltol = …, abstol = nothing) -> linres

Solve the linear system with the LinearSolve.jl cache `linsolve`, optionally
resetting its matrix `A`, unknown `linu`, right-hand side `b`, and tolerances
`reltol` (see `set_linear_reltol!`) and `abstol`. `abstol` is set only for an iterative
solver and left unchanged when `nothing`. Besides a number it may be a function
`abstol(linsolve)`, which is evaluated after `A` is set and the cache's preconditioners
are refreshed for it (see `refresh_precs!`), so that a tolerance measured with
`linsolve.Pl` is in the units of the preconditioner the solve applies. `alg` is the
integrator's concrete member algorithm; callers holding a cache that stores it (e.g.
`cache.alg`) should pass it so the call stays a static dispatch. Charges `stats.nf` for the
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
        if !(abstol isa Number)
            refresh_precs!(linsolve)
            abstol = abstol(linsolve)
        end
        LinearSolve.update_tolerances!(linsolve; abstol)
    end
    linres = solve!(linsolve)

    # TODO: this ignores the add of the `f` count for add_steps!
    drain_jvp_count!(integrator, _alg, linsolve.A)

    return linres
end

"""
    refresh_precs!(linsolve) -> linsolve

Recompute `linsolve.Pl` and `linsolve.Pr` from the algorithm's `precs` if a new `A` or `p`
has marked them stale, as an iterative `solve!` would on entry; that `solve!` then reuses
them.
"""
function refresh_precs!(linsolve)
    precs = hasproperty(linsolve.alg, :precs) ? linsolve.alg.precs : nothing
    if linsolve.precsisfresh && precs !== nothing
        Pl, Pr = precs(linsolve.A, linsolve.p)
        linsolve.Pl = Pl
        linsolve.Pr = Pr
        linsolve.precsisfresh = false
    end
    return linsolve
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

For a DAE Newton matrix `W` (a `JVPCache` or `DAEWOperator`), a user `precs` is wrapped in
the weights instead: the left preconditioner applies `weight .* (Pl \\ x)` and the right one
`Pr \\ (x ./ weight)`. The Krylov solve then measures the weighted norm of the
preconditioned residual, as IDA does, with or without a user `precs`.
"""
function wrapprecs(linsolver, W, weight)
    hasproperty(linsolver, :precs) || return linsolver
    w = _vec(weight)
    if isnothing(linsolver.precs)
        Pl = LinearSolve.InvPreconditioner(Diagonal(w))
        Pr = Diagonal(w)
        precs = Returns((Pl, Pr))
        return remake(linsolver; precs)
    elseif W isa Union{JVPCache, DAEWOperator}
        user_precs = linsolver.precs
        precs = function (A, p)
            Pl, Pr = user_precs(A, p)
            return weight_left(Pl, w), weight_right(Pr, w)
        end
        return remake(linsolver; precs)
    else
        return linsolver
    end
end

_isidentity(P) = P isa UniformScaling && isone(P.λ)

"""
    ComposedLdiv(first, second)

Preconditioner whose `ldiv!` applies `first \\ x` and then `second \\` the result.
"""
struct ComposedLdiv{F, S}
    first::F
    second::S
end

Base.eltype(P::ComposedLdiv) = promote_type(eltype(P.first), eltype(P.second))

function LinearAlgebra.ldiv!(P::ComposedLdiv, x)
    ldiv!(P.first, x)
    return ldiv!(P.second, x)
end

function LinearAlgebra.ldiv!(y, P::ComposedLdiv, x)
    ldiv!(y, P.first, x)
    return ldiv!(P.second, y)
end

function weight_left(Pl, w)
    weight = LinearSolve.InvPreconditioner(Diagonal(w))
    return _isidentity(Pl) ? weight : ComposedLdiv(Pl, weight)
end

weight_right(Pr, w) = _isidentity(Pr) ? Diagonal(w) : ComposedLdiv(Diagonal(w), Pr)

Base.resize!(p::LinearSolve.LinearCache, i) = p
