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
    dolinsolve(integrator, linsolve; A = nothing, linu = nothing, b = nothing, reltol = …) -> linres

Solve the linear system with the LinearSolve.jl cache `linsolve`, optionally
resetting its matrix `A`, unknown `linu`, right-hand side `b`, and tolerance
`reltol`. Updates the RHS-evaluation count for matrix-free/iterative solvers and
returns the LinearSolve result.
"""
function dolinsolve(
        integrator, linsolve; A = nothing, linu = nothing, b = nothing,
        reltol = integrator === nothing ? nothing : integrator.opts.reltol
    )
    b !== nothing && (linsolve.b = b)
    linu !== nothing && (linsolve.u = linu)

    _alg = unwrap_alg(integrator, true)
    if !isnothing(A)
        if integrator isa DEIntegrator
            (; u, p, t) = integrator
            du = hasproperty(integrator, :du) ? integrator.du : nothing
            linsolve.p = (du, u, p, t)
        end
        LinearSolve.reinit!(linsolve; A)
    end

    linres = solve!(linsolve; reltol)

    # TODO: this ignores the add of the `f` count for add_steps!
    if integrator isa SciMLBase.DEIntegrator && _alg.linsolve !== nothing &&
            !LinearSolve.needs_concrete_A(_alg.linsolve) &&
            linsolve.A isa WOperator && linsolve.A.J isa AbstractSciMLOperator
        ad = alg_autodiff(_alg) isa ADTypes.AutoSparse ? ADTypes.dense_ad(alg_autodiff(_alg)) :
            alg_autodiff(_alg)
        if ad isa ADTypes.AutoFiniteDiff || ad isa ADTypes.AutoFiniteDifferences
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2 * linres.iters)
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, linres.iters)
        end
    end

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
