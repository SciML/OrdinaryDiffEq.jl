issuccess_W(W::LinearAlgebra.Factorization) = LinearAlgebra.issuccess(W)
issuccess_W(W::Number) = !iszero(W)
issuccess_W(::Any) = true

function dolinsolve(integrator, linsolve; A = nothing, linu = nothing, b = nothing,
        reltol = integrator === nothing ? nothing : integrator.opts.reltol)
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
