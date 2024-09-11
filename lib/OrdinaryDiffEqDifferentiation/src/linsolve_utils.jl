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
            (;u, p, t) = integrator
            du = hasproperty(integrator, :du) ? integrator.du : nothing
            p = (du, u, p, t)
            reinit!(linsolve; A, p)
        else
            reinit!(linsolve; A)
        end
    end

    linres = solve!(linsolve; reltol)

    # TODO: this ignores the add of the `f` count for add_steps!
    if integrator isa SciMLBase.DEIntegrator && _alg.linsolve !== nothing &&
       !LinearSolve.needs_concrete_A(_alg.linsolve) &&
       linsolve.A isa WOperator && linsolve.A.J isa AbstractSciMLOperator
        if alg_autodiff(_alg) isa AutoForwardDiff
            integrator.stats.nf += linres.iters
        elseif alg_autodiff(_alg) isa AutoFiniteDiff
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2) * linres.iters
        else
            error("$alg_autodiff not yet supported in dolinsolve function")
        end
    end

    return linres
end

#for backward compat delete soon
function wrapprecs(PL, PR, weight, u)
    Pl = _Pl === nothing ? SciMLOperators.IdentityOperator(length(u)) : _Pl
        return linsolver
    Pr = _Pr === nothing ? SciMLOperators.IdentityOperator(length(u)) : _Pr
    end
    Pl, Pr
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
