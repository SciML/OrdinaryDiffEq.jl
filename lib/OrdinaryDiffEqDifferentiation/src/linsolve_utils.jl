issuccess_W(W::LinearAlgebra.Factorization) = LinearAlgebra.issuccess(W)
issuccess_W(W::Number) = !iszero(W)
issuccess_W(::Any) = true

function dolinsolve(integrator, linsolve; A = nothing, linu = nothing, b = nothing,
        reltol = integrator === nothing ? nothing : integrator.opts.reltol)
    A !== nothing && (linsolve.A = A)
    b !== nothing && (linsolve.b = b)
    linu !== nothing && (linsolve.u = linu)

    _alg = unwrap_alg(integrator, true)
    if !isnothing(A)
        (;du, u, p, t) = integrator
        p = isnothing(integrator) ? nothing : (du, u, p, t)
        reinit!(linsolve; A, p)
    end

    linres = solve!(linsolve; reltol)

    # TODO: this ignores the add of the `f` count for add_steps!
    if integrator isa SciMLBase.DEIntegrator && _alg.linsolve !== nothing &&
       !LinearSolve.needs_concrete_A(_alg.linsolve) &&
       linsolve.A isa WOperator && linsolve.A.J isa AbstractSciMLOperator
        if alg_autodiff(_alg) isa AutoForwardDiff
            integrator.stats.nf += linres.iters
        elseif alg_autodiff(_alg) isa AutoFiniteDiff
            integrator.stats.nf += 2 * linres.iters
        else
            error("$alg_autodiff not yet supported in dolinsolve function")
        end
    end

    return linres
end

function wrapprecs(_Pl::Nothing, _Pr::Nothing, weight, u)
    Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight)))
    Pr = Diagonal(_vec(weight))
    Pl, Pr
end

function wrapprecs(_Pl, _Pr, weight, u)
    Pl = _Pl === nothing ? SciMLOperators.IdentityOperator(length(u)) : _Pl
    Pr = _Pr === nothing ? SciMLOperators.IdentityOperator(length(u)) : _Pr
    Pl, Pr
end

Base.resize!(p::LinearSolve.LinearCache, i) = p
