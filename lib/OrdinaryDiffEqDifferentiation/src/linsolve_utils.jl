issuccess_W(W::LinearAlgebra.Factorization) = LinearAlgebra.issuccess(W)
issuccess_W(W::Number) = !iszero(W)
issuccess_W(::Any) = true

function dolinsolve(integrator, linsolve; A = nothing, linu = nothing, b = nothing,
        du = nothing, u = nothing, p = nothing, t = nothing,
        weight = nothing, solverdata = nothing,
        reltol = integrator === nothing ? nothing : integrator.opts.reltol)
    A !== nothing && (linsolve.A = A)
    b !== nothing && (linsolve.b = b)
    linu !== nothing && (linsolve.u = linu)

    Plprev = linsolve.Pl isa LinearSolve.ComposePreconditioner ? linsolve.Pl.outer :
             linsolve.Pl
    Prprev = linsolve.Pr isa LinearSolve.ComposePreconditioner ? linsolve.Pr.outer :
             linsolve.Pr

    _alg = unwrap_alg(integrator, true)

    _Pl, _Pr = _alg.precs(linsolve.A, du, u, p, t, A !== nothing, Plprev, Prprev,
        solverdata)
    if (_Pl !== nothing || _Pr !== nothing)
        __Pl = _Pl === nothing ? SciMLOperators.IdentityOperator(length(integrator.u)) : _Pl
        __Pr = _Pr === nothing ? SciMLOperators.IdentityOperator(length(integrator.u)) : _Pr
        linsolve.Pl = __Pl
        linsolve.Pr = __Pr
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
