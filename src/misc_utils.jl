# Default nlsolve behavior, should move to FiniteDiff.jl

Base.@pure determine_chunksize(u, alg::DiffEqBase.DEAlgorithm) =
    determine_chunksize(u, get_chunksize(alg))
Base.@pure function determine_chunksize(u, CS)
    if CS != 0
        return CS
    else
        return ForwardDiff.pickchunksize(length(u))
    end
end

macro swap!(x, y)
    quote
        local tmp = $(esc(x))
        $(esc(x)) = $(esc(y))
        $(esc(y)) = tmp
    end
end

macro cache(expr)
    name = expr.args[2].args[1].args[1]
    fields = [x for x in expr.args[3].args if typeof(x) != LineNumberNode]
    cache_vars = Expr[]
    jac_vars = Pair{Symbol,Expr}[]
    for x in fields
        if x.args[2] == :uType ||
           x.args[2] == :rateType ||
           x.args[2] == :kType ||
           x.args[2] == :uNoUnitsType
            push!(cache_vars, :(c.$(x.args[1])))
        elseif x.args[2] == :DiffCacheType
            push!(cache_vars, :(c.$(x.args[1]).du))
            push!(cache_vars, :(c.$(x.args[1]).dual_du))
        end
    end
    quote
        $expr
        $(esc(:full_cache))(c::$name) = tuple($(cache_vars...))
    end
end

# Nest one layer of value in order to get rid of possible Dual{Complex} or Complex{Dual} issues
# value should recurse for anything else.
function constvalue(::Type{T}) where {T}
    _T = DiffEqBase.value(T)
    return _T <: Complex ? DiffEqBase.value(real(_T)) : DiffEqBase.value(_T)
end
function constvalue(x)
    _x = DiffEqBase.value(x)
    return _x isa Complex ? DiffEqBase.value(real(_x)) : DiffEqBase.value(_x)
end

function diffdir(integrator::DiffEqBase.DEIntegrator)
    difference = maximum(abs, integrator.uprev) * sqrt(eps(typeof(integrator.t)))
    dir =
        integrator.tdir > zero(integrator.tdir) ?
        integrator.t > integrator.sol.prob.tspan[2] - difference ? -true : true :
        integrator.t < integrator.sol.prob.tspan[2] + difference ? true : -true
end

abstract type AbstractThreadingOption end
struct Sequential <: AbstractThreadingOption end
struct BaseThreads <: AbstractThreadingOption end
struct PolyesterThreads <: AbstractThreadingOption end

isthreaded(b::Bool) = b
isthreaded(::Sequential) = false
isthreaded(::BaseThreads) = true
isthreaded(::PolyesterThreads) = true

macro threaded(option, ex)
    quote
        opt = $(esc(option))
        if (opt === BaseThreads()) || ((opt isa Bool) && opt)
            $(esc(:(Threads.@threads :static $ex)))
        elseif opt === PolyesterThreads()
            $(esc(:(Polyester.@batch $ex)))
        else
            $(esc(ex))
        end
    end
end

function dolinsolve(
    integrator,
    linsolve;
    A = nothing,
    linu = nothing,
    b = nothing,
    du = nothing,
    u = nothing,
    p = nothing,
    t = nothing,
    weight = nothing,
    solverdata = nothing,
    reltol = integrator === nothing ? nothing : integrator.opts.reltol,
)

    A !== nothing && (linsolve = LinearSolve.set_A(linsolve, A))
    b !== nothing && (linsolve = LinearSolve.set_b(linsolve, b))
    linu !== nothing && (linsolve = LinearSolve.set_u(linsolve, linu))

    Plprev =
        linsolve.Pl isa LinearSolve.ComposePreconditioner ? linsolve.Pl.outer : linsolve.Pl
    Prprev =
        linsolve.Pr isa LinearSolve.ComposePreconditioner ? linsolve.Pr.outer : linsolve.Pr

    _alg = unwrap_alg(integrator, true)

    _Pl, _Pr =
        _alg.precs(linsolve.A, du, u, p, t, A !== nothing, Plprev, Prprev, solverdata)
    if (_Pl !== nothing || _Pr !== nothing)
        _weight =
            weight === nothing ?
            (linsolve.Pr isa Diagonal ? linsolve.Pr.diag : linsolve.Pr.inner.diag) : weight
        Pl, Pr = wrapprecs(_Pl, _Pr, _weight)
        linsolve = LinearSolve.set_prec(linsolve, Pl, Pr)
    end

    linres = if reltol === nothing
        solve(linsolve; reltol)
    else
        solve(linsolve; reltol)
    end

    # TODO: this ignores the add of the `f` count for add_steps!
    if integrator isa SciMLBase.DEIntegrator &&
       _alg.linsolve !== nothing &&
       !LinearSolve.needs_concrete_A(_alg.linsolve) &&
       linsolve.A isa WOperator &&
       linsolve.A.J isa SparseDiffTools.JacVec

        if alg_autodiff(_alg)
            integrator.destats.nf += linres.iters
        else
            integrator.destats.nf += 2 * linres.iters
        end
    end

    return linres
end

function wrapprecs(_Pl, _Pr, weight)
    if _Pl !== nothing
        Pl = LinearSolve.ComposePreconditioner(
            LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
            _Pl,
        )
    else
        Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight)))
    end

    if _Pr !== nothing
        Pr = LinearSolve.ComposePreconditioner(Diagonal(_vec(weight)), _Pr)
    else
        Pr = Diagonal(_vec(weight))
    end
    Pl, Pr
end

issuccess_W(W::LinearAlgebra.Factorization) = LinearAlgebra.issuccess(W)
issuccess_W(W::Number) = !iszero(W)
issuccess_W(::Any) = true
