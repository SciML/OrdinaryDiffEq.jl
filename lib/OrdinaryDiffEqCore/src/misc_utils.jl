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
    jac_vars = Pair{Symbol, Expr}[]
    for x in fields
        if x.args[2] == :uType || x.args[2] == :rateType ||
           x.args[2] == :kType || x.args[2] == :uNoUnitsType
            push!(cache_vars, :(c.$(x.args[1])))
        elseif x.args[2] == :DiffCacheType
            push!(cache_vars, :(c.$(x.args[1]).du))
            push!(cache_vars, :(c.$(x.args[1]).dual_du))
        end
    end
    quote
        $(esc(expr))
        $(esc(:full_cache))(c::$(esc(name))) = tuple($(cache_vars...))
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
    dir = integrator.tdir > zero(integrator.tdir) ?
          integrator.t > integrator.sol.prob.tspan[2] - difference ? -1 : 1 :
          integrator.t < integrator.sol.prob.tspan[2] + difference ? 1 : -1
end

error_constant(integrator, order) = error_constant(integrator, integrator.alg, order)

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

macro OnDemandTableauExtract(S_T, T, T2)
    S = getproperty(__module__, S_T)
    s = gensym(:s)
    q = quote
        $s = $S($T, $T2)
    end
    fn = fieldnames(S)
    for n in fn
        push!(q.args, Expr(:(=), n, Expr(:call, :getfield, s, QuoteNode(n))))
    end
    return esc(q)
end
macro OnDemandTableauExtract(S_T, T)
    S = getproperty(__module__, S_T)
    s = gensym(:s)
    q = quote
        $s = $S($T)
    end
    fn = fieldnames(S)
    for n in fn
        push!(q.args, Expr(:(=), n, Expr(:call, :getfield, s, QuoteNode(n))))
    end
    return esc(q)
end

macro fold(arg)
    # https://github.com/JuliaLang/julia/pull/43852
    if VERSION < v"1.8.0-DEV.1484"
        esc(:(@generated $arg))
    else
        esc(:(Base.@assume_effects :foldable $arg))
    end
end

struct DifferentialVarsUndefined end

"""
    get_differential_vars(f, idxs, timeseries::uType)

Returns an array of booleans for which values are the differential variables
vs algebraic variables. Returns `nothing` for the cases where all variables
are differential variables. Returns `DifferentialVarsUndefined` if it cannot
be determined (i.e. the mass matrix is not diagonal).
"""
function get_differential_vars(f, u)
    differential_vars = nothing
    if hasproperty(f, :mass_matrix)
        mm = f.mass_matrix
        mm = mm isa MatrixOperator ? mm.A : mm

        if mm isa UniformScaling || all(!iszero, mm)
            return nothing
        elseif !(mm isa SciMLOperators.AbstractSciMLOperator) && isdiag(mm)
            differential_vars = reshape(diag(mm) .!= 0, size(u))
        else
            return DifferentialVarsUndefined()
        end
    end
end

isnewton(::Any) = false
