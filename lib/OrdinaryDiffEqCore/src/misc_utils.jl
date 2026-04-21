macro swap!(x, y)
    return quote
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
    return quote
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

function diffdir(integrator::SciMLBase.DEIntegrator)
    difference = maximum(abs, integrator.uprev) * sqrt(eps(typeof(integrator.t)))
    return dir = integrator.tdir > zero(integrator.tdir) ?
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

@inline function _threaded_execute(f, ::Union{BaseThreads, Bool}, range)
    return Threads.@threads :static for i in range
        f(i)
    end
end

function _polyester_batch(args...)
    throw(
        ArgumentError(
            LazyString(
                "PolyesterThreads() requires Polyester.jl to be loaded. ",
                "Add `using Polyester` to your code."
            )
        )
    )
end

@inline function _threaded_execute(f, ::PolyesterThreads, range)
    return _polyester_batch(f, range)
end

macro threaded(option, ex)
    ex.head === :for || error("@threaded expects a for loop")
    loop_var = esc(ex.args[1].args[1])
    range_expr = esc(ex.args[1].args[2])
    body = esc(ex.args[2])
    return quote
        opt = $(esc(option))
        if isthreaded(opt)
            _threaded_execute(opt, $range_expr) do $loop_var
                $body
            end
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
    return esc(:(Base.@assume_effects :foldable $arg))
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
    if hasproperty(f, :mass_matrix)
        mm = f.mass_matrix
        mm = mm isa MatrixOperator ? mm.A : mm

        if mm isa UniformScaling
            return nothing
        elseif all(!iszero, mm)
            return trues(size(mm, 1))
        elseif !(mm isa SciMLOperators.AbstractSciMLOperator) && _isdiag(mm)
            return reshape(diag(mm) .!= 0, size(u))
        else
            return DifferentialVarsUndefined()
        end
    else
        return nothing
    end
end

# Fallback for _isdiag - uses LinearAlgebra.isdiag which is O(n²)
# Sparse specialization is provided in OrdinaryDiffEqCoreSparseArraysExt
_isdiag(A::AbstractMatrix) = isdiag(A)

isnewton(::Any) = false

# Extract the chunk size integer from an ADType for use as a type parameter.
# Returns 0 when the chunk size should be automatically determined.
_ad_chunksize_int(::AutoForwardDiff{nothing}) = 0
_ad_chunksize_int(::AutoForwardDiff{CS}) where {CS} = CS
_ad_chunksize_int(::AutoSparse{<:AutoForwardDiff{nothing}}) = 0
_ad_chunksize_int(::AutoSparse{<:AutoForwardDiff{CS}}) where {CS} = CS
_ad_chunksize_int(_) = 0

# Extract the finite difference type from an ADType for use as a type parameter.
_ad_fdtype(::AutoFiniteDiff{FD}) where {FD} = FD
_ad_fdtype(::AutoSparse{<:AutoFiniteDiff{FD}}) where {FD} = FD
_ad_fdtype(_) = Val{:forward}()

# Fix AutoFiniteDiff dir: default dir of true (Bool) makes integration non-reversible
function _fixup_ad(ad_alg::AutoFiniteDiff)
    if ad_alg.dir isa Bool
        @reset ad_alg.dir = Int(ad_alg.dir)
    end
    return ad_alg
end
function _fixup_ad(ad_alg::Bool)
    throw(
        ArgumentError(
            "Passing a `Bool` for keyword argument `autodiff` is no longer supported. " *
                "Use an `ADType` specifier from ADTypes.jl, e.g. `AutoForwardDiff()` or `AutoFiniteDiff()`."
        )
    )
end
_fixup_ad(ad_alg) = ad_alg

# The effective linear solver for an algorithm: prefer the one attached to the
# inner nlsolve object (canonical forward-compatible location, matching the
# NonlinearSolveAlg / NonlinearSolve.jl convention), falling back to the
# algorithm-level `linsolve` field (historical location) when the nlsolve
# doesn't carry one of its own. Rosenbrock-family algorithms don't have a
# nlsolve and are handled by the final fallback.
function effective_linsolve(alg)
    if hasfield(typeof(alg), :nlsolve)
        nl_linsolve = _nlsolve_linsolve(alg.nlsolve)
        nl_linsolve === nothing || return nl_linsolve
    end
    return hasfield(typeof(alg), :linsolve) ? alg.linsolve : nothing
end
_nlsolve_linsolve(nlalg) = hasfield(typeof(nlalg), :linsolve) ? nlalg.linsolve : nothing
