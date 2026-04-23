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

# `isapprox` for timeseries solutions.
#
# The generic `Base.isapprox(::AbstractArray, ::AbstractArray)` computes
# `norm(x - y) <= max(atol, rtol*max(norm(x), norm(y)))`. When an
# `AbstractTimeseriesSolution` (e.g. `ODESolution`) is iterated via the
# `AbstractArray` fallback, iteration can yield an object whose type equals the
# container type (because the solution subtypes `AbstractVectorOfArray{T,N}` and
# its leaf container is a `VectorOfArray` whose `eltype` resolves to a scalar
# while the iterator itself yields an array-like element). Under
# `RecursiveArrayTools` v4 / Julia 1.12, `LinearAlgebra.norm(itr)` first calls
# `norm_recursive_check` which throws
#   ArgumentError: cannot evaluate norm recursively if the type of the initial
#   element is identical to that of the container
# for these nested / `StructArray`-backed storages. Compare snapshot-by-snapshot
# instead, which also sidesteps the cross-container-type issue when comparing
# solutions whose element storage differs (e.g. `Vector{SVector}` vs
# `StructArray{SVector}`).
function Base.isapprox(
        x::SciMLBase.AbstractTimeseriesSolution,
        y::SciMLBase.AbstractTimeseriesSolution;
        atol::Real = 0,
        rtol::Real = Base.rtoldefault(
            _cmp_eltype(x), _cmp_eltype(y), atol
        ),
        nans::Bool = false,
        norm = LinearAlgebra.norm
    )
    length(x.u) == length(y.u) || return false
    isapprox(x.t, y.t; atol = atol, rtol = rtol, nans = nans) || return false
    for i in eachindex(x.u)
        _isapprox_snapshot(x.u[i], y.u[i]; atol, rtol, nans, norm) || return false
    end
    return true
end

# Element-wise comparison of solution snapshots. Snapshots may be plain arrays
# or `VectorOfArray`s. Compare by reducing to a materialised numeric array whose
# `eltype` is a scalar or small static vector, so `LinearAlgebra.norm` does not
# recurse into a container whose element type equals itself.
_isapprox_snapshot(a, b; kwargs...) = isapprox(a, b; kwargs...)
function _isapprox_snapshot(
        a::RecursiveArrayTools.AbstractVectorOfArray,
        b::RecursiveArrayTools.AbstractVectorOfArray;
        kwargs...
    )
    return _isapprox_materialised(_snapshot_array(a), _snapshot_array(b); kwargs...)
end

_snapshot_array(a) = a
_snapshot_array(a::RecursiveArrayTools.AbstractVectorOfArray) = _snapshot_array(parent(a))

# Materialise into `Array` so we don't rely on arithmetic dispatch for exotic
# backing types (e.g. `StructArray{<:SVector}`), which can produce a container
# whose eltype equals its own type under `-`.
function _isapprox_materialised(a::AbstractArray, b::AbstractArray; kwargs...)
    size(a) == size(b) || return false
    return isapprox(collect(a), collect(b); kwargs...)
end

_cmp_eltype(sol::SciMLBase.AbstractTimeseriesSolution) = recursive_bottom_eltype(sol.u)

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
