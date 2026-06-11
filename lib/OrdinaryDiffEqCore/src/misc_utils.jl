macro swap!(x, y)
    return quote
        local tmp = $(esc(x))
        $(esc(x)) = $(esc(y))
        $(esc(y)) = tmp
    end
end

"""
    @cache struct MyCache ... end

Macro used to define a mutable solver cache. It emits the given `struct` definition
and additionally generates a `full_cache(c::MyCache)` method returning the tuple of
its resizable buffer fields (those typed `uType`, `rateType`, `kType`,
`uNoUnitsType`, or the `du`/`dual_du` of a `DiffCacheType`). That `full_cache`
tuple is what the `resize!`/`deleteat!` integrator interface iterates over when the
state length changes.
"""
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
        elseif x.args[1] == :tmp_cache
            # The unified scratch field is detected by name (`tmp_cache`) rather
            # than by type, so the cache may declare it as a concrete
            # `TmpCache{...}` or via a type parameter (to allow opted-out slots).
            # It expands into its sub-buffers so they show up in `full_cache`
            # (used by resize!, etc.) just like inline scratch fields used to.
            # Opted-out slots are `nothing`; `full_cache` consumers skip those.
            push!(cache_vars, :(c.$(x.args[1]).tmp))
            push!(cache_vars, :(c.$(x.args[1]).tmp2))
            push!(cache_vars, :(c.$(x.args[1]).atmp))
            push!(cache_vars, :(c.$(x.args[1]).rate_tmp))
            push!(cache_vars, :(c.$(x.args[1]).rate_tmp2))
        end
    end
    return quote
        $(esc(expr))
        $(GlobalRef(SciMLBase, :full_cache))(c::$(esc(name))) = tuple($(cache_vars...))
    end
end

# Nest one layer of value in order to get rid of possible Dual{Complex} or Complex{Dual} issues
# value should recurse for anything else.
"""
    constvalue(x)

Strip any ForwardDiff/unit wrapper from `x` (or a type `T`) down to its underlying
numeric value, taking the real part for `Complex` so that a scalar constant can be
compared/used unambiguously. Used e.g. for eigenvalue estimates.
"""
function constvalue(::Type{T}) where {T}
    _T = SciMLBase.value(T)
    return _T <: Complex ? SciMLBase.value(real(_T)) : SciMLBase.value(_T)
end
function constvalue(x)
    _x = SciMLBase.value(x)
    return _x isa Complex ? SciMLBase.value(real(_x)) : SciMLBase.value(_x)
end

"""
    diffdir(integrator) -> Int

Return the finite-difference direction (`+1` or `-1`) to use for time
derivatives, chosen so the stencil stays inside the integration interval near an
endpoint.
"""
function diffdir(integrator::SciMLBase.DEIntegrator)
    difference = maximum(abs, integrator.uprev) * sqrt(eps(typeof(integrator.t)))
    return dir = integrator.tdir > zero(integrator.tdir) ?
        integrator.t > integrator.sol.prob.tspan[2] - difference ? -1 : 1 :
        integrator.t < integrator.sol.prob.tspan[2] + difference ? 1 : -1
end

"""
    error_constant(integrator, order) -> Real

Return the leading error constant of the current method at the given `order`, used
when scaling the local error estimate. Dispatches on `integrator.alg`.
"""
error_constant(integrator, order) = error_constant(integrator, integrator.alg, order)

"""
    AbstractThreadingOption

Abstract supertype for the `thread = …` option controlling internal broadcasting.
The concrete choices are [`Sequential`](@ref), [`BaseThreads`](@ref), and
[`PolyesterThreads`](@ref); [`isthreaded`](@ref) reports whether a given option
enables multithreading.
"""
abstract type AbstractThreadingOption end
"""
    Sequential() <: AbstractThreadingOption

Threading option that disables internal multithreading — all per-element work runs
on a single thread. [`isthreaded`](@ref)`(Sequential())` is `false`.
"""
struct Sequential <: AbstractThreadingOption end
"""
    BaseThreads() <: AbstractThreadingOption

Threading option that parallelizes internal broadcasting with Julia's built-in
`Threads.@threads`. [`isthreaded`](@ref)`(BaseThreads())` is `true`.
"""
struct BaseThreads <: AbstractThreadingOption end
"""
    PolyesterThreads() <: AbstractThreadingOption

Threading option that parallelizes internal broadcasting with Polyester.jl's
low-overhead `@batch`. [`isthreaded`](@ref)`(PolyesterThreads())` is `true`.
"""
struct PolyesterThreads <: AbstractThreadingOption end

"""
    isthreaded(opt) -> Bool

Return whether the threading option `opt` enables multithreaded internal
broadcasting. `true` for [`BaseThreads`](@ref)/[`PolyesterThreads`](@ref) (and
`true` for a `Bool` `opt` equal to `true`), `false` for [`Sequential`](@ref).
"""
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

"""
    DifferentialVarsUndefined

Sentinel returned by [`get_differential_vars`](@ref) when the differential vs
algebraic split cannot be determined (the mass matrix is not diagonal). In that
case dense output falls back to linear interpolation.
"""
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

        if mm isa UniformScaling || mm isa SciMLOperators.IdentityOperator
            return nothing
        elseif mm isa SciMLOperators.ScalarOperator
            # λ·I: every variable is algebraic when λ is zero.
            return iszero(mm.val) ? falses(size(u)) : nothing
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

"""
    find_algebraic_vars_eqs(M)

Find algebraic variables (zero columns) and algebraic equations (zero rows) from mass matrix.
Returns `(algebraic_vars, algebraic_eqs)` as boolean arrays (true = algebraic).

Works on CPU and GPU arrays. Sparse specialization (O(nnz)) is provided in
OrdinaryDiffEqCoreSparseArraysExt.
"""
function find_algebraic_vars_eqs(M::Diagonal)
    _idxs = map(iszero, diag(M))
    return _idxs, _idxs
end

function find_algebraic_vars_eqs(M::AbstractMatrix)
    algebraic_vars = vec(all(iszero, M, dims = 1))
    algebraic_eqs = vec(all(iszero, M, dims = 2))
    return algebraic_vars, algebraic_eqs
end

function find_algebraic_vars_eqs(M::SciMLOperators.AbstractSciMLOperator)
    return find_algebraic_vars_eqs(convert(AbstractMatrix, M))
end

"""
    isnewton(nlsolver) -> Bool

Return whether the nonlinear solver `nlsolver` is a Newton-type solver (as opposed
to a fixed-point / functional iteration), which determines whether a W-matrix is
formed and updated.
"""
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
"""
    _fixup_ad(ad, args...)

Internal helper that adjusts an autodiff choice `ad` to be consistent with the
problem/solver context (e.g. disabling AD when it is not applicable). Returns the
possibly-modified autodiff choice.
"""
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

# Warm-start state for the scalar interval searches in `ode_interpolation` /
# `ode_interpolation!`. The strategy is stored as a mutable
# `FindFirstFunctions.StrategyKind` enum field, so it can be re-selected as
# the time grid's structure becomes known (at solve init from
# `saveat`/`adaptive`, on grid growth, and at the ending phase) without
# changing the container's type. All strategies return exact
# `searchsortedfirst`/`searchsortedlast` results; the kind only affects
# lookup speed. Races on the mutable fields under concurrent interpolation
# only degrade the starting guess, never correctness.
mutable struct TsSearchHint{T <: AbstractVector}
    const ts::T
    # 0 means "no query yet": the first search then starts from `lastindex(ts)`
    # at query time. `ts` grows after construction, and mid-solve consumers
    # (delay-equation history lookups especially) query near the current end,
    # so a guess frozen at construction would go stale.
    idx_prev::Int
    kind::StrategyKind
    # `length(ts)` at the last grid probe. `typemax(Int)` disables re-probing
    # (the grid's uniformity is already known from the solve options).
    probed_len::Int
end

TsSearchHint(ts::AbstractVector) = TsSearchHint(ts, 0, KIND_BRACKET_GALLOP, 0)

@inline function ts_hint_start(h::TsSearchHint, v)
    prev = h.idx_prev
    return ifelse(prev == 0, lastindex(v), prev)
end

# Sampled uniformity probe: O(64) regardless of grid size. Interpolation
# search wins on near-uniform grids (`saveat` ranges, fixed-dt stepping,
# smooth adaptive solves); the hinted bracketing gallop is the robust choice
# for irregular grids under the correlated access patterns of adjoints and
# history lookups.
function _ts_grid_kind(ts::AbstractVector)
    n = length(ts)
    n < 4 && return KIND_BRACKET_GALLOP
    tf = first(ts)
    tl = last(ts)
    mean_dt = (tl - tf) / (n - 1)
    iszero(mean_dt) && return KIND_BRACKET_GALLOP
    # Skip the first and last few intervals: the adaptive initial-dt ramp-up
    # and the final truncated step are structural boundary artifacts, not
    # indicative of the interior grid.
    skip = min(4, (n - 1) ÷ 4)
    navail = n - 1 - 2 * skip
    navail < 1 && return KIND_BRACKET_GALLOP
    nsamples = min(navail, 64)
    stride = navail ÷ nsamples
    @inbounds for k in 0:(nsamples - 1)
        i = firstindex(ts) + skip + k * stride
        r = (ts[i + 1] - ts[i]) / mean_dt
        # sign flips, plateaus, and strong local stretching all disqualify
        # the linear index guess that interpolation search relies on
        (0.25 <= r <= 4.0) || return KIND_BRACKET_GALLOP
    end
    return KIND_INTERPOLATION_SEARCH
end

@inline function reprobe_ts_hint!(h::TsSearchHint)
    h.kind = _ts_grid_kind(h.ts)
    h.probed_len = length(h.ts)
    return nothing
end

# Cheap growth check performed by interpolation consumers (never by the step
# loop): re-probe once the grid has doubled since the last look.
@inline function maybe_reprobe_ts_hint!(h::TsSearchHint)
    p = h.probed_len
    p == typemax(Int) && return nothing
    length(h.ts) >= 2 * max(p, 32) && reprobe_ts_hint!(h)
    return nothing
end

# Interpolation objects that carry a `TsSearchHint` return it here; the
# fallback keeps foreign interpolation types on the plain binary search. The
# typed method for `InterpolationData` lives in interp_func.jl.
@inline _ts_hint(id) = nothing

# Ending-phase re-probe, called from `_postamble!` once the time grid is
# final. Skipped when the kind was fixed from the solve options.
function _finalize_ts_hint!(sol)
    interp = hasproperty(sol, :interp) ? sol.interp : nothing
    h = _ts_hint(interp)
    h === nothing && return nothing
    h.probed_len == typemax(Int) && return nothing
    reprobe_ts_hint!(h)
    return nothing
end
