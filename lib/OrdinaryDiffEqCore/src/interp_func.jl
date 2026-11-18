"""
    OrdinaryDiffEqInterpolation{cacheType} <: SciMLBase.AbstractDiffEqInterpolation

Abstract supertype for the dense-output interpolation object attached to a
solution. Given a saved timeseries plus derivative (`k`) history it evaluates the
continuous extension. See [`InterpolationData`](@ref) for the concrete type.
"""
abstract type OrdinaryDiffEqInterpolation{cacheType} <:
SciMLBase.AbstractDiffEqInterpolation end

"""
    InterpolationData(f, timeseries, ts, ks, alg_choice, dense, cache, differential_vars, sensitivitymode)

Concrete [`OrdinaryDiffEqInterpolation`](@ref) storing everything needed to
evaluate the continuous solution: the RHS `f`, the saved states `timeseries` at
times `ts`, the stage-derivative history `ks`, the per-step `alg_choice` (for
composite algorithms), whether `dense` output is available, the solver `cache`,
the `differential_vars` mask (for DAEs), and a `sensitivitymode` flag. Calling
`(interp)(tvals, idxs, deriv, p, continuity)` performs the interpolation.
"""
struct InterpolationData{
        F, uType, tType, kType, algType <: Union{Nothing, Vector{Int}}, cacheType, DV,
    } <:
    OrdinaryDiffEqInterpolation{cacheType}
    f::F
    timeseries::uType
    ts::tType
    ks::kType
    alg_choice::algType
    dense::Bool
    cache::cacheType
    differential_vars::DV
    sensitivitymode::Bool
    # Warm-start hint for the interval search in scalar `ode_interpolation`;
    # see `TsSearchHint`.
    ts_hint::TsSearchHint{tType}
end

# Downstream packages (e.g. StochasticDiffEq) construct `InterpolationData`
# positionally with these nine arguments; the hint then starts on the robust
# gallop strategy and re-selects itself as the grid is probed.
function InterpolationData(
        f, timeseries, ts, ks, alg_choice, dense, cache,
        differential_vars, sensitivitymode
    )
    return InterpolationData(
        f, timeseries, ts, ks, alg_choice, dense, cache,
        differential_vars, sensitivitymode, TsSearchHint(ts)
    )
end

@inline _ts_hint(id::InterpolationData) = id.ts_hint

@static if isdefined(SciMLBase, :enable_interpolation_sensitivitymode)
    function SciMLBase.enable_interpolation_sensitivitymode(interp::InterpolationData)
        InterpolationData(
            interp.f, interp.timeseries, interp.ts, interp.ks,
            interp.alg_choice, interp.dense, interp.cache,
            interp.differential_vars, true
        )
    end
end

function SciMLBase.interp_summary(
        interp::OrdinaryDiffEqInterpolation{
            cacheType,
        }
    ) where {
        cacheType,
    }
    return SciMLBase.interp_summary(cacheType, interp.dense)
end
function SciMLBase.interp_summary(::Type{cacheType}, dense::Bool) where {cacheType}
    return dense ? "3rd order Hermite" : "1st order linear"
end
function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {cacheType <: CompositeCache}
    if !dense
        return "1st order linear"
    end
    caches = fieldtype(cacheType, :caches)
    return join([SciMLBase.interp_summary(ct, dense) for ct in fieldtypes(caches)], ", ")
end

function (interp::InterpolationData)(tvals, idxs, deriv, p, continuity::Symbol = :left)
    return ode_interpolation(tvals, interp, idxs, deriv, p, continuity)
end
function (interp::InterpolationData)(val, tvals, idxs, deriv, p, continuity::Symbol = :left)
    return ode_interpolation!(val, tvals, interp, idxs, deriv, p, continuity)
end

function InterpolationData(id::InterpolationData, f)
    return InterpolationData(
        f, id.timeseries,
        id.ts,
        id.ks,
        id.alg_choice,
        id.dense,
        id.cache,
        id.differential_vars,
        id.sensitivitymode
    )
end

# strip interpolation of function information
function SciMLBase.strip_interpolation(id::InterpolationData)
    cache = strip_cache(id.cache)

    return InterpolationData(
        nothing, id.timeseries,
        id.ts,
        id.ks,
        id.alg_choice,
        id.dense,
        cache,
        id.differential_vars,
        id.sensitivitymode
    )
end

"""
    strip_cache(cache)

Return a lightweight copy of `cache` for `SciMLBase.strip_interpolation`, which
drops function information and working buffers from a solution's interpolation
object before serialization. Fields are set to `nothing` except the tableau
(`tab`, `tabs`), numbers and fields whose type cannot hold `nothing`,
since the interpolant reads these; kept arrays of non-numbers are emptied.
[`DefaultCache`](@ref) and `CompositeCache` strip each of their sub-caches.

# Developer API

Solver extensions may specialize this hook for their cache types. End-user code
should call `solve` and use solution APIs, rather than construct caches or
depend on cache fields and their serialized representation.
"""
function strip_cache(cache)
    T = typeof(cache)
    declared = _declared_fieldtypes(T)
    eltypevars = map(i -> isdefined(cache, i) ? _nothing_eltype(declared[i], getfield(cache, i)) : nothing, 1:fieldcount(T))
    keep = map(1:fieldcount(T)) do i
        isdefined(cache, i) || return !_accepts_nothing(declared[i])
        value = getfield(cache, i)
        return fieldname(T, i) in (:tab, :tabs) || value isa Number ||
            !(_accepts_nothing(declared[i]) || eltypevars[i] !== nothing)
    end
    while true
        pinned = findfirst(eachindex(keep)) do j
            var = declared[j] isa TypeVar ? declared[j] : eltypevars[j]
            return !keep[j] && var !== nothing &&
                any(i -> keep[i] && _mentions(declared[i], var), eachindex(keep))
        end
        pinned === nothing && break
        keep[pinned] = true
    end
    fields = map(1:fieldcount(T)) do i
        keep[i] && isdefined(cache, i) && return _strip_kept(fieldname(T, i), getfield(cache, i))
        eltypevars[i] === nothing && return nothing
        return Array{Nothing}(undef, map(zero, size(getfield(cache, i))))
    end
    return ConstructionBase.constructorof(T)(fields...)
end

# A field declared `Array{E}` holding non-numbers, with `E` a parameter that admits `nothing`,
# can become an empty `Array{Nothing}` so that `E` no longer names user function types.
_nothing_eltype(declared, value) = value isa Array && !(eltype(value) <: Number) &&
    declared isa DataType && declared.name.wrapper === Array &&
    declared.parameters[1] isa TypeVar && _accepts_nothing(declared.parameters[1]) ?
    declared.parameters[1] : nothing

# Arrays the interpolant does not read, such as per-stage linear solvers, are emptied.
_strip_kept(name, value) = name in (:tab, :tabs) || !(value isa Array) || eltype(value) <: Number ?
    value : similar(value, ntuple(_ -> 0, ndims(value)))

function _declared_fieldtypes(T)
    wrapper = T.name.wrapper
    while wrapper isa UnionAll
        wrapper = wrapper.body
    end
    return fieldtypes(wrapper)
end

_accepts_nothing(declared::TypeVar) = Nothing <: declared.ub
_accepts_nothing(declared) = Nothing <: declared

_mentions(declared, var) = declared === var ||
    (declared isa DataType && any(p -> _mentions(p, var), declared.parameters)) ||
    (declared isa Union && (_mentions(declared.a, var) || _mentions(declared.b, var)))

function strip_cache(cache::DefaultCache)
    names = (:cache1, :cache2, :cache3, :cache4, :cache5, :cache6)
    caches = map(n -> isdefined(cache, n) ? strip_cache(getfield(cache, n)) : nothing, names)
    types = map((c, n) -> isdefined(cache, n) ? typeof(c) : Nothing, caches, names)
    stripped = DefaultCache{types..., Nothing, Nothing}(nothing, nothing, cache.current, nothing)
    foreach((c, n) -> isdefined(cache, n) && setfield!(stripped, n, c), caches, names)
    return stripped
end

function strip_cache(cache::CompositeCache)
    return CompositeCache(map(strip_cache, cache.caches), nothing, cache.current)
end
