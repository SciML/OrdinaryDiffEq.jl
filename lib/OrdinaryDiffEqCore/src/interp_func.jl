abstract type OrdinaryDiffEqInterpolation{cacheType} <:
SciMLBase.AbstractDiffEqInterpolation end

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
    # Warm-start hint for the interval search in scalar `ode_interpolation`:
    # repeated evaluations (adjoints sweeping backwards, saveat post-processing,
    # delay-equation history lookups) hit nearby intervals, so galloping from
    # the previous hit beats a fresh binary search over all of `ts`. Races on
    # the hint under concurrent interpolation only degrade the starting guess,
    # never correctness.
    ts_hint::Guesser{tType}
end

# The guesser is constructed in linear-extrapolation mode unconditionally:
# `ts` is generally still being grown when the interpolation is constructed,
# so `Guesser(ts)`'s evenly-spaced probe would run over an incomplete grid.
# The linear guess only seeds the bracketing gallop — a poor guess on a very
# non-uniform grid degrades the start point, never the result — and solver
# time grids are close enough to even that it beats both a plain binary
# search and the previous-hit guess on sweeps, clustered stage times, and
# random access alike. Downstream packages construct `InterpolationData`
# positionally with these nine arguments.
function InterpolationData(
        f, timeseries, ts, ks, alg_choice, dense, cache,
        differential_vars, sensitivitymode
    )
    return InterpolationData(
        f, timeseries, ts, ks, alg_choice, dense, cache,
        differential_vars, sensitivitymode, Guesser(ts, Ref(1), true)
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

function strip_cache(cache)
    if !(cache isa OrdinaryDiffEqCore.DefaultCache)
        cache = ConstructionBase.constructorof(typeof(cache))(
            [
                nothing
                    for name in
                    fieldnames(typeof(cache))
            ]...
        )
    else
        # need to do something special for default cache
        cache = OrdinaryDiffEqCore.DefaultCache{
            Nothing, Nothing, Nothing, Nothing,
            Nothing, Nothing, Nothing, Nothing,
        }(nothing, nothing, 0, nothing)
    end

    return cache
end
