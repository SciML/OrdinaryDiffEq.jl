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
end

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

Return a lightweight copy of `cache` with all fields set to `nothing`, used by
`SciMLBase.strip_interpolation` to drop the (potentially large) working buffers
from a solution's interpolation object before serialization. Has a special path
for [`DefaultCache`](@ref).
"""
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
