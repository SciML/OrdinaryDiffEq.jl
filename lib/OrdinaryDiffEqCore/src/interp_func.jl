abstract type OrdinaryDiffEqInterpolation{cacheType} <:
              DiffEqBase.AbstractDiffEqInterpolation end

struct InterpolationData{
    F, uType, tType, kType, algType <: Union{Nothing, Vector{Int}}, cacheType, DV} <:
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
        InterpolationData(interp.f, interp.timeseries, interp.ts, interp.ks,
            interp.alg_choice, interp.dense, interp.cache,
            interp.differential_vars, true)
    end
end

function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{
        cacheType,
}) where {
        cacheType,
}
    DiffEqBase.interp_summary(cacheType, interp.dense)
end
function DiffEqBase.interp_summary(::Type{cacheType}, dense::Bool) where {cacheType}
    dense ? "3rd order Hermite" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {cacheType <: CompositeCache}
    if !dense
        return "1st order linear"
    end
    caches = fieldtype(cacheType, :caches)
    join([DiffEqBase.interp_summary(ct, dense) for ct in fieldtypes(caches)], ", ")
end

function (interp::InterpolationData)(tvals, idxs, deriv, p, continuity::Symbol = :left)
    ode_interpolation(tvals, interp, idxs, deriv, p, continuity)
end
function (interp::InterpolationData)(val, tvals, idxs, deriv, p, continuity::Symbol = :left)
    ode_interpolation!(val, tvals, interp, idxs, deriv, p, continuity)
end

function InterpolationData(id::InterpolationData, f)
    InterpolationData(f, id.timeseries,
        id.ts,
        id.ks,
        id.alg_choice,
        id.dense,
        id.cache,
        id.differential_vars,
        id.sensitivitymode)
end

# strip interpolation of function information
function SciMLBase.strip_interpolation(id::InterpolationData)
    cache = id.cache
    
    cache = strip_cache(cache)

    InterpolationData(nothing, id.timeseries,
        id.ts,
        id.ks,
        id.alg_choice,
        id.dense,
        cache,
        id.differential_vars,
        id.sensitivitymode)
end

function strip_cache(cache::OrdinaryDiffEqCache)
    # most caches don't have jac_config or grad_config
    cache
end