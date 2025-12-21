abstract type OrdinaryDiffEqInterpolation{cacheType} <:
              SciMLBase.AbstractDiffEqInterpolation end

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

function SciMLBase.interp_summary(interp::OrdinaryDiffEqInterpolation{
        cacheType,
}) where {
        cacheType,
}
    SciMLBase.interp_summary(cacheType, interp.dense)
end
function SciMLBase.interp_summary(::Type{cacheType}, dense::Bool) where {cacheType}
    dense ? "3rd order Hermite" : "1st order linear"
end
function SciMLBase.interp_summary(::Type{cacheType},
        dense::Bool) where {cacheType <: CompositeCache}
    if !dense
        return "1st order linear"
    end
    caches = fieldtype(cacheType, :caches)
    join([SciMLBase.interp_summary(ct, dense) for ct in fieldtypes(caches)], ", ")
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
    cache = strip_cache(id.cache)

    InterpolationData(nothing, id.timeseries,
        id.ts,
        id.ks,
        id.alg_choice,
        id.dense,
        cache,
        id.differential_vars,
        id.sensitivitymode)
end

function strip_cache(cache)
    if cache isa OrdinaryDiffEqCore.DefaultCache
        cache = OrdinaryDiffEqCore.DefaultCache{Nothing, Nothing, Nothing, Nothing,
            Nothing, Nothing, Nothing, Nothing}(nothing, nothing, 0, nothing)
    else
        try
            cache = SciMLBase.constructorof(typeof(cache))([nothing
                                                            for name in
                                                                fieldnames(typeof(cache))]...)
        catch
            cache = (; (name => nothing for name in fieldnames(typeof(cache)))...)
        end
    end

    cache
end
