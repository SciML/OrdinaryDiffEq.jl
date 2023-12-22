abstract type OrdinaryDiffEqInterpolation{cacheType} <:
              DiffEqBase.AbstractDiffEqInterpolation end

struct InterpolationData{F, uType, tType, kType, cacheType} <:
       OrdinaryDiffEqInterpolation{cacheType}
    f::F
    timeseries::uType
    ts::tType
    ks::kType
    dense::Bool
    cache::cacheType
end

struct CompositeInterpolationData{F, uType, tType, kType, cacheType} <:
       OrdinaryDiffEqInterpolation{cacheType}
    f::F
    timeseries::uType
    ts::tType
    ks::kType
    alg_choice::Vector{Int}
    dense::Bool
    cache::cacheType
end

function DiffEqBase.interp_summary(interp::OrdinaryDiffEqInterpolation{
    cacheType,
}) where {
    cacheType,
}
    DiffEqBase.interp_summary(cacheType, interp.dense)
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {cacheType <:
         FunctionMapConstantCache}
    "left-endpoint piecewise constant"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {cacheType <: FunctionMapCache}
    "left-endpoint piecewise constant"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{DP5ConstantCache, DP5Cache}}
    dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{Rosenbrock23ConstantCache,
        Rosenbrock32ConstantCache,
        Rosenbrock23Cache,
        Rosenbrock32Cache}}
    dense ? "specialized 2nd order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{Rodas4ConstantCache, Rodas23WConstantCache, Rodas3PConstantCache,
        Rodas4Cache, Rodas23WCache, Rodas3PCache}}
    dense ? "specialized 3rd order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{Rosenbrock5ConstantCache,
        Rosenbrock5Cache}}
    dense ? "specialized 4rd order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{SSPRK22, SSPRK22ConstantCache,
        SSPRK33, SSPRK33ConstantCache,
        SSPRK43, SSPRK43ConstantCache,
        SSPRK432, SSPRK432ConstantCache,
    }}
    dense ? "2nd order \"free\" SSP interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <: Union{OwrenZen3Cache,
        OwrenZen3ConstantCache}}
    dense ? "specialized 3rd order \"free\" interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <: Union{OwrenZen4Cache,
        OwrenZen4ConstantCache}}
    dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <: Union{OwrenZen5Cache,
        OwrenZen5ConstantCache}}
    dense ? "specialized 5th order \"free\" interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{Tsit5Cache, Tsit5ConstantCache
    }}
    dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{BS5ConstantCache, BS5Cache}}
    dense ? "specialized 5th order lazy interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{Vern6Cache, Vern6ConstantCache
    }}
    dense ? "specialized 6th order lazy interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{Vern7Cache, Vern7ConstantCache
    }}
    dense ? "specialized 7th order lazy interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{Vern8Cache, Vern8ConstantCache
    }}
    dense ? "specialized 8th order lazy interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{Vern9Cache, Vern9ConstantCache
    }}
    dense ? "specialized 9th order lazy interpolation" : "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{DP8ConstantCache, DP8Cache}}
    dense ? "specialized 7th order interpolation" : "1st order linear"
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
function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{DPRKN6ConstantCache,
        DPRKN6Cache}}
    dense ? "specialized 6th order interpolation" : "1st order linear"
end

function (interp::InterpolationData)(tvals, idxs, deriv, p, continuity::Symbol = :left)
    ode_interpolation(tvals, interp, idxs, deriv, p, continuity)
end
function (interp::CompositeInterpolationData)(tvals, idxs, deriv, p,
    continuity::Symbol = :left)
    ode_interpolation(tvals, interp, idxs, deriv, p, continuity)
end
function (interp::InterpolationData)(val, tvals, idxs, deriv, p, continuity::Symbol = :left)
    ode_interpolation!(val, tvals, interp, idxs, deriv, p, continuity)
end
function (interp::CompositeInterpolationData)(val, tvals, idxs, deriv, p,
    continuity::Symbol = :left)
    ode_interpolation!(val, tvals, interp, idxs, deriv, p, continuity)
end

function InterpolationData(id::InterpolationData, f)
    InterpolationData(f, id.timeseries,
        id.ts,
        id.ks,
        id.dense,
        id.cache)
end

function CompositeInterpolationData(id::CompositeInterpolationData, f)
    CompositeInterpolationData(f, id.timeseries,
        id.ts,
        id.ks,
        id.alg_choice,
        id.dense,
        id.cache)
end
