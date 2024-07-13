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
        Union{BS5ConstantCache, BS5Cache}}
    dense ? "specialized 5th order lazy interpolation" : "1st order linear"
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
        Union{DP5ConstantCache, DP5Cache}}
    dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
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