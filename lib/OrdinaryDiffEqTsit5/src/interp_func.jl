function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{Tsit5Cache, Tsit5ConstantCache
}}
    dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
end