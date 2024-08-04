function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{DP8ConstantCache, DP8Cache}}
    dense ? "specialized 7th order interpolation" : "1st order linear"
end