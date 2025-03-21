function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{ExplicitTaylor2Cache, ExplicitTaylor2ConstantCache, DAETSCache
}}
    dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
end
