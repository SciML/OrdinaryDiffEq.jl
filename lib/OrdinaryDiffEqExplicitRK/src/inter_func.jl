function SciMLBase.interp_summary(::Type{cacheType},
        dense::Bool) where {cacheType <: Union{ExplicitRKCache, ExplicitRKConstantCache}}
    dense ? "specialized RK interpolation" : "1st order linear"
end