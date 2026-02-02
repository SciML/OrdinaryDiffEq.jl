function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{DP8ConstantCache, DP8Cache},
    }
    return dense ? "specialized 7th order interpolation" : "1st order linear"
end
