function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {cacheType <: DPRKN6Caches}
    return dense ? "specialized 6th order interpolation" : "1st order linear"
end
