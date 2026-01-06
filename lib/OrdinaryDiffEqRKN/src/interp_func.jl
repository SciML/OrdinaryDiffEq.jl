function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            DPRKN6ConstantCache,
            DPRKN6Cache,
        },
    }
    return dense ? "specialized 6th order interpolation" : "1st order linear"
end
