function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            Tsit5Cache, Tsit5ConstantCache,
        },
    }
    return dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
end
