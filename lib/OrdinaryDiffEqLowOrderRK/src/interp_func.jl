function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <: Union{
            OwrenZen3Cache,
            OwrenZen3ConstantCache,
        },
    }
    return dense ? "specialized 3rd order \"free\" interpolation" : "1st order linear"
end
function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <: Union{
            OwrenZen4Cache,
            OwrenZen4ConstantCache,
        },
    }
    return dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
end
function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <: Union{
            OwrenZen5Cache,
            OwrenZen5ConstantCache,
        },
    }
    return dense ? "specialized 5th order \"free\" interpolation" : "1st order linear"
end
function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{DP5ConstantCache, DP5Cache},
    }
    return dense ? "specialized 4th order \"free\" interpolation" : "1st order linear"
end
function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{BS5ConstantCache, BS5Cache},
    }
    return dense ? "specialized 5th order lazy interpolation" : "1st order linear"
end
