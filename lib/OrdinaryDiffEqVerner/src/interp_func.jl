function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            Vern6Cache, Vern6ConstantCache,
        },
    }
    return dense ? "specialized 6th order lazy interpolation" : "1st order linear"
end

function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            Vern7Cache, Vern7ConstantCache,
        },
    }
    return dense ? "specialized 7th order lazy interpolation" : "1st order linear"
end

function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            Vern8Cache, Vern8ConstantCache,
        },
    }
    return dense ? "specialized 8th order lazy interpolation" : "1st order linear"
end

function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            Vern9Cache, Vern9ConstantCache,
        },
    }
    return dense ? "specialized 9th order lazy interpolation" : "1st order linear"
end

function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            RKV76IIaCache, RKV76IIaConstantCache,
        },
    }
    return dense ? "specialized 7th order lazy interpolation" : "1st order linear"
end
