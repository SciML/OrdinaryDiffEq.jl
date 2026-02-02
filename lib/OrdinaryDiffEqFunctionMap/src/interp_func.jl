function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        FunctionMapConstantCache,
    }
    return "left-endpoint piecewise constant"
end
function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {cacheType <: FunctionMapCache}
    return "left-endpoint piecewise constant"
end
