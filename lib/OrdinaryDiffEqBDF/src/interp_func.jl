function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            QNDFConstantCache, QNDFCache,
            FBDFConstantCache, FBDFCache,
            DFBDFConstantCache, DFBDFCache,
            MOOSE234ConstantCache, MOOSE234Cache,
        },
    }
    return dense ? "specialized backward-difference stiffness-aware interpolation" :
        "1st order linear"
end
