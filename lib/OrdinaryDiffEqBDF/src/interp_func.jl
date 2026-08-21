function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            QNDFConstantCache, QNDFCache,
            FBDFConstantCache, FBDFCache,
            DFBDFConstantCache, DFBDFCache,
        },
    }
    return dense ? "specialized backward-difference stiffness-aware interpolation" :
        "1st order linear"
end

function SciMLBase.interp_summary(
        ::Type{cacheType},
        dense::Bool
    ) where {
        cacheType <:
        Union{
            NordsieckBDFConstantCache, NordsieckBDFCache,
            DNordsieckBDFConstantCache, DNordsieckBDFCache,
        },
    }
    return dense ? "Nordsieck polynomial interpolation" : "1st order linear"
end
