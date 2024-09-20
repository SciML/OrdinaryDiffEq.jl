function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{RosenbrockCombinedConstantCache,
        RosenbrockCache}}
    dense ? "specialized $(cache.interp_order) order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end
