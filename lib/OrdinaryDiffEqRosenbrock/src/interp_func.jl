function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{RosenbrockConstantCache,
        RosenbrockCache}}
    dense ? "specialized 2nd order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end
function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{RosenbrockCache, RosenbrockConstantCache}}
    dense ? "specialized 3rd order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end

function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{RosenbrockConstantCache}}
    dense ? "specialized 4rd order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end
