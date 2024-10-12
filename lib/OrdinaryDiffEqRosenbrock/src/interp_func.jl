function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{Rosenbrock23ConstantCache,
        Rosenbrock32ConstantCache,
        Rosenbrock23Cache,
        Rosenbrock32Cache}}
    dense ? "specialized 2nd order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end

function DiffEqBase.interp_summary(cache::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{RosenbrockCombinedConstantCache,
        RosenbrockCache}}
    dense ? "specialized ? order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end
