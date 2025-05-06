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
function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{RosenbrockCombinedConstantCache, Rodas23WConstantCache, Rodas3PConstantCache,
        RosenbrockCache, Rodas23WCache, Rodas3PCache}}
    dense ? "specialized 3rd order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end

function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{RosenbrockCombinedConstantCache,
        RosenbrockCache}}
    dense ? "specialized 4rd order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end
