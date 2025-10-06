function SciMLBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{Rosenbrock23ConstantCache,
        Rosenbrock32ConstantCache,
        RosenbrockCombinedCache}}
    dense ? "specialized 2nd order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end
function SciMLBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{RosenbrockCombinedConstantCache, Rodas23WConstantCache, Rodas3PConstantCache,
        RosenbrockCache, Rodas23WCache, Rodas3PCache}}
    dense ? "specialized 3rd order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end

function SciMLBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{RosenbrockCombinedConstantCache,
        RosenbrockCache}}
    dense ? "specialized 4th (Rodas6P = 5th) order \"free\" stiffness-aware interpolation" :
    "1st order linear"
end
