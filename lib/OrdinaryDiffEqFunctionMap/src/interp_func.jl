function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {cacheType <:
                            FunctionMapConstantCache}
    "left-endpoint piecewise constant"
end
function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {cacheType <: FunctionMapCache}
    "left-endpoint piecewise constant"
end