function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{SSPRK22, SSPRK22ConstantCache,
    SSPRK33, SSPRK33ConstantCache,
    SSPRK43, SSPRK43ConstantCache,
    SSPRK432, SSPRK432ConstantCache
}}
dense ? "2nd order \"free\" SSP interpolation" : "1st order linear"
end