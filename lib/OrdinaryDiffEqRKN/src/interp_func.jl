function DiffEqBase.interp_summary(::Type{cacheType},
    dense::Bool) where {
    cacheType <:
    Union{DPRKN6ConstantCache,
    DPRKN6Cache}}
dense ? "specialized 6th order interpolation" : "1st order linear"
end