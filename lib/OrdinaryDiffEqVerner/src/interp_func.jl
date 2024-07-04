function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{Vern6Cache, Vern6ConstantCache
}}
    dense ? "specialized 6th order lazy interpolation" : "1st order linear"
end

function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{Vern7Cache, Vern7ConstantCache
}}
    dense ? "specialized 7th order lazy interpolation" : "1st order linear"
end

function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{Vern8Cache, Vern8ConstantCache
}}
    dense ? "specialized 8th order lazy interpolation" : "1st order linear"
end

function DiffEqBase.interp_summary(::Type{cacheType},
        dense::Bool) where {
        cacheType <:
        Union{Vern9Cache, Vern9ConstantCache
}}
    dense ? "specialized 9th order lazy interpolation" : "1st order linear"
end