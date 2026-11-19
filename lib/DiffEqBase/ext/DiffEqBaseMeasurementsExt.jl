module DiffEqBaseMeasurementsExt

using DiffEqBase: DiffEqBase
using Measurements: Measurements

# Support adaptive steps should be errorless
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{
            <:Measurements.Measurement,
            N,
        },
        t
    ) where {N}
    return sqrt(sum(x -> abs2(x.val), u) / max(length(u), 1))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Array{<:Measurements.Measurement, N},
        t
    ) where {N}
    return sqrt(sum(x -> abs2(x.val), u) / max(length(u), 1))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::Measurements.Measurement, t)
    return abs(u.val)
end

end
