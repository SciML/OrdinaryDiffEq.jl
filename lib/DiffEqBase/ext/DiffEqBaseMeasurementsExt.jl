module DiffEqBaseMeasurementsExt

using DiffEqBase
import DiffEqBase: value
using Measurements

# Strip the ± uncertainty when DiffEqBase.value is called on a Measurement.
# Required so generic algorithms (e.g. CVHin initdt) can convert intermediate
# Measurement-valued quantities back to plain `_tType` step sizes.
value(x::Measurements.Measurement) = Measurements.value(x)
value(::Type{Measurements.Measurement{T}}) where {T} = T

# Support adaptive steps should be errorless
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{
            <:Measurements.Measurement,
            N,
        },
        t
    ) where {N}
    return sqrt(
        sum(
            x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
            zip((value(x) for x in u), Iterators.repeated(t))
        ) / length(u)
    )
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Array{<:Measurements.Measurement, N},
        t
    ) where {N}
    return sqrt(
        sum(
            x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
            zip((value(x) for x in u), Iterators.repeated(t))
        ) / length(u)
    )
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::Measurements.Measurement, t)
    return abs(Measurements.value(u))
end

end
