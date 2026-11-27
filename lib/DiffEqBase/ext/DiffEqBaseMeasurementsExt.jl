module DiffEqBaseMeasurementsExt

import DiffEqBase
import SciMLBase: value
import Measurements

# Support adaptive steps should be errorless
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{
            <:Measurements.Measurement,
            N,
        },
        t
    ) where {N}
    return sqrt(sum(x -> abs2(value(x)), u) / max(length(u), 1))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Array{<:Measurements.Measurement, N},
        t
    ) where {N}
    return sqrt(sum(x -> abs2(value(x)), u) / max(length(u), 1))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::Measurements.Measurement, t)
    return abs(value(u))
end

end
