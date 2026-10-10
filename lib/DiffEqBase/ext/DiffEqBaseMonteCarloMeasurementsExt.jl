module DiffEqBaseMonteCarloMeasurementsExt

using DiffEqBase
import DiffEqBase: value
using MonteCarloMeasurements

# Support adaptive steps should be errorless
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{
            <:MonteCarloMeasurements.AbstractParticles,
            N,
        }, t
    ) where {N}
    return sqrt(sum(x -> abs2(value(x)), u) / max(length(u), 1))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{
            <:MonteCarloMeasurements.AbstractParticles,
            N,
        },
        t::AbstractArray{
            <:MonteCarloMeasurements.AbstractParticles,
            N,
        }
    ) where {N}
    return sqrt(sum(x -> abs2(value(x)), u) / max(length(u), 1))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::MonteCarloMeasurements.AbstractParticles, t)
    return abs(value(u))
end

end
