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
    return sqrt(
        mean(
            x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
            zip((value(x) for x in u), Iterators.repeated(t))
        )
    )
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
    return sqrt(
        mean(
            x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
            zip((value(x) for x in u), Iterators.repeated(value.(t)))
        )
    )
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::MonteCarloMeasurements.AbstractParticles, t)
    return abs(value(u))
end

end
