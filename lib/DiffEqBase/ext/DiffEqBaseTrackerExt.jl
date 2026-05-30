module DiffEqBaseTrackerExt

using DiffEqBase
import DiffEqBase: value
import Tracker
import RecursiveArrayTools

# Support adaptive with non-tracked time
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::Tracker.TrackedArray, t)
    return sqrt(sum(abs2, DiffEqBase.value(u)) / length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{<:Tracker.TrackedReal, N},
        t
    ) where {N}
    return sqrt(
        sum(
            x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
            zip((DiffEqBase.value(x) for x in u), Iterators.repeated(t))
        ) / length(u)
    )
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Array{<:Tracker.TrackedReal, N},
        t
    ) where {N}
    return sqrt(
        sum(
            x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
            zip((DiffEqBase.value(x) for x in u), Iterators.repeated(t))
        ) / length(u)
    )
end
@inline DiffEqBase.ODE_DEFAULT_NORM(u::Tracker.TrackedReal, t) = abs(DiffEqBase.value(u))

# Support TrackedReal time, don't drop tracking on the adaptivity there
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Tracker.TrackedArray,
        t::Tracker.TrackedReal
    )
    return sqrt(sum(abs2, u) / length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{<:Tracker.TrackedReal, N},
        t::Tracker.TrackedReal
    ) where {N}
    return sqrt(
        sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]), zip(u, Iterators.repeated(t))) /
            length(u)
    )
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Array{<:Tracker.TrackedReal, N},
        t::Tracker.TrackedReal
    ) where {N}
    return sqrt(
        sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]), zip(u, Iterators.repeated(t))) /
            length(u)
    )
end
@inline DiffEqBase.ODE_DEFAULT_NORM(u::Tracker.TrackedReal, t::Tracker.TrackedReal) = abs(u)

function DiffEqBase.solve_up(
        prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing,
        }, u0::Tracker.TrackedArray,
        p::Tracker.TrackedArray, args...; kwargs...
    )
    return Tracker.track(DiffEqBase.solve_up, prob, sensealg, u0, p, args...; kwargs...)
end

function DiffEqBase.solve_up(
        prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing,
        }, u0::Tracker.TrackedArray, p, args...;
        kwargs...
    )
    return Tracker.track(DiffEqBase.solve_up, prob, sensealg, u0, p, args...; kwargs...)
end

function DiffEqBase.solve_up(
        prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing,
        }, u0, p::Tracker.TrackedArray, args...;
        kwargs...
    )
    return Tracker.track(DiffEqBase.solve_up, prob, sensealg, u0, p, args...; kwargs...)
end

Tracker.@grad function DiffEqBase.solve_up(
        prob,
        sensealg::Union{
            Nothing,
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
        },
        u0, p, args...;
        kwargs...
    )
    sol,
        pb_f = DiffEqBase._solve_adjoint(
        prob, sensealg, Tracker.data(u0), Tracker.data(p),
        SciMLBase.TrackerOriginator(), args...; kwargs...
    )

    # In RecursiveArrayTools v4 `AbstractVectorOfArray <: AbstractArray`, so
    # the `sol isa AbstractArray` branch below now matches ODESolution and
    # returns the nested `Vector{Vector{Float64}}` in `sol.u`. Downstream
    # `sum(solve(...))` then reduces the outer vector element-wise and Tracker
    # errors with "Function output is not scalar". Return the wrapper itself
    # for AbstractVectorOfArray so the caller's reduction goes through the
    # RAT v4 AbstractArray interface and produces a scalar as before.
    if sol isa RecursiveArrayTools.AbstractVectorOfArray
        return sol, pb_f
    end
    if sol isa AbstractArray
        !hasfield(typeof(sol), :u) && return sol, pb_f # being safe here
        return sol.u, pb_f # AbstractNoTimeSolution isa AbstractArray
    end
    return sol, pb_f
end

end
