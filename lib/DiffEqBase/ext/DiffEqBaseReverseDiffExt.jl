module DiffEqBaseReverseDiffExt

using DiffEqBase
import DiffEqBase: value
import ReverseDiff
import DiffEqBase.ArrayInterface

# Support adaptive with non-tracked time
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::ReverseDiff.TrackedArray, t)
    return sqrt(sum(abs2, DiffEqBase.value(u)) / length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{<:ReverseDiff.TrackedReal, N},
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
        u::Array{<:ReverseDiff.TrackedReal, N},
        t
    ) where {N}
    return sqrt(
        sum(
            x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
            zip((DiffEqBase.value(x) for x in u), Iterators.repeated(t))
        ) / length(u)
    )
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::ReverseDiff.TrackedReal, t)
    return abs(DiffEqBase.value(u))
end

# Support TrackedReal time, don't drop tracking on the adaptivity there
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::ReverseDiff.TrackedArray,
        t::ReverseDiff.TrackedReal
    )
    return sqrt(sum(abs2, u) / length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{<:ReverseDiff.TrackedReal, N},
        t::ReverseDiff.TrackedReal
    ) where {N}
    return sqrt(
        sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]), zip(u, Iterators.repeated(t))) /
            length(u)
    )
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Array{<:ReverseDiff.TrackedReal, N},
        t::ReverseDiff.TrackedReal
    ) where {N}
    return sqrt(
        sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]), zip(u, Iterators.repeated(t))) /
            length(u)
    )
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::ReverseDiff.TrackedReal,
        t::ReverseDiff.TrackedReal
    )
    return abs(u)
end

# `ReverseDiff.TrackedArray`
function DiffEqBase.solve_up(
        prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing,
        }, u0::ReverseDiff.TrackedArray,
        p::ReverseDiff.TrackedArray, args...; kwargs...
    )
    return ReverseDiff.track(DiffEqBase.solve_up, prob, sensealg, u0, p, args...; kwargs...)
end

function DiffEqBase.solve_up(
        prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing,
        }, u0, p::ReverseDiff.TrackedArray,
        args...; kwargs...
    )
    return ReverseDiff.track(DiffEqBase.solve_up, prob, sensealg, u0, p, args...; kwargs...)
end

function DiffEqBase.solve_up(
        prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing,
        }, u0::ReverseDiff.TrackedArray, p,
        args...; kwargs...
    )
    return ReverseDiff.track(DiffEqBase.solve_up, prob, sensealg, u0, p, args...; kwargs...)
end

# `AbstractArray{<:ReverseDiff.TrackedReal}`
function DiffEqBase.solve_up(
        prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing,
        },
        u0::AbstractArray{<:ReverseDiff.TrackedReal},
        p::AbstractArray{<:ReverseDiff.TrackedReal}, args...;
        kwargs...
    )
    return DiffEqBase.solve_up(
        prob, sensealg, ArrayInterface.aos_to_soa(u0),
        ArrayInterface.aos_to_soa(p), args...;
        kwargs...
    )
end

function DiffEqBase.solve_up(
        prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing,
        }, u0,
        p::AbstractArray{<:ReverseDiff.TrackedReal},
        args...; kwargs...
    )
    return DiffEqBase.solve_up(
        prob, sensealg, u0, ArrayInterface.aos_to_soa(p), args...; kwargs...
    )
end

function DiffEqBase.solve_up(
        prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing,
        }, u0::ReverseDiff.TrackedArray,
        p::AbstractArray{<:ReverseDiff.TrackedReal},
        args...; kwargs...
    )
    return DiffEqBase.solve_up(
        prob, sensealg, u0, ArrayInterface.aos_to_soa(p), args...; kwargs...
    )
end

function DiffEqBase.solve_up(
        prob::DiffEqBase.DEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing,
        },
        u0::AbstractArray{<:ReverseDiff.TrackedReal}, p,
        args...; kwargs...
    )
    return DiffEqBase.solve_up(
        prob, sensealg, ArrayInterface.aos_to_soa(u0), p, args...; kwargs...
    )
end

function DiffEqBase.solve_up(
        prob::DiffEqBase.DEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing,
        },
        u0::AbstractArray{<:ReverseDiff.TrackedReal}, p::ReverseDiff.TrackedArray,
        args...; kwargs...
    )
    return DiffEqBase.solve_up(
        prob, sensealg, ArrayInterface.aos_to_soa(u0), p, args...; kwargs...
    )
end

# Required becase ReverseDiff.@grad function DiffEqBase.solve_up is not supported!
import DiffEqBase: solve_up
ReverseDiff.@grad function solve_up(prob, sensealg, u0, p, args...; kwargs...)
    out = DiffEqBase._solve_adjoint(
        prob, sensealg, ReverseDiff.value(u0),
        ReverseDiff.value(p),
        SciMLBase.ReverseDiffOriginator(), args...; kwargs...
    )
    function actual_adjoint(_args...)
        original_adjoint = out[2](_args...)
        if isempty(args) # alg is missing
            tuple(original_adjoint[1:4]..., original_adjoint[6:end]...)
        else
            original_adjoint
        end
    end
    Array(out[1]), actual_adjoint
end

end
