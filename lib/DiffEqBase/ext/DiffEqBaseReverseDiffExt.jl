module DiffEqBaseReverseDiffExt

using DiffEqBase
import DiffEqBase: value
import ReverseDiff
import SciMLBase
import DiffEqBase.ArrayInterface

function implicit_event_time(t, residual::ReverseDiff.TrackedReal, slope)
    return ReverseDiff.track(implicit_event_time, t, residual, slope)
end

ReverseDiff.@grad function implicit_event_time(t, residual, slope)
    return t, dt -> (nothing, iszero(dt) ? zero(dt) : -dt / slope, nothing)
end

function event_condition_slope(condition, root, residual, bounds)
    lo, hi = minmax(bounds...)
    step = cbrt(eps(typeof(root))) * max(one(root), abs(root))
    left, right = root - step, root + step
    if lo <= left && right <= hi
        return (condition(right) - condition(left)) / (right - left)
    end

    # Keep a second-order stencil inside the bracket when the root is near an endpoint.
    edge = root - lo < hi - root ? hi : lo
    h = copysign(min(step, abs(edge - root) / 2), edge - root)
    inner = root + h
    outer = clamp(root + 2h, lo, hi)
    near, far = inner - root, outer - root
    secant = (condition(outer) - residual) / far
    if iszero(near) || near == far
        return secant
    end
    inner_secant = (condition(inner) - residual) / near
    return inner_secant + near * (inner_secant - secant) / (far - near)
end

function DiffEqBase.find_root(
        f::F, bracket::Tuple{T, T}, rootfind::SciMLBase.RootfindOpt; kwargs...
    ) where {F, T <: ReverseDiff.TrackedReal}
    bounds = map(ReverseDiff.value, bracket)
    primal_condition(t, p) = ReverseDiff.value(f(t, p))
    root = DiffEqBase.find_root(primal_condition, bounds, rootfind; kwargs...)
    p = SciMLBase.NullParameters()
    residual = f(root, p)
    residual isa ReverseDiff.TrackedReal || return root

    # At a regular crossing, the event time derivative is -∂ₚf/∂ₜf.
    slope = event_condition_slope(
        t -> primal_condition(t, p), root, ReverseDiff.value(residual), bounds
    )
    return implicit_event_time(root, residual, slope)
end

# Support adaptive with non-tracked time
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::ReverseDiff.TrackedArray, t)
    return sqrt(sum(abs2, DiffEqBase.value(u)) / length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{<:ReverseDiff.TrackedReal, N},
        t
    ) where {N}
    return sqrt(sum(x -> abs2(DiffEqBase.value(x)), u) / max(length(u), 1))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Array{<:ReverseDiff.TrackedReal, N},
        t
    ) where {N}
    return sqrt(sum(x -> abs2(DiffEqBase.value(x)), u) / max(length(u), 1))
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
    return sqrt(sum(abs2, u) / max(length(u), 1))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::Array{<:ReverseDiff.TrackedReal, N},
        t::ReverseDiff.TrackedReal
    ) where {N}
    return sqrt(sum(abs2, u) / max(length(u), 1))
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
        prob::DiffEqBase.AbstractDEProblem,
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
        prob::DiffEqBase.AbstractDEProblem,
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

# Required because ReverseDiff.@grad function DiffEqBase.solve_up is not supported!
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
