module DiffEqBaseReverseDiffExt

using DiffEqBase
import DiffEqBase: value
import ReverseDiff
import SciMLBase
import DiffEqBase.ArrayInterface

function fixed_event_time(t::ReverseDiff.TrackedReal, transform)
    return ReverseDiff.track(fixed_event_time, t, transform)
end

ReverseDiff.@grad function fixed_event_time(t, transform)
    return transform(ReverseDiff.value(t)), _ -> (nothing, nothing)
end

function implicit_event_time(t, residual::ReverseDiff.TrackedReal, slope)
    return ReverseDiff.track(implicit_event_time, t, residual, slope)
end

ReverseDiff.@grad function implicit_event_time(t, residual, slope)
    return ReverseDiff.value(t), dt -> (
            nothing, iszero(dt) ? zero(dt) : -dt / ReverseDiff.value(slope), nothing,
        )
end

function sample_event_slope(condition, t, residual, h, bounds, centered)
    lo, hi = bounds
    if centered
        left, right = clamp(t - h, lo, hi), clamp(t + h, lo, hi)
        left < right || return nothing
        return (condition(right) - condition(left)) / (right - left)
    end

    inner, outer = clamp(t + h, lo, hi), clamp(t + 2h, lo, hi)
    near, far = inner - t, outer - t
    iszero(far) && return nothing
    secant = (condition(outer) - residual) / far
    if iszero(near) || near == far
        return secant
    end
    inner_secant = (condition(inner) - residual) / near
    return inner_secant + near * (inner_secant - secant) / (far - near)
end

function converged_event_slope(condition, t, residual, bounds, ::Type{V}) where {V}
    lo, hi = bounds
    left, right = t - lo, hi - t
    h = cbrt(eps(V)) * max(one(V), abs(ReverseDiff.value(t)))
    centered = min(left, right) >= h
    if !centered
        h = if left < right
            min(h, right / 2)
        else
            -min(h, left / 2)
        end
    end

    slope = sample_event_slope(condition, t, residual, h, bounds, centered)
    slope === nothing && return nothing
    estimates = [slope]
    reltol = sqrt(eps(V))
    for _ in 1:precision(V)
        h /= 2
        estimate = sample_event_slope(condition, t, residual, h, bounds, centered)
        estimate === nothing && return nothing
        for (level, previous) in pairs(estimates)
            v = ReverseDiff.value(estimate)
            agreement = level == 1 ? reltol : 4reltol
            if isfinite(v) && !iszero(v) &&
                    abs(v - ReverseDiff.value(previous)) <= agreement * abs(v)
                return estimate
            end

            # Richardson extrapolation removes the leading stencil error.
            estimates[level] = estimate
            order = centered ? 2level : level + 1
            estimate += (estimate - previous) / (2^order - 1)
        end
        length(estimates) < 3 && push!(estimates, estimate)
    end
    return nothing
end

function DiffEqBase.find_root(
        f::F, bracket::Tuple{T, T}, rootfind::SciMLBase.RootfindOpt; kwargs...
    ) where {F, V <: AbstractFloat, T <: ReverseDiff.TrackedReal{V}}
    root = invoke(DiffEqBase.find_root, Tuple{F, Any, SciMLBase.RootfindOpt}, f, bracket, rootfind; kwargs...)
    lo, hi = minmax(bracket...)
    inner_lo, inner_hi = fixed_event_time(lo, nextfloat), fixed_event_time(hi, prevfloat)
    inner_lo < inner_hi || return root

    # Keep the root live for replay and hold the interpolation time fixed in the pullback.
    sample = fixed_event_time(clamp(root, inner_lo, inner_hi), identity)
    p = SciMLBase.NullParameters()
    condition(t) = f(t, p)
    residual = condition(sample)
    residual isa ReverseDiff.TrackedReal || return root
    slope = converged_event_slope(condition, sample, residual, (inner_lo, inner_hi), V)
    slope === nothing && return root
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
