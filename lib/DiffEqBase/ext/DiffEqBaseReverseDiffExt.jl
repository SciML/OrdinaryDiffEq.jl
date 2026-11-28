module DiffEqBaseReverseDiffExt

using DiffEqBase
import DiffEqBase: value
import ReverseDiff, ForwardDiff
import SciMLBase
import DiffEqBase.ArrayInterface

function fixed_event_time(t::ReverseDiff.TrackedReal)
    return ReverseDiff.track(fixed_event_time, t)
end

ReverseDiff.@grad function fixed_event_time(t)
    return ReverseDiff.value(t), dt -> (zero(dt),)
end

function implicit_event_time(t, residual::ReverseDiff.TrackedReal, slope)
    return ReverseDiff.track(implicit_event_time, t, residual, slope)
end

ReverseDiff.@grad function implicit_event_time(t, residual, slope)
    function pullback(dt)
        s = ReverseDiff.value(slope)
        if !isfinite(s) || iszero(s)
            return dt, zero(dt), zero(dt)
        end
        return zero(dt), -dt / s, zero(dt)
    end
    return ReverseDiff.value(t), pullback
end

struct EventTimeDerivative end
const EventTimeTag = ForwardDiff.Tag{EventTimeDerivative}

function dual_callback_state(integrator, callback, t::ForwardDiff.Dual{Tag}) where {Tag}
    time = ForwardDiff.value(t)
    kwargs = callback.idxs === nothing ? (;) : (; idxs = callback.idxs)
    state = integrator(time, Val{0}; kwargs...)
    tangent = integrator(time, Val{1}; kwargs...)
    primal_with_tangent(u, du) = ForwardDiff.Dual{Tag}(u, du)
    return map(primal_with_tangent, state, tangent)
end

function DiffEqBase.get_condition(
        integrator::SciMLBase.DEIntegrator, callback::SciMLBase.ContinuousCallback,
        t::ForwardDiff.Dual{<:EventTimeTag}
    )
    state = dual_callback_state(integrator, callback, t)
    integrator.sol.stats.ncondition += 1
    return callback.condition(state, t, integrator)
end

function DiffEqBase.get_condition(
        integrator::SciMLBase.DEIntegrator, callback::SciMLBase.VectorContinuousCallback,
        t::ForwardDiff.Dual{<:EventTimeTag}
    )
    state = dual_callback_state(integrator, callback, t)
    result = similar(integrator.callback_cache.tmp_condition, typeof(t), callback.len)
    integrator.sol.stats.ncondition += 1
    callback.condition(result, state, t, integrator)
    return result
end

function DiffEqBase.find_root(
        f::F, bracket::Tuple{T, T}, rootfind::SciMLBase.RootfindOpt; kwargs...
    ) where {F, V <: AbstractFloat, T <: ReverseDiff.TrackedReal{V}}
    root = invoke(DiffEqBase.find_root, Tuple{F, Any, SciMLBase.RootfindOpt}, f, bracket, rootfind; kwargs...)
    # The residual sees a fixed time in ReverseDiff and a unit time tangent in ForwardDiff.
    sample = fixed_event_time(root)
    tag = typeof(ForwardDiff.Tag(EventTimeDerivative(), typeof(sample)))
    time = ForwardDiff.Dual{tag}(sample, one(V))
    p = SciMLBase.NullParameters()
    applicable(f, time, p) || return root
    result = f(time, p)
    result isa ForwardDiff.Dual{tag} || return root
    residual = ForwardDiff.value(result)
    residual isa ReverseDiff.TrackedReal || return root
    slope = ForwardDiff.partials(result, 1)
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
