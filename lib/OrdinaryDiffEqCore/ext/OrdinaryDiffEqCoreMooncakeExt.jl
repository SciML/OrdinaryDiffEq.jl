module OrdinaryDiffEqCoreMooncakeExt

using OrdinaryDiffEqCore, Mooncake
using Mooncake: @zero_adjoint, @zero_derivative, @is_primitive, MinimalCtx, CoDual, Dual,
    NoRData, MutableTangent, NoTangent, primal, lazy_zero_rdata, instantiate

# Most of these rules mirror the inactive_noinl rules in
# OrdinaryDiffEqCoreEnzymeCoreExt. They cover bookkeeping/logging/error
# checking that returns nothing tangent-bearing. The two notable exceptions
# are documented inline:
#   - fixed_t_for_tstop_error! is intentionally NOT marked here (returns
#     `ttmp` which carries the proposed new time tangent — Enzyme's
#     `inactive_noinl` lets the tangent flow through, but Mooncake's
#     @zero_adjoint zeros it).
#   - ode_determine_initdt IS marked because Mooncake cannot trace through
#     the underlying _ode_initdt_iip (try/catch UpsilonNodes), even though
#     dt0 mathematically depends on u0. The dropped contribution is small
#     and dominated by the rest of the integration.

Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.increment_nf!), Vararg,
}

# NOTE: fixed_t_for_tstop_error! is intentionally NOT marked @zero_adjoint.
# In the no-tstop branch the function returns its `ttmp` argument unchanged,
# and the caller assigns the result back to `integrator.t`. Mooncake's
# @zero_adjoint replaces the return with `zero_fcodual`, which silently
# zeros the gradient flow through `t` and (because subsequent integration
# steps and saving/interpolation depend on `t`) produces wrong gradients
# downstream. The function has no try/catch so Mooncake can trace through
# it directly. Enzyme's `inactive_noinl` has different semantics (input
# tangents are still allowed to flow through the return), so the same rule
# is correct on the Enzyme side.
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.increment_accept!), Vararg,
}
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.increment_reject!), Vararg,
}
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.check_error!), Vararg,
}
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.log_step!), Vararg,
}
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.final_progress), Vararg,
}
# NOTE: ode_determine_initdt depends on `u0`, so its return value (dt0) IS
# mathematically a function of `u0`. Marking it zero-derivative is therefore
# not strictly correct — it drops the (small) contribution of dt0 to the
# gradient. We keep it because the underlying `_ode_initdt_iip` contains
# try/catch blocks that Mooncake cannot trace through (UpsilonNode error).
# In practice the introduced error is at the level of the controller's
# first-step correction and the resulting gradient still matches ForwardDiff
# to ~10 digits in the tests below.
#
# Both modes, not just @zero_adjoint (ReverseMode only): HVP forward-differentiates
# whatever reverse rule gets built, so a ReverseMode-only primitive still lets the outer
# forward pass recurse into this function's body under HVP, hitting the same
# UpsilonNode/try-catch error one level up instead of avoiding it.
Mooncake.@zero_derivative Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.ode_determine_initdt), Vararg,
}

# Keep the SciMLBase.check_error rule as well (different from check_error!)
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.SciMLBase.check_error), Vararg,
}

# initialize_saveat builds its result from TwicePrecision ranges, and there's no
# forward-mode getfield rule for that type yet (SciML/SciMLSensitivity.jl#1427). Safe to
# zero here: the returned heap is only ever used for control flow (isempty/length/in,
# heap-top comparisons), never a value that reaches the solution or loss. reinit_saveat!
# is left unmarked since it's off this code path.
#
# Written by hand instead of via @zero_derivative: that macro's zero_tangent recurses into
# the heap's valtree::Vector{T} and fails to build an empty array's placeholder under
# forward-mode _new_. Plain zeros(eltype(y.valtree), n) sidesteps it while still matching T
# (e.g. Float32 saveat).
@is_primitive MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.initialize_saveat), Type, Any, Any,
}

_zero_saveat_tangent(y) = MutableTangent((ordering = NoTangent(), valtree = zeros(eltype(y.valtree), length(y.valtree))))

function Mooncake.frule!!(
        ::Dual{typeof(OrdinaryDiffEqCore.initialize_saveat)}, T::Dual{<:Type},
        saveat::Dual, tspan::Dual,
    )
    y = OrdinaryDiffEqCore.initialize_saveat(primal(T), primal(saveat), primal(tspan))
    return Dual(y, _zero_saveat_tangent(y))
end

function Mooncake.rrule!!(
        pf::CoDual{typeof(OrdinaryDiffEqCore.initialize_saveat)}, T::CoDual{<:Type},
        saveat::CoDual, tspan::CoDual,
    )
    y = OrdinaryDiffEqCore.initialize_saveat(primal(T), primal(saveat), primal(tspan))
    # saveat/tspan need their own-typed zero rdata (e.g. 0.0), not a hardcoded NoRData(),
    # same as promote_f's t argument -- otherwise accumulation crashes if they alias a
    # genuinely-differentiated value elsewhere (e.g. saveat = p[1]).
    lazy_pf, lazy_T, lazy_saveat, lazy_tspan = map(
        lazy_zero_rdata, (primal(pf), primal(T), primal(saveat), primal(tspan))
    )
    function initialize_saveat_pb!!(::NoRData)
        return (
            instantiate(lazy_pf), instantiate(lazy_T), instantiate(lazy_saveat),
            instantiate(lazy_tspan),
        )
    end
    return CoDual(y, _zero_saveat_tangent(y)), initialize_saveat_pb!!
end

end
