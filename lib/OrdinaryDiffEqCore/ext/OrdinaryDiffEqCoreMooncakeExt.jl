module OrdinaryDiffEqCoreMooncakeExt

using OrdinaryDiffEqCore, Mooncake
using Mooncake: @zero_adjoint, @is_primitive, MinimalCtx, CoDual, Dual, NoRData,
    MutableTangent, NoTangent, primal, lazy_zero_rdata, instantiate

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
# mathematically a function of `u0`. Marking it @zero_adjoint is therefore
# not strictly correct — it drops the (small) contribution of dt0 to the
# gradient. We keep it because the underlying `_ode_initdt_iip` contains
# try/catch blocks that Mooncake cannot trace through (UpsilonNode error).
# In practice the introduced error is at the level of the controller's
# first-step correction and the resulting gradient still matches ForwardDiff
# to ~10 digits in the tests below.
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.ode_determine_initdt), Vararg,
}

# Keep the SciMLBase.check_error rule as well (different from check_error!)
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.SciMLBase.check_error), Vararg,
}

# initialize_saveat computes the set of times to save the solution at from
# `saveat`/`tspan`, using `TwicePrecision`-based ranges internally (no forward-mode
# `lgetfield` rule exists for that yet, chalk-lab/Mooncake.jl -- SciML/SciMLSensitivity.jl
# #1427). Unlike `fixed_t_for_tstop_error!` above, this is safe to zero: every actual use
# of the returned heap (`isempty`, `length`, `in`, and the stepping loop's heap-top
# comparisons in solve.jl) is pure control-flow/membership logic, never a value read that
# flows into the solution or the loss. `reinit_saveat!` (used only when reinitializing an
# existing integrator, not on this code path) is intentionally left unmarked -- not
# verified here.
#
# Written by hand rather than via `@zero_derivative`: that macro's generated rule calls
# `zero_tangent`/`zero_dual` on the *actual* returned `BinaryHeap`, which recursively needs
# `zero_tangent` of its internal `valtree::Vector{Float64}` -- and building that under
# forward-mode `_new_` hits a separate, genuine Mooncake-core gap (`Vector`'s `:ref` field
# isn't marked `PossiblyUninitTangent`, so the placeholder-construction fallback tries a
# nonexistent no-arg `MemoryRef{Float64}()`). Building the zero tangent by hand with
# `zeros(Float64, n)` (plain array allocation, not Mooncake's `_new_`-based path) sidesteps
# that bug entirely.
@is_primitive MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.initialize_saveat), Type, Any, Any,
}

_zero_saveat_tangent(y) = MutableTangent((ordering=NoTangent(), valtree=zeros(Float64, length(y.valtree))))

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
    # `saveat`/`tspan` are ordinary Float64-like values, not structurally NoRData like the
    # heap's own valtree -- their zero rdata must match their *own* type (e.g. `0.0`, not
    # `NoRData()`), same reasoning as `promote_f`'s `t` argument fix above. Getting this
    # wrong crashes accumulation whenever `saveat`/`tspan` alias a genuinely-differentiated
    # value elsewhere (e.g. `saveat = p[1]`).
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
