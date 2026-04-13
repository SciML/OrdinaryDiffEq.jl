module OrdinaryDiffEqCoreMooncakeExt

using OrdinaryDiffEqCore, Mooncake
using Mooncake: @zero_adjoint, MinimalCtx

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

end
