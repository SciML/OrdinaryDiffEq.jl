module OrdinaryDiffEqCoreMooncakeExt

using OrdinaryDiffEqCore, Mooncake
using Mooncake: @zero_adjoint, MinimalCtx

# These rules mirror the inactive_noinl rules in OrdinaryDiffEqCoreEnzymeCoreExt.
# These functions handle bookkeeping, logging, and error checking rather than
# numerical computations that need gradient information.

Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.increment_nf!), Vararg,
}
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.fixed_t_for_tstop_error!), Vararg,
}
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
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.ode_determine_initdt), Vararg,
}

# Keep the SciMLBase.check_error rule as well (different from check_error!)
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.SciMLBase.check_error), Vararg,
}

end
