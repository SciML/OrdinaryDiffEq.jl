module OrdinaryDiffEqCoreMooncakeExt

using OrdinaryDiffEqCore, Mooncake
using Mooncake: @zero_adjoint, MinimalCtx
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.ode_determine_initdt),
    Any, Any, Any, Any, Any, Any, Any, Any, Any,
}
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.SciMLBase.check_error), Any,
}
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.fixed_t_for_floatingpoint_error!), Any, Any,
}
Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{
    typeof(OrdinaryDiffEqCore.final_progress), Any,
}

end
