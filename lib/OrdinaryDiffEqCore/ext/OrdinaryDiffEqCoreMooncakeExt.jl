module OrdinaryDiffEqCoreMooncakeExt

using OrdinaryDiffEqCore, Mooncake
using Mooncake: @zero_adjoint, MinimalCtx
@zero_adjoint MinimalCtx Tuple{typeof(OrdinaryDiffEqCore.ode_determine_initdt), Any, Any, Any, Any, Any, Any, Any, Any, Any}

end