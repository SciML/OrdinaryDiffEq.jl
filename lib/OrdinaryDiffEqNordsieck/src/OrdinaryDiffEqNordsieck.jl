module OrdinaryDiffEqNordsieck

import OrdinaryDiffEqCore: alg_order, alg_adaptive_order,
    qmin_default, qmax_default, qsteady_min_default, qsteady_max_default,
    get_current_alg_order, DummyController,
    AbstractController, AbstractControllerCache, OrdinaryDiffEqAdaptiveAlgorithm,
    OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm,
    alg_cache, OrdinaryDiffEqMutableCache,
    OrdinaryDiffEqConstantCache, initialize!,
    initialize!, perform_step!, stepsize_controller!,
    step_accept_controller!, step_reject_controller!,
    accept_step_controller, post_newton_controller!,
    setup_controller_cache, get_qmin, get_qmax, get_qsteady_min,
    get_qsteady_max, get_EEst, CommonControllerOptions, resolve_basic, _resolved_QT,
    calculate_residuals, calculate_residuals!,
    get_current_adaptive_order, get_fsalfirstlast,
    ode_interpolant, ode_interpolant!, trivial_limiter!,
    generic_solver_docstring, default_controller
using MuladdMacro, FastBroadcast, RecursiveArrayTools
import LinearAlgebra: rmul!
using FastBroadcast: Serial
using OrdinaryDiffEqTsit5: Tsit5ConstantCache, Tsit5Cache
import OrdinaryDiffEqCore

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("controllers.jl")
include("alg_utils.jl")
include("nordsieck_utils.jl")
include("nordsieck_caches.jl")
include("nordsieck_perform_step.jl")

export AN5, JVODE, JVODE_Adams, JVODE_BDF

end
