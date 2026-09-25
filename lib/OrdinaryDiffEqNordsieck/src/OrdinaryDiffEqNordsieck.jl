module OrdinaryDiffEqNordsieck

# `alg_order` is owned (and public) by SciMLBase; import it from its owner.
import SciMLBase: alg_order
# `initialize!`, `calculate_residuals`, `calculate_residuals!` are owned by
# DiffEqBase (re-provided through OrdinaryDiffEqCore); import them from the owner
# so the ExplicitImports owner check passes. `initialize!` is public in
# DiffEqBase; the `calculate_residuals` pair is DiffEqBase-internal.
import DiffEqBase: initialize!, calculate_residuals, calculate_residuals!
import OrdinaryDiffEqCore: alg_adaptive_order,
    qmin_default, qmax_default, qsteady_min_default, qsteady_max_default,
    get_current_alg_order,
    AbstractController, AbstractControllerCache, OrdinaryDiffEqAdaptiveAlgorithm,
    OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm,
    alg_cache, OrdinaryDiffEqMutableCache,
    OrdinaryDiffEqConstantCache,
    perform_step!, stepsize_controller!,
    step_accept_controller!, step_reject_controller!,
    accept_step_controller, post_newton_controller!,
    setup_controller_cache, get_qmin, get_qmax, get_qsteady_min,
    get_qsteady_max, get_EEst, set_EEst!, increment_nf!,
    CommonControllerOptions, resolve_basic, _resolved_QT,
    get_current_adaptive_order, get_fsalfirstlast,
    ode_interpolant, ode_interpolant!, trivial_limiter!,
    generic_solver_docstring, default_controller
using MuladdMacro: @muladd
using FastBroadcast: @.., Serial
using RecursiveArrayTools: recursivefill!
import LinearAlgebra: rmul!
using OrdinaryDiffEqTsit5: Tsit5ConstantCache, Tsit5Cache

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
