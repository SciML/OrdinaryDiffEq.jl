module OrdinaryDiffEqStabilizedRK

import OrdinaryDiffEqCore: alg_adaptive_order,
    gamma_default, qmax_default, alg_extrapolates,
    fac_default_gamma, has_dtnew_modification,
    perform_step!, unwrap_alg,
    default_controller, PredictiveController,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqAdaptiveAlgorithm,
    alg_cache, @cache,
    constvalue, full_cache, get_fsalfirstlast,
    generic_solver_docstring, trivial_limiter!
import OrdinaryDiffEqCore
using FastBroadcast: FastBroadcast, @..
using MuladdMacro: MuladdMacro, @muladd
using RecursiveArrayTools: RecursiveArrayTools, recursivefill!
# `calculate_residuals`/`calculate_residuals!`/`initialize!` are owned by DiffEqBase;
# `alg_order`/`_vec`/`value` are owned by SciMLBase. Import from the owner.
# `initialize!`/`alg_order` are extended here, so they must be `import`ed.
import DiffEqBase: initialize!
using DiffEqBase: calculate_residuals, calculate_residuals!

using SciMLBase: SciMLBase
import SciMLBase: alg_order, _vec, value
using Reexport: Reexport, @reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("rkc_utils.jl")
include("rkc_caches.jl")
include("rkc_perform_step.jl")
include("rkc_tableaus_serk2.jl")
include("rkc_tableaus_rock4.jl")
include("rkc_tableaus_rock2.jl")
include("rkc_tableaus_eserk5.jl")
include("rkc_tableaus_eserk4.jl")

export ROCK2, ROCK4, RKC, RKMC2, ESERK4, ESERK5, SERK2, TSRKC2, TSRKC3, RKL1, RKL2, RKG1, RKG2

end
