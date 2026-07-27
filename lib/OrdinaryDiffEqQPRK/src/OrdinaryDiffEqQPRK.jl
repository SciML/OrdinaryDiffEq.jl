module OrdinaryDiffEqQPRK

import OrdinaryDiffEqCore: TmpCache, build_tmp_cache,
    OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqConstantCache,
    explicit_rk_docstring, @cache,
    OrdinaryDiffEqMutableCache,
    @fold, @OnDemandTableauExtract,
    trivial_limiter!, alg_cache,
    perform_step!, get_fsalfirstlast,
    constvalue
import SciMLBase: alg_order
import DiffEqBase: initialize!, calculate_residuals!, calculate_residuals
using FastBroadcast: FastBroadcast, Serial, @..
using MuladdMacro: MuladdMacro, @muladd
using RecursiveArrayTools: recursive_unitless_bottom_eltype, recursivefill!
import OrdinaryDiffEqCore

using SciMLBase: SciMLBase
using Reexport: Reexport, @reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("qprk_caches.jl")
include("qprk_tableaus.jl")
include("qprk_perform_step.jl")

export QPRK98

end
