module OrdinaryDiffEqPDIRK

import OrdinaryDiffEqCore: isfsal, alg_order, _unwrap_val,
                           OrdinaryDiffEqNewtonAlgorithm, OrdinaryDiffEqConstantCache,
                           OrdinaryDiffEqMutableCache, constvalue, alg_cache,
                           uses_uprev, @unpack, unwrap_alg, @cache, DEFAULT_PRECS,
                           @threaded, initialize!, perform_step!, isthreaded,
                           full_cache, get_fsalfirstlast, differentiation_rk_docstring,
                           bool_to_ADType
import StaticArrays: SVector
import MuladdMacro: @muladd
import FastBroadcast: @..
using Polyester

using Reexport
@reexport using DiffEqBase

using OrdinaryDiffEqDifferentiation: dolinsolve
using OrdinaryDiffEqNonlinearSolve: NLNewton, build_nlsolver, nlsolve!, nlsolvefail,
                                    markfirststage!

import ADTypes: AutoForwardDiff, AbstractADType

include("algorithms.jl")
include("alg_utils.jl")
include("pdirk_caches.jl")
include("pdirk_perform_step.jl")

export PDIRK44

end
