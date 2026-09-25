module StochasticDiffEqIIF

using Reexport: Reexport, @reexport
@reexport using StochasticDiffEqCore
using StochasticDiffEqCore: StochasticDiffEqCore

import OrdinaryDiffEqCore
# `perform_step!`, `issplit`, `current_extrapolant`, `current_extrapolant!` are part
# of OrdinaryDiffEqCore's solver-author interface but are not (yet) declared `public`,
# so they are tightly ignored in the QA explicit-imports checks. (`initialize!` is
# owned by DiffEqBase and imported from its public owner there.)
import OrdinaryDiffEqCore: perform_step!, issplit,
    current_extrapolant, current_extrapolant!

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    NLSOLVEJL_SETUP, IIFNLSolveFunc, DiffCache, get_du, get_chunksize,
    determine_chunksize, unwrap_alg,
    @cache

import DiffEqBase
import DiffEqBase: initialize!, full_cache, rand_cache, ratenoise_cache
import SciMLBase: is_diagonal_noise
import FastBroadcast: @..

import MuladdMacro: @muladd
import SciMLBase

using LinearAlgebra: mul!, rmul!

include("algorithms.jl")
include("alg_utils.jl")
include("caches/iif_caches.jl")
include("perform_step/iif.jl")

export IIF1M, IIF2M, IIF1Mil

end # module
