module StochasticDiffEqIIF

using Reexport
@reexport using StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!, initialize!, issplit,
    current_extrapolant, current_extrapolant!

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    NLSOLVEJL_SETUP, IIFNLSolveFunc, DiffCache, get_du, get_chunksize,
    determine_chunksize, unwrap_alg,
    @cache

import DiffEqBase: is_diagonal_noise, @..
import DiffEqBase: full_cache, rand_cache, ratenoise_cache

import MuladdMacro: @muladd
import SciMLBase

using LinearAlgebra
using RecursiveArrayTools

include("algorithms.jl")
include("alg_utils.jl")
include("caches/iif_caches.jl")
include("perform_step/iif.jl")

export IIF1M, IIF2M, IIF1Mil

end # module
