module StochasticDiffEqRODE

using Reexport
@reexport using StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!, initialize!, issplit

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, is_split_step,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqRODEAlgorithm, StochasticDiffEqRODEAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    @cache

import DiffEqBase: is_diagonal_noise, @..
import DiffEqBase: full_cache, rand_cache, ratenoise_cache

import MuladdMacro: @muladd

import SciMLBase

using LinearAlgebra
using StaticArrays
using RecursiveArrayTools

include("algorithms.jl")
include("alg_utils.jl")

include("caches/basic_method_caches.jl")
include("caches/dynamical_caches.jl")

include("perform_step/low_order.jl")
include("perform_step/dynamical.jl")

export RandomEM, RandomHeun, RandomTamedEM, BAOAB

end # module
