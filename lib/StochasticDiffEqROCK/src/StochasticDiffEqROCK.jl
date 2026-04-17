module StochasticDiffEqROCK

using Reexport
@reexport using StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!, initialize!, issplit

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, is_split_step,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    @cache

import DiffEqBase: is_diagonal_noise, @..
import DiffEqBase: calculate_residuals, calculate_residuals!
import DiffEqBase: full_cache, rand_cache, ratenoise_cache

import MuladdMacro: @muladd
import SciMLBase

using LinearAlgebra
using StaticArrays
using RecursiveArrayTools
import DiffEqNoiseProcess
using Random: rand!

include("algorithms.jl")
include("alg_utils.jl")
include("caches/SROCK_caches.jl")
include("SROCK_tableaus.jl")
include("SROCK_utils.jl")
include("perform_step/SROCK_perform_step.jl")

export SROCK1, SROCK2, KomBurSROCK2, SROCKC2, SROCKEM, SKSROCK, TangXiaoSROCK2

end # module
