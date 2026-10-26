module StochasticDiffEqROCK

using Reexport: Reexport, @reexport
@reexport using StochasticDiffEqCore
using StochasticDiffEqCore: StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, is_split_step, unwrap_alg,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    @cache

import DiffEqBase
import DiffEqBase: full_cache, rand_cache, ratenoise_cache
import SciMLBase: is_diagonal_noise
import FastBroadcast: @..

import MuladdMacro: @muladd
import SciMLBase

import Random

using LinearAlgebra: LinearAlgebra, mul!
using StaticArrays: StaticArrays, SArray
using RecursiveArrayTools: RecursiveArrayTools

include("algorithms.jl")
include("alg_utils.jl")
include("caches/SROCK_caches.jl")
include("SROCK_tableaus.jl")
include("SROCK_utils.jl")
include("perform_step/SROCK_perform_step.jl")

export SROCK1, SROCK2, KomBurSROCK2, SROCKC2, SROCKEM, SKSROCK, TangXiaoSROCK2

end # module
