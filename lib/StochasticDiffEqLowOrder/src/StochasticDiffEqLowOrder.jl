module StochasticDiffEqLowOrder

using Reexport: Reexport, @reexport
@reexport using StochasticDiffEqCore
using StochasticDiffEqCore: StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!, issplit

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, is_split_step, supports_regular_jumps,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    get_Jalg, get_iterated_I, get_iterated_I!,
    IICommutative, IILevyArea,
    calc_twopoint_random!, calc_twopoint_random,
    @cache

import DiffEqBase
import DiffEqBase: calculate_residuals, calculate_residuals!
import DiffEqBase: full_cache, rand_cache, ratenoise_cache
import SciMLBase: is_diagonal_noise
import FastBroadcast: @..

import MuladdMacro: @muladd

import SciMLBase

using ArrayInterface: ArrayInterface
using LinearAlgebra: LinearAlgebra, UniformScaling, diagind, mul!, norm, rmul!
using StaticArrays: StaticArrays, SArray
using RecursiveArrayTools: RecursiveArrayTools

include("algorithms.jl")
include("alg_utils.jl")

include("caches/basic_method_caches.jl")
include("caches/lamba_caches.jl")
include("caches/predcorr_caches.jl")

include("perform_step/low_order.jl")
include("perform_step/lamba.jl")
include("perform_step/predcorr.jl")
include("perform_step/split.jl")

export EM, SplitEM, EulerHeun, LambaEM, LambaEulerHeun,
    SimplifiedEM, RKMil, RKMilCommute, PCEuler

end # module
