module StochasticDiffEqMilstein

using Reexport: Reexport, @reexport
@reexport using StochasticDiffEqCore
using StochasticDiffEqCore: StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, is_split_step,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    get_Jalg, get_iterated_I, get_iterated_I!,
    IICommutative, IILevyArea,
    _z_prototype,
    @cache

import DiffEqBase
import DiffEqBase: calculate_residuals, calculate_residuals!
import DiffEqBase: full_cache, rand_cache, ratenoise_cache
import SciMLBase: is_diagonal_noise
import FastBroadcast: @..

import MuladdMacro: @muladd

import SciMLBase

using ArrayInterface: ArrayInterface
using LinearAlgebra: LinearAlgebra, UniformScaling, diagind, mul!
using StaticArrays: StaticArrays
using RecursiveArrayTools: RecursiveArrayTools
using StochasticDiffEqLevyArea: StochasticDiffEqLevyArea, MronRoe, levyarea,
    LevyAreaCoefficients
using DiffEqNoiseProcess: DiffEqNoiseProcess

include("algorithms.jl")
include("alg_utils.jl")

include("caches/basic_method_caches.jl")
include("caches/explicit_3s_mil_methods.jl")

include("perform_step/rkmilgeneral.jl")
include("perform_step/explicit_3s_mil_methods.jl")

export RKMilGeneral,
    WangLi3SMil_A, WangLi3SMil_B, WangLi3SMil_C,
    WangLi3SMil_D, WangLi3SMil_E, WangLi3SMil_F

end # module
