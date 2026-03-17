module StochasticDiffEqHighOrder

using Reexport
@reexport using StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!, initialize!

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, alg_stability_size, delta_default,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    StochasticCompositeAlgorithm,
    get_iterated_I, get_iterated_I!,
    IICommutative, IILevyArea,
    calc_twopoint_random!, calc_twopoint_random,
    @cache

import DiffEqBase: is_diagonal_noise, @..
import DiffEqBase: calculate_residuals, calculate_residuals!
import DiffEqBase: full_cache, rand_cache, ratenoise_cache
import DiffEqBase: Tableau

import MuladdMacro: @muladd

import SciMLBase

using LinearAlgebra
using StaticArrays
using RecursiveArrayTools

include("tableaus.jl")
include("algorithms.jl")
include("alg_utils.jl")

include("caches/rossler_caches.jl")
include("caches/sra_caches.jl")

include("perform_step/sri.jl")
include("perform_step/sra.jl")

export SRI, SRIW1, SRIW2, SOSRI, SOSRI2, SRA, SRA1, SRA2, SRA3, SOSRA, SOSRA2

end # module
