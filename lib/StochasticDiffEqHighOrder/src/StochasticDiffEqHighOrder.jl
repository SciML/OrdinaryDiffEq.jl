module StochasticDiffEqHighOrder

using Reexport: Reexport, @reexport
import StochasticDiffEqCore
@reexport using StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!

import DiffEqBase

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, alg_stability_size, delta_default,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    StochasticCompositeAlgorithm,
    get_iterated_I, get_iterated_I!,
    IICommutative, IILevyArea,
    calc_twopoint_random!, calc_twopoint_random,
    @cache

import DiffEqBase: @..
import DiffEqBase: calculate_residuals, calculate_residuals!
import DiffEqBase: full_cache, rand_cache, ratenoise_cache
import DiffEqBase: Tableau

import MuladdMacro: @muladd

import SciMLBase
import SciMLBase: is_diagonal_noise

using LinearAlgebra: LinearAlgebra, dot, mul!
using StaticArrays: StaticArrays, SArray
import RecursiveArrayTools

include("tableaus.jl")
include("algorithms.jl")
include("alg_utils.jl")

include("caches/rossler_caches.jl")
include("caches/sra_caches.jl")

include("perform_step/sri.jl")
include("perform_step/sra.jl")

export SRI, SRIW1, SRIW2, SOSRI, SOSRI2, SRA, SRA1, SRA2, SRA3, SOSRA, SOSRA2

# Tableau types, constructors and order checkers
export RosslerSRI, RosslerSRA,
    constructSRIW1, constructSRIW2, constructSRIOpt1, constructSRIOpt2,
    constructSRA1, constructSRA2, constructSRA3, constructSOSRA, constructSOSRA2,
    constructSKenCarp, constructExplicitSKenCarp,
    checkSRIOrder, checkSRAOrder

end # module
