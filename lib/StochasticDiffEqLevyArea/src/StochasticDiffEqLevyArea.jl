module StochasticDiffEqLevyArea

import LinearAlgebra: mul!
import Random: randn!, default_rng, AbstractRNG
using Random: Xoshiro
import SpecialFunctions: trigamma

"""
    abstract type AbstractIteratedIntegralAlgorithm end

Abstract type for algorithms for the simulation of iterated stochastic integrals.
"""
abstract type AbstractIteratedIntegralAlgorithm end

# Error norms
include("error_norms.jl")
export MaxL2, FrobeniusL2, AbstractErrorNorm

# Coefficient struct (no algorithm-type dependencies)
include("coefficients.jl")
export LevyAreaCoefficients, coefficient_length

# Algorithm implementations (define levyarea, convorder, errcoeff, norv)
include("fourier.jl")
include("milstein.jl")
include("wiktorsson.jl")
include("mronroe.jl")
export AbstractIteratedIntegralAlgorithm
export Fourier, Milstein, Wiktorsson, MronRoe
export levyarea

const ITER_INT_ALGS = [Fourier(), Milstein(), Wiktorsson(), MronRoe()]

# Coefficient generation (needs algorithm types)
include("generate.jl")
export generate_coefficients

# Algorithm selection utilities (needs ITER_INT_ALGS)
include("alg_utils.jl")
export terms_needed, optimal_algorithm

# Public API for iterated integrals
include("iterated_integrals.jl")
export iterated_integrals

# Path reconstruction from Fourier coefficients
include("reconstruct.jl")
export reconstruct_path, iterated_integrals_subinterval

end # module
