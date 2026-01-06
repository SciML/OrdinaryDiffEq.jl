module OrdinaryDiffEqRKIP

using Compat: logrange
using LinearAlgebra: ldiv!, exp, axpy!, norm, mul!
using SciMLOperators: AbstractSciMLOperator
using MaybeInplace: @bb
using SciMLBase: isinplace
using DiffEqBase: ExplicitRKTableau
using DiffEqDevTools: constructVerner6

import OrdinaryDiffEqCore: OrdinaryDiffEqAdaptiveExponentialAlgorithm, alg_adaptive_order,
    alg_order, alg_cache, @cache, SplitFunction, get_fsalfirstlast,
    initialize!, perform_step!,
    has_dtnew_modification, calculate_residuals,
    calculate_residuals!, increment_nf!,
    OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqMutableCache,
    dtnew_modification, generic_solver_docstring

include("rkip_cache.jl")
include("algorithms.jl")
include("rkip_utils.jl")
include("rkip_perform_step.jl")

export RKIP

end
