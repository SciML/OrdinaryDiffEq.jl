module OrdinaryDiffEqLinear

import OrdinaryDiffEqCore: alg_order, alg_extrapolates, dt_required,
                           OrdinaryDiffEqLinearExponentialAlgorithm,
                           OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqAlgorithm,
                           OrdinaryDiffEqExponentialAlgorithm,
                           OrdinaryDiffEqMutableCache, @cache, alg_cache,
                           OrdinaryDiffEqConstantCache,
                           initialize!, perform_step!, unwrap_alg,
                           calculate_residuals!, get_fsalfirstlast,
                           _vec, isdtchangeable, full_cache,
                           generic_solver_docstring, @SciMLMessage
using LinearAlgebra: mul!, I
using SciMLOperators: AbstractSciMLOperator
using ExponentialUtilities
using RecursiveArrayTools
import OrdinaryDiffEqCore
import DiffEqBase

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("linear_caches.jl")
include("integrator_interface.jl")
include("linear_perform_step.jl")

export MagnusMidpoint, LinearExponential, MagnusLeapfrog, LieEuler, CayleyEuler,
       MagnusGauss4, MagnusNC6, MagnusGL6, MagnusGL8, MagnusNC8, MagnusGL4,
       MagnusAdapt4, RKMK2, RKMK4, LieRK4, CG2, CG3, CG4a

end
