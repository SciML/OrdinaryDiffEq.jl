module OrdinaryDiffEqLinear

import OrdinaryDiffEqCore: TmpCache, build_tmp_cache,
    alg_extrapolates, dt_required,
    OrdinaryDiffEqLinearExponentialAlgorithm,
    OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqAlgorithm,
    OrdinaryDiffEqMutableCache, @cache, alg_cache,
    OrdinaryDiffEqConstantCache,
    perform_step!, unwrap_alg,
    get_fsalfirstlast,
    isdtchangeable, full_cache,
    generic_solver_docstring
using LinearAlgebra: mul!, I
using SciMLOperators: AbstractSciMLOperator, update_coefficients, update_coefficients!
using ExponentialUtilities: ExponentialUtilities, ExpMethodGeneric, ExpvCache,
    KrylovSubspace, PhivCache, arnoldi!, exponential!,
    expv, expv!, expv_timestep, expv_timestep!
using RecursiveArrayTools: recursivefill!
import OrdinaryDiffEqCore
import DiffEqBase
import DiffEqBase: calculate_residuals!, initialize!
import SciMLBase: alg_order, _vec

using Reexport: Reexport, @reexport
@reexport using SciMLBase
using SciMLBase: SciMLBase, SciMLOperators, SplitFunction

include("algorithms.jl")
include("alg_utils.jl")
include("linear_caches.jl")
include("integrator_interface.jl")
include("linear_perform_step.jl")

export MagnusMidpoint, LinearExponential, MagnusLeapfrog, LieEuler, CayleyEuler,
    MagnusGauss4, MagnusNC6, MagnusGL6, MagnusGL8, MagnusNC8, MagnusGL4,
    MagnusAdapt4, RKMK2, RKMK4, LieRK4, CG2, CG3, CG4a

end
