module OrdinaryDiffEqDifferentiation

import ADTypes: AutoFiniteDiff, AutoForwardDiff

import SparseDiffTools: SparseDiffTools, matrix_colors, forwarddiff_color_jacobian!,
                forwarddiff_color_jacobian, ForwardColorJacCache,
                default_chunk_size, getsize, JacVec

import ForwardDiff, FiniteDiff
import ForwardDiff.Dual
import LinearSolve
import LinearSolve: OperatorAssumptions
using DiffEqBase
import LinearAlgebra
import InteractiveUtils

using DiffEqBase: TimeGradientWrapper,
                  UJacobianWrapper, TimeDerivativeWrapper,
                  UDerivativeWrapper
using SciMLBase: AbstractSciMLOperator
using OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveImplicitAlgorithm, DAEAlgorithm,
OrdinaryDiffEqImplicitAlgorithm, CompositeAlgorithm, OrdinaryDiffEqExponentialAlgorithm,
OrdinaryDiffEqAdaptiveExponentialAlgorithm, @unpack, AbstractNLSolver

import OrdinaryDiffEq: alg_autodiff

using FastBroadcast: @..

@static if isdefined(DiffEqBase, :OrdinaryDiffEqTag)
    import DiffEqBase: OrdinaryDiffEqTag
else
    struct OrdinaryDiffEqTag end
end

include("alg_utils.jl")
include("linsolve_utils.jl")
include("derivative_utils.jl")
include("derivative_wrappers.jl")

end
