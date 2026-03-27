# Abstract types are defined in OrdinaryDiffEqCore and imported via StochasticDiffEqCore.jl:
# StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
# StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODEAlgorithm,
# StochasticDiffEqRODEAdaptiveAlgorithm, StochasticDiffEqRODECompositeAlgorithm

abstract type StochasticDiffEqNewtonAdaptiveAlgorithm <:
StochasticDiffEqAdaptiveAlgorithm end
abstract type StochasticDiffEqNewtonAlgorithm <:
StochasticDiffEqAlgorithm end

abstract type StochasticDiffEqJumpAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpNewtonAdaptiveAlgorithm <: StochasticDiffEqJumpAdaptiveAlgorithm end

abstract type StochasticDiffEqJumpDiffusionAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpDiffusionAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm <: StochasticDiffEqJumpDiffusionAdaptiveAlgorithm end

abstract type IteratedIntegralApprox end

struct IICommutative <: IteratedIntegralApprox end
struct IILevyArea <: IteratedIntegralApprox end

struct StochasticCompositeAlgorithm{T, F} <: StochasticDiffEqCompositeAlgorithm
    algs::T
    choice_function::F
end

const SplitSDEAlgorithms = Union{}  # Extended by solver subpackages via Union
