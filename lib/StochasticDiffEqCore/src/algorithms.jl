# Abstract types are defined in OrdinaryDiffEqCore and imported via StochasticDiffEqCore.jl:
# StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
# StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODEAlgorithm,
# StochasticDiffEqRODEAdaptiveAlgorithm, StochasticDiffEqRODECompositeAlgorithm

abstract type StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller} <:
StochasticDiffEqAdaptiveAlgorithm end
abstract type StochasticDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ, Controller} <:
StochasticDiffEqAlgorithm end

abstract type StochasticDiffEqJumpAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpNewtonAdaptiveAlgorithm{
    CS, AD, FDT, ST, CJ, Controller,
} <: StochasticDiffEqJumpAdaptiveAlgorithm end

abstract type StochasticDiffEqJumpDiffusionAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpDiffusionAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm{
    CS, AD, FDT, ST, CJ, Controller,
} <: StochasticDiffEqJumpDiffusionAdaptiveAlgorithm end

abstract type IteratedIntegralApprox end

struct IICommutative <: IteratedIntegralApprox end
struct IILevyArea <: IteratedIntegralApprox end

struct StochasticCompositeAlgorithm{T, F} <: StochasticDiffEqCompositeAlgorithm
    algs::T
    choice_function::F
end

const SplitSDEAlgorithms = Union{}  # Extended by solver subpackages via Union
