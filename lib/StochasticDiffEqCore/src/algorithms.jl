# All abstract algorithm types are now defined in OrdinaryDiffEqCore and imported
# via StochasticDiffEqCore.jl (including Newton/Jump subtypes).

abstract type IteratedIntegralApprox end

struct IICommutative <: IteratedIntegralApprox end
struct IILevyArea <: IteratedIntegralApprox end

struct StochasticCompositeAlgorithm{T, F} <: StochasticDiffEqCompositeAlgorithm
    algs::T
    choice_function::F
end

const SplitSDEAlgorithms = Union{}  # Extended by solver subpackages via Union
