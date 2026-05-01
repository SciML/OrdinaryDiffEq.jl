abstract type IteratedIntegralApprox end

struct IICommutative <: IteratedIntegralApprox end
struct IILevyArea <: IteratedIntegralApprox end

struct StochasticCompositeAlgorithm{T, F} <: StochasticDiffEqCompositeAlgorithm
    algs::T
    choice_function::F
end

const SplitSDEAlgorithms = Union{}  # Extended by solver subpackages via Union
