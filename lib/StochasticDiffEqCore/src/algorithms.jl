"""
    IteratedIntegralApprox

Supertype for the strategies used to approximate the iterated stochastic integrals
(Lévy areas) that higher order SDE solvers need for non-diagonal noise.

An algorithm carries its choice in the `ii_approx` field, and
[`get_Jalg`](@ref) turns that choice plus the problem's noise structure into the
concrete `AbstractJ` object that computes the integrals.

Subtypes: [`IICommutative`](@ref), [`IILevyArea`](@ref).
"""
abstract type IteratedIntegralApprox end

"""
    IICommutative()

Iterated-integral strategy that assumes the noise is commutative.

Under commutative noise the Lévy area terms cancel and the iterated integrals reduce
to `1/2 ΔW ΔWᵀ`, which is cheap and exact for that problem class. Selecting this for
a genuinely non-commutative problem silently loses the method's strong order, so use
[`IILevyArea`](@ref) unless commutativity is known to hold.
"""
struct IICommutative <: IteratedIntegralApprox end

"""
    IILevyArea()

Iterated-integral strategy that approximates the Lévy area numerically.

This is the general-purpose choice: the number of terms is selected automatically
from the step size and the requested accuracy, and an appropriate simulation
algorithm is picked by `StochasticDiffEqLevyArea.optimal_algorithm`. It is more
expensive than [`IICommutative`](@ref) but does not assume anything about the noise
structure.
"""
struct IILevyArea <: IteratedIntegralApprox end

"""
    StochasticCompositeAlgorithm(algs, choice_function)

Algorithm that switches between the members of `algs` from step to step.

`choice_function(integrator)` returns the index into `algs` of the member to use for
the next step; the matching cache is held in a
[`StochasticCompositeCache`](@ref) and selected through its `current` field.

This is the mechanism behind the automatic stiffness-switching solvers — see
[`AutoAlgSwitch`](@ref), which pairs a nonstiff and a stiff algorithm with an
[`AutoSwitch`](@ref) choice function.
"""
struct StochasticCompositeAlgorithm{T, F} <: StochasticDiffEqCompositeAlgorithm
    algs::T
    choice_function::F
end

const SplitSDEAlgorithms = Union{}  # Extended by solver subpackages via Union
