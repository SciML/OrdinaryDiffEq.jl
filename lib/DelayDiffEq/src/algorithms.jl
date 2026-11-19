abstract type DelayDiffEqAlgorithm <: AbstractDDEAlgorithm end
abstract type AbstractMethodOfStepsAlgorithm{constrained} <: DelayDiffEqAlgorithm end

struct MethodOfSteps{algType, F, constrained} <: AbstractMethodOfStepsAlgorithm{constrained}
    alg::algType
    fpsolve::F
end

"""
    MethodOfSteps(alg; constrained = false, fpsolve = NLFunctional())

Construct an algorithm that solves delay differential equations by the method of steps,
where `alg` is an ODE algorithm from OrdinaryDiffEq.jl upon which the calculation of
steps is based.

If the algorithm is `constrained`, steps are at most the smallest nonzero constant lag, both for adaptive and fixed time steps.
Fixed-point iteration `fpsolve` is applied to every step that needs the solution after the start of the step, for example a step longer than a dependent lag.
`fpsolve` can be `NLFunctional()` or `NLAnderson()`.


Citations:

General Approach

Zivari-Piran, Hossein, and Wayne H. Enright. "An efficient unified approach for the numerical 
solution of delay differential equations." Numerical Algorithms 53.2-3 (2010): 397-417.

State-Dependent Delays

S. P. Corwin, D. Sarafyan and S. Thompson in "DKLAG6: a code based on continuously embedded
sixth-order Runge-Kutta methods for the solution of state-dependent functional differential
equations", Applied Numerical Mathematics, 1997.
"""
function MethodOfSteps(alg; constrained = false, fpsolve = NLFunctional())
    return MethodOfSteps{typeof(alg), typeof(fpsolve), constrained}(alg, fpsolve)
end
