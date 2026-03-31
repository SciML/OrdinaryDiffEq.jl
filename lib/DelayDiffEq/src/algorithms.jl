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

If the algorithm is `constrained` only steps of size at most the minimal delay will be
taken. If it is unconstrained, fixed-point iteration `fpsolve` is applied for step sizes
that exceed the minimal delay.


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
