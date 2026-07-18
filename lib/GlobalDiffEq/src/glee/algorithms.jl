# Explicit general linear methods with built-in global error estimation
# (GLEE methods) from Constantinescu (2016), doi:10.1137/15M1014633.
abstract type AbstractGLEEAlgorithm <:
OrdinaryDiffEqCore.OrdinaryDiffEqAdaptiveAlgorithm end

const _GLEE_DOCS_SHARED = """
GLEE methods propagate the solution `y` together with an asymptotically
correct estimate `ε` of its global (accumulated) error, at the cost of a few
extra stages per step. Solving any `ODEProblem` with a GLEE method produces a
solution whose states are `ArrayPartition`s: `sol.u[i].x[1]` is the solution
and `sol.u[i].x[2]` is the global error estimate at `sol.t[i]` (see
[`global_error_estimate`](@ref)). The per-step increment of `ε` is an
asymptotically correct local error estimate, which drives standard step-size
adaptivity, so local tolerances behave exactly as for ordinary adaptive
Runge-Kutta methods while the global error is estimated for free.

Only explicit, mass-matrix-free ODEs are supported. The reference for the
methods and their theory is:

Emil M. Constantinescu, *Estimating Global Errors in Time Stepping*, SIAM
Journal on Numerical Analysis 54(6), 2016. [arXiv:1503.05166](https://arxiv.org/abs/1503.05166)
"""

"""
    GLEE23()

3-stage, second-order explicit general linear method with global error
estimation (Constantinescu 2016, eq. (4.6); PETSc's `TSGLEE23`). The cheapest
GLEE method: only the first decoupling condition holds, so prefer
[`GLEE24`](@ref) for long-time integration or mildly stiff problems.

$(_GLEE_DOCS_SHARED)
"""
struct GLEE23 <: AbstractGLEEAlgorithm end

"""
    GLEE24()

4-stage, second-order explicit general linear method with global error
estimation (Constantinescu 2016, eq. (A.3); PETSc's `TSGLEE24`). Satisfies
both decoupling conditions (`B·U` and `B·A·U` diagonal), which keeps the error
estimate faithful in long-time integration; this is the recommended
second-order GLEE method.

$(_GLEE_DOCS_SHARED)
"""
struct GLEE24 <: AbstractGLEEAlgorithm end

"""
    GLEE35()

5-stage, third-order explicit general linear method with global error
estimation (Constantinescu 2016, eq. (4.9); PETSc's `TSGLEE35`). Satisfies
both decoupling conditions and has a large negative-real-axis stability
region; this is the recommended third-order GLEE method.

$(_GLEE_DOCS_SHARED)
"""
struct GLEE35 <: AbstractGLEEAlgorithm end

SciMLBase.alg_order(::GLEE23) = 2
SciMLBase.alg_order(::GLEE24) = 2
SciMLBase.alg_order(::GLEE35) = 3

OrdinaryDiffEqCore.alg_adaptive_order(alg::AbstractGLEEAlgorithm) =
    SciMLBase.alg_order(alg)

_glee_tableau_for(::GLEE23, ::Type{T}, ::Type{T2}) where {T, T2} = GLEE23Tableau(T, T2)
_glee_tableau_for(::GLEE24, ::Type{T}, ::Type{T2}) where {T, T2} = GLEE24Tableau(T, T2)
_glee_tableau_for(::GLEE35, ::Type{T}, ::Type{T2}) where {T, T2} = GLEE35Tableau(T, T2)
