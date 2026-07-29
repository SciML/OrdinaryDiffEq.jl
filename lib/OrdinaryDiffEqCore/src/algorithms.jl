"""
    OrdinaryDiffEqAlgorithm <: SciMLBase.AbstractODEAlgorithm

Abstract supertype of every ODE algorithm defined in the OrdinaryDiffEq.jl
ecosystem. A solver sublibrary defines a concrete algorithm by subtyping one of
the more specific abstract types below (adaptive / implicit / exponential / …)
rather than this root directly. Trait functions such as [`isadaptive`](@ref),
[`isimplicit`](@ref), [`isfsal`](@ref), and [`alg_order`](@ref) dispatch on this
hierarchy.
"""
abstract type OrdinaryDiffEqAlgorithm <: SciMLBase.AbstractODEAlgorithm end
"""
    OrdinaryDiffEqAdaptiveAlgorithm <: OrdinaryDiffEqAlgorithm

Abstract supertype for explicit ODE algorithms that support adaptive step-size
control (they carry an embedded error estimate). Subtyping this makes
[`isadaptive`](@ref) return `true`.
"""
abstract type OrdinaryDiffEqAdaptiveAlgorithm <: OrdinaryDiffEqAlgorithm end
"""
    OrdinaryDiffEqCompositeAlgorithm <: OrdinaryDiffEqAlgorithm

Abstract supertype for algorithms that dispatch between several sub-algorithms at
runtime (see [`CompositeAlgorithm`](@ref)). Used by automatic stiffness switching
and by the default solver.
"""
abstract type OrdinaryDiffEqCompositeAlgorithm <: OrdinaryDiffEqAlgorithm end

# SDE/RODE algorithm type hierarchy (used by StochasticDiffEq)
"""
    StochasticDiffEqAlgorithm <: SciMLBase.AbstractSDEAlgorithm

Abstract supertype of every SDE algorithm in the StochasticDiffEq.jl ecosystem.
Defined here in OrdinaryDiffEqCore so that the shared integrator machinery can
dispatch on it.
"""
abstract type StochasticDiffEqAlgorithm <: SciMLBase.AbstractSDEAlgorithm end
"""
    StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm

Abstract supertype for SDE algorithms supporting adaptive step-size control.
"""
abstract type StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
"""
    StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm

Abstract supertype for composite (runtime-switching) SDE algorithms.
"""
abstract type StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm end

"""
    StochasticDiffEqRODEAlgorithm <: SciMLBase.AbstractRODEAlgorithm

Abstract supertype for random ODE (RODE) algorithms.
"""
abstract type StochasticDiffEqRODEAlgorithm <: SciMLBase.AbstractRODEAlgorithm end
"""
    StochasticDiffEqRODEAdaptiveAlgorithm <: StochasticDiffEqRODEAlgorithm

Abstract supertype for adaptive RODE algorithms.
"""
abstract type StochasticDiffEqRODEAdaptiveAlgorithm <: StochasticDiffEqRODEAlgorithm end
"""
    StochasticDiffEqRODECompositeAlgorithm <: StochasticDiffEqRODEAlgorithm

Abstract supertype for composite (runtime-switching) RODE algorithms.
"""
abstract type StochasticDiffEqRODECompositeAlgorithm <: StochasticDiffEqRODEAlgorithm end

# SDE Newton/Jump algorithm subtypes (used by StochasticDiffEq implicit solvers)
"""
    StochasticDiffEqNewtonAdaptiveAlgorithm <: StochasticDiffEqAdaptiveAlgorithm

Abstract supertype for adaptive implicit SDE algorithms solved with a Newton
nonlinear solver.
"""
abstract type StochasticDiffEqNewtonAdaptiveAlgorithm <:
StochasticDiffEqAdaptiveAlgorithm end
"""
    StochasticDiffEqNewtonAlgorithm <: StochasticDiffEqAlgorithm

Abstract supertype for fixed-step implicit SDE algorithms solved with a Newton
nonlinear solver.
"""
abstract type StochasticDiffEqNewtonAlgorithm <:
StochasticDiffEqAlgorithm end

"""
    StochasticDiffEqJumpAlgorithm <: StochasticDiffEqAlgorithm

Abstract supertype for SDE algorithms that additionally integrate jump terms.
"""
abstract type StochasticDiffEqJumpAlgorithm <: StochasticDiffEqAlgorithm end
"""
    StochasticDiffEqJumpAdaptiveAlgorithm <: StochasticDiffEqAlgorithm

Abstract supertype for adaptive jump-SDE algorithms.
"""
abstract type StochasticDiffEqJumpAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
"""
    StochasticDiffEqJumpNewtonAdaptiveAlgorithm <: StochasticDiffEqJumpAdaptiveAlgorithm

Abstract supertype for adaptive implicit (Newton) jump-SDE algorithms.
"""
abstract type StochasticDiffEqJumpNewtonAdaptiveAlgorithm <: StochasticDiffEqJumpAdaptiveAlgorithm end

"""
    StochasticDiffEqJumpDiffusionAlgorithm <: StochasticDiffEqAlgorithm

Abstract supertype for jump-diffusion algorithms.
"""
abstract type StochasticDiffEqJumpDiffusionAlgorithm <: StochasticDiffEqAlgorithm end
"""
    StochasticDiffEqJumpDiffusionAdaptiveAlgorithm <: StochasticDiffEqAlgorithm

Abstract supertype for adaptive jump-diffusion algorithms.
"""
abstract type StochasticDiffEqJumpDiffusionAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
"""
    StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm <: StochasticDiffEqJumpDiffusionAdaptiveAlgorithm

Abstract supertype for adaptive implicit (Newton) jump-diffusion algorithms.
"""
abstract type StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm <: StochasticDiffEqJumpDiffusionAdaptiveAlgorithm end

# SDE/RODE cache type hierarchy (used by StochasticDiffEq)
"""
    StochasticDiffEqCache <: SciMLBase.DECache

Abstract supertype of every SDE/RODE solver cache. Analogue of
[`OrdinaryDiffEqCache`](@ref) for the StochasticDiffEq.jl solvers; concrete caches
subtype [`StochasticDiffEqConstantCache`](@ref) or
[`StochasticDiffEqMutableCache`](@ref).
"""
abstract type StochasticDiffEqCache <: SciMLBase.DECache end
"""
    StochasticDiffEqConstantCache <: StochasticDiffEqCache

Abstract supertype for out-of-place SDE/RODE caches (immutable state).
"""
abstract type StochasticDiffEqConstantCache <: StochasticDiffEqCache end
"""
    StochasticDiffEqMutableCache <: StochasticDiffEqCache

Abstract supertype for in-place SDE/RODE caches (preallocated, mutated in place).
"""
abstract type StochasticDiffEqMutableCache <: StochasticDiffEqCache end

"""
    OrdinaryDiffEqAdaptiveImplicitAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm

Abstract supertype for implicit ODE algorithms with adaptive step-size control.
"""
abstract type OrdinaryDiffEqAdaptiveImplicitAlgorithm <:
OrdinaryDiffEqAdaptiveAlgorithm end
"""
    OrdinaryDiffEqNewtonAdaptiveAlgorithm <: OrdinaryDiffEqAdaptiveImplicitAlgorithm

Abstract supertype for adaptive implicit algorithms solved with a Newton
nonlinear solver. Distinguished from [`OrdinaryDiffEqNewtonAlgorithm`](@ref) by
supporting Jacobian reuse across steps ([`alg_can_repeat_jac`](@ref) is `true`).
"""
abstract type OrdinaryDiffEqNewtonAdaptiveAlgorithm <:
OrdinaryDiffEqAdaptiveImplicitAlgorithm end
"""
    OrdinaryDiffEqRosenbrockAdaptiveAlgorithm <: OrdinaryDiffEqAdaptiveImplicitAlgorithm

Abstract supertype for adaptive Rosenbrock / Rosenbrock-W methods.
"""
abstract type OrdinaryDiffEqRosenbrockAdaptiveAlgorithm <:
OrdinaryDiffEqAdaptiveImplicitAlgorithm end

"""
    OrdinaryDiffEqImplicitAlgorithm <: OrdinaryDiffEqAlgorithm

Abstract supertype for implicit ODE algorithms (those that solve a nonlinear
system each step). Subtyping this makes [`isimplicit`](@ref) return `true`.
"""
abstract type OrdinaryDiffEqImplicitAlgorithm <:
OrdinaryDiffEqAlgorithm end
"""
    OrdinaryDiffEqNewtonAlgorithm <: OrdinaryDiffEqImplicitAlgorithm

Abstract supertype for fixed-step implicit algorithms whose implicit stages are
solved with a (quasi-)Newton nonlinear solver (see [`AbstractNLSolver`](@ref)).
"""
abstract type OrdinaryDiffEqNewtonAlgorithm <:
OrdinaryDiffEqImplicitAlgorithm end
"""
    OrdinaryDiffEqRosenbrockAlgorithm <: OrdinaryDiffEqImplicitAlgorithm

Abstract supertype for fixed-step Rosenbrock (and Rosenbrock-W) methods, which
use the Jacobian directly through linear solves rather than a Newton iteration.
"""
abstract type OrdinaryDiffEqRosenbrockAlgorithm <:
OrdinaryDiffEqImplicitAlgorithm end
"""
    NewtonAlgorithm

`Union` of [`OrdinaryDiffEqNewtonAlgorithm`](@ref) and
[`OrdinaryDiffEqNewtonAdaptiveAlgorithm`](@ref), i.e. every implicit algorithm
that drives its stages with a Newton nonlinear solver. Used as a dispatch handle
for the Newton/W-matrix machinery.
"""
const NewtonAlgorithm = Union{
    OrdinaryDiffEqNewtonAlgorithm,
    OrdinaryDiffEqNewtonAdaptiveAlgorithm,
}
"""
    RosenbrockAlgorithm

`Union` of [`OrdinaryDiffEqRosenbrockAlgorithm`](@ref) and
[`OrdinaryDiffEqRosenbrockAdaptiveAlgorithm`](@ref), i.e. every Rosenbrock /
Rosenbrock-W method.
"""
const RosenbrockAlgorithm = Union{
    OrdinaryDiffEqRosenbrockAlgorithm,
    OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
}

"""
    OrdinaryDiffEqExponentialAlgorithm <: OrdinaryDiffEqAlgorithm

Abstract supertype for exponential integrators that evaluate matrix-function
(e.g. `exp`, `φ`) actions of the (linearized) operator each step.
"""
abstract type OrdinaryDiffEqExponentialAlgorithm <:
OrdinaryDiffEqAlgorithm end
"""
    OrdinaryDiffEqAdaptiveExponentialAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm

Abstract supertype for exponential integrators with adaptive step-size control.
"""
abstract type OrdinaryDiffEqAdaptiveExponentialAlgorithm <:
OrdinaryDiffEqAdaptiveAlgorithm end
"""
    OrdinaryDiffEqLinearExponentialAlgorithm <: OrdinaryDiffEqExponentialAlgorithm

Abstract supertype for exponential integrators specialized to (semi)linear
problems `u' = A u (+ B(t))` where the operator action can be applied directly.
"""
abstract type OrdinaryDiffEqLinearExponentialAlgorithm <:
OrdinaryDiffEqExponentialAlgorithm end
"""
    ExponentialAlgorithm

`Union` of [`OrdinaryDiffEqExponentialAlgorithm`](@ref) and
[`OrdinaryDiffEqAdaptiveExponentialAlgorithm`](@ref).
"""
const ExponentialAlgorithm = Union{
    OrdinaryDiffEqExponentialAlgorithm,
    OrdinaryDiffEqAdaptiveExponentialAlgorithm,
}

"""
    OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm

Abstract supertype for variable-order variable-step Adams (multistep) methods.
For these algorithms the current order lives on the cache, so
[`get_current_alg_order`](@ref) / [`get_current_adaptive_order`](@ref) read
`cache.order` rather than a fixed algorithm order.
"""
abstract type OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end

# DAE Specific Algorithms
"""
    DAEAlgorithm <: SciMLBase.AbstractDAEAlgorithm

Abstract supertype for fully-implicit DAE algorithms (solving
`f(du, u, p, t) = 0`), e.g. `DFBDF`, `DImplicitEuler`. Distinct from the mass-matrix
DAE path, which reuses the ODE algorithm types with a singular mass matrix.
"""
abstract type DAEAlgorithm <: SciMLBase.AbstractDAEAlgorithm end

# Partitioned ODE Specific Algorithms
"""
    OrdinaryDiffEqPartitionedAlgorithm <: OrdinaryDiffEqAlgorithm

Abstract supertype for partitioned / dynamical-ODE algorithms (e.g. symplectic
and Nyström methods) that split the state into position/velocity partitions.
"""
abstract type OrdinaryDiffEqPartitionedAlgorithm <: OrdinaryDiffEqAlgorithm end
"""
    OrdinaryDiffEqAdaptivePartitionedAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm

Abstract supertype for partitioned / dynamical-ODE algorithms with adaptive
step-size control.
"""
abstract type OrdinaryDiffEqAdaptivePartitionedAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
"""
    PartitionedAlgorithm

`Union` of [`OrdinaryDiffEqPartitionedAlgorithm`](@ref) and
[`OrdinaryDiffEqAdaptivePartitionedAlgorithm`](@ref).
"""
const PartitionedAlgorithm = Union{
    OrdinaryDiffEqPartitionedAlgorithm,
    OrdinaryDiffEqAdaptivePartitionedAlgorithm,
}

# Second order ODE Specific Algorithms
"""
    OrdinaryDiffEqImplicitSecondOrderAlgorithm <: OrdinaryDiffEqImplicitAlgorithm

Abstract supertype for implicit second-order (dynamical) ODE algorithms.
"""
abstract type OrdinaryDiffEqImplicitSecondOrderAlgorithm <:
OrdinaryDiffEqImplicitAlgorithm end
"""
    OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm <: OrdinaryDiffEqAdaptiveImplicitAlgorithm

Abstract supertype for adaptive implicit second-order (dynamical) ODE algorithms.
"""
abstract type OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm <:
OrdinaryDiffEqAdaptiveImplicitAlgorithm end
"""
    ImplicitSecondOrderAlgorithm

`Union` of [`OrdinaryDiffEqImplicitSecondOrderAlgorithm`](@ref) and
[`OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm`](@ref).
"""
const ImplicitSecondOrderAlgorithm = Union{
    OrdinaryDiffEqImplicitSecondOrderAlgorithm,
    OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm,
}

function SciMLBase.remake(thing::OrdinaryDiffEqAlgorithm; kwargs...)
    T = SciMLBase.remaker_of(thing)
    return T(; SciMLBase.struct_as_namedtuple(thing)..., kwargs...)
end


###############################################################################

################################################################################

@inline trivial_limiter!(u, integrator, p, t) = nothing

################################################################################

################################################################################

################################################################################

######################################

#########################################

#########################################

"""
    CompositeAlgorithm(algs, choice_function)

A composite algorithm that chooses between multiple ODE solvers based on a user-defined choice function.
This allows for adaptive algorithm switching based on problem characteristics or performance metrics.

# Arguments

  - `algs`: Tuple or array of ODE algorithms to choose from
  - `choice_function`: Function that determines which algorithm to use at each step

The choice function receives the integrator and should return an index indicating which algorithm to use.
This enables sophisticated algorithm switching strategies based on solution behavior, step size, or other criteria.
"""
struct CompositeAlgorithm{CS, T, F} <: OrdinaryDiffEqCompositeAlgorithm
    algs::T
    choice_function::F
    function CompositeAlgorithm(algs::T, choice_function::F) where {T, F}
        CS = mapreduce(alg -> 0, max, algs)
        return new{CS, T, F}(algs, choice_function)
    end
end

@truncate_stacktrace CompositeAlgorithm 1

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, :silence!)
    Base.Experimental.silence!(CompositeAlgorithm)
end

"""
    AutoSwitchCache

Mutable state backing a stiffness-based [`AutoSwitch`](@ref) choice function. It
counts consecutive stiff/nonstiff step diagnostics (`count`,
`successive_switches`), stores the two candidate algorithms
(`nonstiffalg`/`stiffalg`), the switching thresholds
(`maxstiffstep`/`maxnonstiffstep`, `nonstifftol`/`stifftol`), and which branch is
currently active (`is_stiffalg`, `current`). Called on the integrator to return
the index of the algorithm to use for the next step.
"""
mutable struct AutoSwitchCache{nAlg, sAlg, tolType, T}
    count::Int
    successive_switches::Int
    nonstiffalg::nAlg
    stiffalg::sAlg
    is_stiffalg::Bool
    maxstiffstep::Int
    maxnonstiffstep::Int
    nonstifftol::tolType
    stifftol::tolType
    dtfac::T
    stiffalgfirst::Bool
    switch_max::Int
    current::Int
    function AutoSwitchCache(
            count::Int,
            successive_switches::Int,
            nonstiffalg::nAlg,
            stiffalg::sAlg,
            is_stiffalg::Bool,
            maxstiffstep::Int,
            maxnonstiffstep::Int,
            nonstifftol::tolType,
            stifftol::tolType,
            dtfac::T,
            stiffalgfirst::Bool,
            switch_max::Int,
            current::Int = 0
        ) where {nAlg, sAlg, tolType, T}
        return new{nAlg, sAlg, tolType, T}(
            count,
            successive_switches,
            nonstiffalg,
            stiffalg,
            is_stiffalg,
            maxstiffstep,
            maxnonstiffstep,
            nonstifftol,
            stifftol,
            dtfac,
            stiffalgfirst,
            switch_max,
            current
        )
    end
end

"""
    AutoSwitch(nonstiffalg, stiffalg; kwargs...)

An automatic algorithm switching method that dynamically chooses between a nonstiff and stiff solver
based on the problem's stiffness detection. This provides robust performance across a wide range of problems
without requiring the user to know the problem's stiffness characteristics a priori.

# Arguments

  - `nonstiffalg`: Algorithm to use for nonstiff regions (default: Tsit5())
  - `stiffalg`: Algorithm to use for stiff regions (default: Rodas5P())

# Keywords

  - `maxstiffstep`: Maximum number of consecutive steps before switching from nonstiff to stiff (default: 10)
  - `maxnonstiffstep`: Maximum number of consecutive steps before switching from stiff to nonstiff (default: 3)
  - `nonstifftol`: Tolerance for detecting nonstiff behavior (default: 3//4)
  - `stifftol`: Tolerance for detecting stiff behavior (default: 9//10)
  - `dtfac`: Factor for step size adjustment during switches (default: 2.0)
  - `stiffalgfirst`: Whether to start with the stiff algorithm (default: false)
  - `switch_max`: Maximum number of algorithm switches allowed (default: 10)

The switching decision is based on step size rejections and stability estimates.
"""
struct AutoSwitch{nAlg, sAlg, tolType, T}
    nonstiffalg::nAlg
    stiffalg::sAlg
    maxstiffstep::Int
    maxnonstiffstep::Int
    nonstifftol::tolType
    stifftol::tolType
    dtfac::T
    stiffalgfirst::Bool
    switch_max::Int
end
