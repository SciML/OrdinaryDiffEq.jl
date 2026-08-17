"""
$(DocStringExtensions.README)
"""
module OrdinaryDiffEq

import DocStringExtensions

# Load packages (no blanket @reexport)

# Import specific algorithms from OrdinaryDiffEqDefault's dependencies
using OrdinaryDiffEqTsit5: Tsit5, AutoTsit5
using OrdinaryDiffEqVerner: Vern6, Vern7, Vern8, Vern9,
    AutoVern6, AutoVern7, AutoVern8, AutoVern9
using OrdinaryDiffEqRosenbrock: Rosenbrock23, Rodas5P
using OrdinaryDiffEqBDF: FBDF

# Import ODE-relevant types from SciMLBase
using SciMLBase: SciMLBase,
    ODEProblem, ODEFunction, ODESolution,
    SplitODEProblem, SplitFunction,
    SecondOrderODEProblem, DynamicalODEProblem,
    DAEProblem, DAEFunction, DAESolution, EnsembleProblem,
    CallbackSet, ContinuousCallback, DiscreteCallback, VectorContinuousCallback,
    ReturnCode, set_proposed_dt!,
    remake, successful_retcode, reinit!, terminate!,
    derivative_discontinuity!, add_tstop!, ODEAliasSpecifier

using SciMLLogging: SciMLLogging

"""
    AutoDespecialize

Use a stable dynamic parameter container so compiled solver code can be shared across
different concrete parameter layouts. The original parameter is recovered at the SciML
function call barrier.
"""
const AutoDespecialize = SciMLBase.AutoDespecialize

"""
    AutoRespecialize

Use the constrained non-dynamic parameter policy, which recovers the original concrete
parameter type from an opaque container on supported solver paths.
"""
const AutoRespecialize = SciMLBase.AutoRespecialize

"""
    AutoDePSpecialize

Deprecated alias for [`AutoRespecialize`](@ref).
"""
const AutoDePSpecialize = SciMLBase.AutoDePSpecialize

# Verbosity specifier owned by DiffEqBase (exported there). The umbrella package
# imports and exports it directly so `verbose = DEVerbosity(...)` is reachable from
# a plain `using OrdinaryDiffEq`.
using DiffEqBase: DEVerbosity

# Import ADTypes for autodiff specification
using ADTypes: ADTypes, AutoForwardDiff, AutoFiniteDiff, AutoSparse

# Import from OrdinaryDiffEqCore
using OrdinaryDiffEqCore: OrdinaryDiffEqCore

"""
    Sequential() <: OrdinaryDiffEqCore.AbstractThreadingOption

Use one thread for solver work that can otherwise be executed in parallel.
This is the default threading option for deterministic execution.
"""
const Sequential = OrdinaryDiffEqCore.Sequential

"""
    BaseThreads() <: OrdinaryDiffEqCore.AbstractThreadingOption

Use Julia's built-in `Threads.@threads` for solver work that can be executed in
parallel. The active Julia process must have more than one thread for this to
provide parallelism.
"""
const BaseThreads = OrdinaryDiffEqCore.BaseThreads

"""
    PolyesterThreads() <: OrdinaryDiffEqCore.AbstractThreadingOption

Use Polyester.jl's low-overhead threaded execution for solver work that can be
executed in parallel. Load `Polyester` before selecting this option.
"""
const PolyesterThreads = OrdinaryDiffEqCore.PolyesterThreads

# Import from OrdinaryDiffEqDefault
using OrdinaryDiffEqDefault: DefaultODEAlgorithm

import CommonSolve: init, solve, solve!, step!

# --- Exports ---

# General Functions
export solve, solve!, init, step!

# Problem types
export ODEProblem, ODEFunction, ODESolution
export SplitODEProblem, SplitFunction
export SecondOrderODEProblem, DynamicalODEProblem
export DAEProblem, DAEFunction, EnsembleProblem, DAESolution

# Callbacks
export CallbackSet, ContinuousCallback, DiscreteCallback, VectorContinuousCallback

# Utilities
export ReturnCode, derivative_discontinuity!, add_tstop!, ODEAliasSpecifier
export SciMLBase, SciMLLogging, remake, successful_retcode, reinit!, set_proposed_dt!, terminate!

# Verbosity
export DEVerbosity

# Specialization levels
export AutoDespecialize, AutoRespecialize, AutoDePSpecialize

# ADTypes
export AutoForwardDiff, AutoFiniteDiff, AutoSparse

# Default algorithm
export DefaultODEAlgorithm

# Widely-used algorithms
export Tsit5, AutoTsit5
export Vern6, Vern7, Vern8, Vern9, AutoVern6, AutoVern7, AutoVern8, AutoVern9
export Rosenbrock23, Rodas5P
export FBDF

# Reachable as `OrdinaryDiffEq.PolyesterThreads()`, the spelling v6 supported.
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Expr(:public, :Sequential, :BaseThreads, :PolyesterThreads))
end

end # module
