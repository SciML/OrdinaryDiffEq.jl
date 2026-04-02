"""
$(DocStringExtensions.README)
"""
module OrdinaryDiffEq

import DocStringExtensions

# Load packages (no blanket @reexport)
using OrdinaryDiffEqCore
using OrdinaryDiffEqDefault

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
    DAEProblem, DAEFunction, DAESolution,
    CallbackSet, ContinuousCallback, DiscreteCallback, VectorContinuousCallback,
    ReturnCode, set_proposed_dt!,
    remake, successful_retcode, reinit!

# Import ADTypes for autodiff specification
using ADTypes: ADTypes, AutoForwardDiff, AutoFiniteDiff, AutoSparse

# Import from OrdinaryDiffEqCore
using OrdinaryDiffEqCore: OrdinaryDiffEqCore,
    CompositeAlgorithm, AutoSwitch

# Import from OrdinaryDiffEqDefault
using OrdinaryDiffEqDefault: DefaultODEAlgorithm, DefaultImplicitODEAlgorithm

import CommonSolve: init, solve, solve!, step!

# --- Exports ---

# General Functions
export solve, solve!, init, step!

# Problem types
export ODEProblem, ODEFunction, ODESolution
export SplitODEProblem, SplitFunction
export SecondOrderODEProblem, DynamicalODEProblem
export DAEProblem, DAEFunction, DAESolution

# Callbacks
export CallbackSet, ContinuousCallback, DiscreteCallback, VectorContinuousCallback

# Utilities
export ReturnCode, u_modified!, add_tstop!, ODEAliasSpecifier
export remake, successful_retcode, reinit!, set_proposed_dt!

# ADTypes
export AutoForwardDiff, AutoFiniteDiff, AutoSparse

# Default algorithm
export DefaultODEAlgorithm

# Widely-used algorithms
export Tsit5, AutoTsit5
export Vern6, Vern7, Vern8, Vern9, AutoVern6, AutoVern7, AutoVern8, AutoVern9
export Rosenbrock23, Rodas5P
export FBDF

end # module
