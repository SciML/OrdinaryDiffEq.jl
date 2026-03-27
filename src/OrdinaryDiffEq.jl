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

# Import ODE-relevant types from SciMLBase (NOT blanket reexport)
using SciMLBase: SciMLBase,
    ODEProblem, ODEFunction, ODESolution,
    SplitODEProblem, SplitFunction,
    SecondOrderODEProblem, DynamicalODEProblem,
    DAEProblem, DAEFunction, DAESolution,
    DiscreteProblem, DiscreteFunction,
    EnsembleProblem, EnsembleAlgorithm, EnsembleSolution,
    EnsembleThreads, EnsembleDistributed, EnsembleSerial, EnsembleSplitThreads,
    CallbackSet, ContinuousCallback, DiscreteCallback, VectorContinuousCallback,
    ODEAliasSpecifier,
    ReturnCode,
    remake, successful_retcode,
    addsteps!, savevalues!, terminate!, reinit!,
    du_cache, full_cache, u_cache, get_tmp_cache,
    add_saveat!

# Import from OrdinaryDiffEqCore
using OrdinaryDiffEqCore: OrdinaryDiffEqCore,
    CompositeAlgorithm, AutoSwitch,
    ShampineCollocationInit, BrownFullBasicInit, NoInit,
    isfsal, ode_interpolant,
    IController, PIController, PIDController

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
export DiscreteProblem, DiscreteFunction

# Ensemble
export EnsembleProblem, EnsembleAlgorithm, EnsembleSolution
export EnsembleThreads, EnsembleDistributed, EnsembleSerial, EnsembleSplitThreads

# Callbacks
export CallbackSet, ContinuousCallback, DiscreteCallback, VectorContinuousCallback

# Utilities
export ODEAliasSpecifier, ReturnCode
export remake, successful_retcode, reinit!
export addsteps!, ode_interpolant, terminate!, savevalues!, isfsal
export du_cache, full_cache, u_cache, get_tmp_cache, add_saveat!

# Core types
export CompositeAlgorithm, AutoSwitch
export ShampineCollocationInit, BrownFullBasicInit, NoInit
export IController, PIController, PIDController

# Default algorithm
export DefaultODEAlgorithm, DefaultImplicitODEAlgorithm

# Widely-used algorithms
export Tsit5, AutoTsit5
export Vern6, Vern7, Vern8, Vern9, AutoVern6, AutoVern7, AutoVern8, AutoVern9
export Rosenbrock23, Rodas5P
export FBDF

end # module
