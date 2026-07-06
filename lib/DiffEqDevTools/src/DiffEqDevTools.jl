module DiffEqDevTools

import DiffEqBase
using DiffEqBase: ExplicitRKTableau, ImplicitRKTableau
import RecipesBase
using RecipesBase: @recipe, @series
import RecursiveArrayTools
using RecursiveArrayTools: recursive_mean, vecvecapply
using DiffEqNoiseProcess: NoiseGrid, NoiseWrapper
import StructArrays
using StructArrays: StructArray
import SciMLBase
using SciMLBase: AbstractODEAlgorithm, AbstractODEProblem,
    AbstractSDEProblem, AbstractEnsembleProblem, AbstractDAEProblem,
    AbstractDEAlgorithm, AbstractTimeseriesSolution,
    DAEProblem, EnsembleProblem, EnsembleSolution, EnsembleThreads,
    NonlinearProblem, ODEProblem, ReturnCode, SDDEProblem, SDEProblem, remake
using CommonSolve: init, solve, step!
import SimpleNonlinearSolve
using SimpleNonlinearSolve: SimpleTrustRegion, AutoFiniteDiff
import LinearAlgebra
import RootedTrees
using RootedTrees: RootedTreeIterator, RungeKuttaMethod, residual_order_condition
import Distributed
import Statistics
using Statistics: mean, std

import Base: length

# These problem/solution/algorithm abstracts and `@def` are owned by SciMLBase but
# not yet declared `public` there; accessed via SciMLBase (their owner).
using SciMLBase: AbstractDDEAlgorithm, AbstractODESolution, AbstractRODEProblem,
    AbstractSDDEProblem, AbstractBVProblem, @def
# `ConvergenceSetup` and `ODERKTableau` are defined and owned only in DiffEqBase
# (not re-exported by SciMLBase) and are not declared `public` there.
using DiffEqBase: ConvergenceSetup, ODERKTableau

import LinearAlgebra: norm, I

const TIMESERIES_ERRORS = Set([:l2, :l∞, :L2, :L∞])
const DENSE_ERRORS = Set([:L2, :L∞])
const WEAK_TIMESERIES_ERRORS = Set([:weak_l2, :weak_l∞])
const WEAK_DENSE_ERRORS = Set([:weak_L2, :weak_L∞])
const WEAK_ERRORS = union(
    Set([:weak_final]),
    WEAK_TIMESERIES_ERRORS, WEAK_DENSE_ERRORS
)
const ALL_ERRORS = union(
    [:final],
    TIMESERIES_ERRORS, DENSE_ERRORS, WEAK_TIMESERIES_ERRORS, WEAK_DENSE_ERRORS, WEAK_ERRORS
)

include("benchmark.jl")
include("convergence.jl")
include("plotrecipes.jl")
include("test_solution.jl")
include("ode_tableaus.jl")
include("tableau_info.jl")

export ConvergenceSimulation, Shootout, ShootoutSet, TestSolution

#Benchmark Functions
export Shootout, ShootoutSet, WorkPrecision, WorkPrecisionSet

export test_convergence, analyticless_test_convergence, appxtrue

export get_sample_errors

#Tab Functions
export stability_region, residual_order_condition, check_tableau, imaginary_stability_interval

export deduce_Butcher_tableau

end # module
