module DiffEqDevTools

using DiffEqBase: AbstractODEAlgorithm
using DiffEqBase, RecipesBase, RecursiveArrayTools, DiffEqNoiseProcess, StructArrays
using NLsolve, LinearAlgebra, RootedTrees

using LinearAlgebra, Distributed

using Statistics

import Base: length

import DiffEqBase: AbstractODEProblem, AbstractDDEProblem, AbstractDDEAlgorithm,
    AbstractODESolution, AbstractRODEProblem, AbstractSDEProblem,
    AbstractSDDEProblem, AbstractEnsembleProblem,
    AbstractDAEProblem, AbstractBVProblem, @def, ConvergenceSetup,
    AbstractDEAlgorithm,
    ODERKTableau, AbstractTimeseriesSolution, ExplicitRKTableau,
    ImplicitRKTableau

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
