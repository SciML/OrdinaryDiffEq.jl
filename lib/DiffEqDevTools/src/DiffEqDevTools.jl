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
                   DEAlgorithm,
                   ODERKTableau, AbstractTimeseriesSolution, ExplicitRKTableau,
                   ImplicitRKTableau

import LinearAlgebra: norm, I

const TIMESERIES_ERRORS = Set([:l2, :l∞, :L2, :L∞])
const DENSE_ERRORS = Set([:L2, :L∞])
const WEAK_TIMESERIES_ERRORS = Set([:weak_l2, :weak_l∞])
const WEAK_DENSE_ERRORS = Set([:weak_L2, :weak_L∞])
const WEAK_ERRORS = union(Set([:weak_final]),
    WEAK_TIMESERIES_ERRORS, WEAK_DENSE_ERRORS)
const ALL_ERRORS = union([:final],
    TIMESERIES_ERRORS, DENSE_ERRORS, WEAK_TIMESERIES_ERRORS, WEAK_DENSE_ERRORS, WEAK_ERRORS)

include("benchmark.jl")
include("convergence.jl")
include("plotrecipes.jl")
include("test_solution.jl")
include("ode_tableaus.jl")
include("tableau_info.jl")

export ConvergenceSimulation, Shootout, ShootoutSet, TestSolution

#Benchmark Functions
export Shootout, ShootoutSet, WorkPrecision, WorkPrecisionSet

export test_convergence, analyticless_test_convergence, appxtrue!, appxtrue

export get_sample_errors

#Tab Functions
export stability_region, residual_order_condition, check_tableau, imaginary_stability_interval

#Tableaus
export deduce_Butcher_tableau
export constructEuler, constructKutta3, constructRK4, constructRK438Rule,
       constructImplicitEuler, constructMidpointRule, constructTrapezoidalRule,
       constructLobattoIIIA4, constructLobattoIIIB2, constructLobattoIIIB4,
       constructLobattoIIIC2, constructLobattoIIIC4, constructLobattoIIICStar2,
       constructLobattoIIICStar4, constructLobattoIIID2, constructLobattoIIID4,
       constructRadauIA3, constructRadauIA5,
       constructRadauIIA3, constructRadauIIA5,
       constructRalston, constructRalston4, constructHeun, constructRKF5,
       constructBogakiShampine3,
       constructCashKarp, constructRKF8, constructDormandPrince8,
       constructFeagin10, constructFeagin12, constructFeagin14,
       constructDormandPrince8_64bit, constructRKF5, constructRungeFirst5,
       constructCassity5, constructLawson5,
       constructLutherKonen5, constructLutherKonen52,
       constructLutherKonen53, constructPapakostasPapaGeorgiou5,
       constructPapakostasPapaGeorgiou52,
       constructTsitouras5, constructBogakiShampine5, constructSharpSmart5,
       constructButcher6, constructButcher7,
       constructDverk, constructClassicVerner6,
       constructClassicVerner7, constructClassicVerner8,
       constructVernerRobust7, constructEnrightVerner7,
       constructTanakaYamashitaStable7,
       constructTanakaYamashitaEfficient7, constructSharpSmart7,
       constructSharpVerner7,
       constructVernerEfficient7, constructCooperVerner8,
       constructCooperVerner82,
       constructTsitourasPapakostas8, constructdverk78, constructEnrightVerner8,
       constructCurtis8, constructVernerRobust9, constructVernerEfficient9,
       constructSharp9, constructTsitouras9,
       constructTsitouras92,
       constructOno12, constructCurtis10,
       constructOno10,
       constructCurtis10, constructBaker10,
       constructHairer10, constructButcher63,
       constructButcher6, constructButcher62,
       constructVerner6, constructDormandPrince6,
       constructSharpVerner6, constructVerner9162,
       constructVerner916, constructVernerRobust6,
       constructVernerEfficient6, constructPapakostas6, constructLawson6,
       constructTsitourasPapakostas6, constructDormandLockyerMcCorriganPrince6,
       constructTanakaKasugaYamashitaYazaki6D,
       constructTanakaKasugaYamashitaYazaki6C,
       constructTanakaKasugaYamashitaYazaki6B,
       constructTanakaKasugaYamashitaYazaki6A,
       constructMikkawyEisa, constructChummund6, constructChummund62,
       constructHuta62, constructHuta6, constructRKF4,
       constructVerner7, constructVerner8,
       constructVerner6, constructSSPRK22, constructSSPRK33,
       constructSSPRK43, constructSSPRK104,
       constructRKO65

end # module
