module DiffEqBase
if isdefined(Base, :Experimental) &&
        isdefined(Base.Experimental, Symbol("@max_methods"))
    @eval Base.Experimental.@max_methods 1
end

import PrecompileTools

import FastPower
@deprecate fastpow(x, y) FastPower.fastpower(x, y)

using ArrayInterface

using StaticArraysCore # data arrays

using LinearAlgebra, Printf

using DocStringExtensions

using FunctionWrappers: FunctionWrapper

using MuladdMacro


using FastBroadcast: @.., True, False

using Static: reduce_tup

import RecursiveArrayTools
import TruncatedStacktraces

using Setfield


using Markdown

using ConcreteStructs: @concrete
using FastClosures: @closure

import FunctionWrappersWrappers

using SciMLBase

using SciMLOperators: AbstractSciMLOperator, AbstractSciMLScalarOperator, DEFAULT_UPDATE_FUNC

using SciMLBase: @def, DEIntegrator, AbstractDEProblem,
    AbstractDiffEqInterpolation,
    DECallback, AbstractDEOptions, DECache, AbstractContinuousCallback,
    AbstractDiscreteCallback, AbstractLinearProblem,
    AbstractNonlinearProblem,
    AbstractOptimizationProblem, AbstractSteadyStateProblem,
    AbstractJumpProblem,
    AbstractNoiseProblem, AbstractEnsembleProblem,
    AbstractDynamicalODEProblem,
    AbstractDEAlgorithm, StandardODEProblem, AbstractIntegralProblem,
    AbstractSensitivityAlgorithm, AbstractODEAlgorithm,
    AbstractSDEAlgorithm, AbstractDDEAlgorithm, AbstractDAEAlgorithm,
    AbstractSDDEAlgorithm, AbstractRODEAlgorithm, AbstractBVPAlgorithm,
    DAEInitializationAlgorithm,
    AbstractSteadyStateAlgorithm, AbstractODEProblem,
    AbstractDiscreteProblem, AbstractNonlinearAlgorithm,
    AbstractSDEProblem, AbstractRODEProblem, AbstractDDEProblem,
    AbstractDAEProblem, AbstractSDDEProblem, AbstractBVProblem,
    AbstractTimeseriesSolution, AbstractNoTimeSolution, numargs,
    AbstractODEFunction, AbstractSDEFunction, AbstractRODEFunction,
    AbstractDDEFunction, AbstractSDDEFunction, AbstractDAEFunction,
    AbstractNonlinearFunction, AbstractEnsembleSolution,
    AbstractODESolution, AbstractRODESolution, AbstractDAESolution,
    AbstractDDESolution,
    EnsembleAlgorithm, EnsembleSolution, EnsembleSummary,
    NonlinearSolution,
    TimeGradientWrapper, TimeDerivativeWrapper, UDerivativeWrapper,
    UJacobianWrapper, ParamJacobianWrapper, JacobianWrapper,
    check_error!, has_jac, has_tgrad, has_Wfact, has_Wfact_t, has_paramjac,
    AbstractODEIntegrator, AbstractSDEIntegrator, AbstractRODEIntegrator,
    AbstractDDEIntegrator, AbstractSDDEIntegrator,
    AbstractDAEIntegrator, unwrap_cache, has_reinit, reinit!,
    postamble!, last_step_failed, islinear, has_stats,
    initialize_dae!, build_solution, solution_new_retcode,
    solution_new_tslocation, plot_indices, NonlinearAliasSpecifier,
    NullParameters, isinplace, AbstractADType, AbstractDiscretization,
    DISCRETE_OUTOFPLACE_DEFAULT, DISCRETE_INPLACE_DEFAULT,
    has_analytic, calculate_solution_errors!, AbstractNoiseProcess,
    has_colorvec, parameterless_type, undefined_exports,
    is_diagonal_noise, AbstractDiffEqFunction, sensitivity_solution,
    interp_summary, AbstractHistoryFunction, LinearInterpolation,
    ConstantInterpolation, HermiteInterpolation, SensitivityInterpolation,
    NoAD, @add_kwonly,
    calculate_ensemble_errors, isconstant,
    DEFAULT_REDUCTION, isautodifferentiable,
    isadaptive, isdiscrete, has_syms, AbstractAnalyticalSolution,
    RECOMPILE_BY_DEFAULT, wrap_sol, has_destats

import SciMLBase: solve, init, step!, solve!, __init, __solve, update_coefficients!,
    update_coefficients, isadaptive, wrapfun_oop, wrapfun_iip,
    unwrap_fw, promote_tspan, set_u!, set_t!, set_ut!,
    extract_alg, checkkwargs, has_kwargs, _concrete_solve_adjoint, _concrete_solve_forward,
    eltypedual, get_updated_symbolic_problem, get_concrete_p, get_concrete_u0, promote_u0,
    isconcreteu0, isconcretedu0, get_concrete_du0, _reshape, value, unitfulvalue, anyeltypedual, allowedkeywords,
    sse, totallength, __sum, DualEltypeChecker, KeywordArgError, KeywordArgWarn, KeywordArgSilent, KWARGWARN_MESSAGE, KWARGERROR_MESSAGE,
    CommonKwargError, IncompatibleInitialConditionError, NO_DEFAULT_ALGORITHM_MESSAGE, NoDefaultAlgorithmError, NO_TSPAN_MESSAGE, NoTspanError,
    NAN_TSPAN_MESSAGE, NaNTspanError, NON_SOLVER_MESSAGE, NonSolverError, NOISE_SIZE_MESSAGE, NoiseSizeIncompatibilityError, PROBSOLVER_PAIRING_MESSAGE,
    ProblemSolverPairingError, compatible_problem_types, DIRECT_AUTODIFF_INCOMPATIBILITY_MESSAGE, DirectAutodiffError, NONNUMBER_ELTYPE_MESSAGE, NonNumberEltypeError,
    GENERIC_NUMBER_TYPE_ERROR_MESSAGE, GenericNumberTypeError, COMPLEX_SUPPORT_ERROR_MESSAGE, ComplexSupportError, COMPLEX_TSPAN_ERROR_MESSAGE, ComplexTspanError,
    TUPLE_STATE_ERROR_MESSAGE, TupleStateError, MASS_MATRIX_ERROR_MESSAGE, IncompatibleMassMatrixError, LATE_BINDING_TSTOPS_ERROR_MESSAGE, LateBindingTstopsNotSupportedError,
    NONCONCRETE_ELTYPE_MESSAGE, NonConcreteEltypeError, _vec

import SciMLStructures

using Reexport
Reexport.@reexport using SciMLBase

SciMLBase.isfunctionwrapper(x::FunctionWrapper) = true

# Rootfinder for callbacks
using BracketingNonlinearSolve: ModAB

import SymbolicIndexingInterface as SII

## Extension Functions

## Types

"""
$(TYPEDEF)
"""
abstract type Tableau end

"""
$(TYPEDEF)
"""
abstract type ODERKTableau <: Tableau end

"""
$(TYPEDEF)
"""
abstract type DECostFunction end

import SciMLBase: Void, unwrapped_f

include("utils.jl")
include("stats.jl")
include("calculate_residuals.jl")
include("tableaus.jl")
include("dae_initialization.jl")

include("callbacks.jl")
include("common_defaults.jl")
include("solve.jl")
include("internal_euler.jl")
include("norecompile.jl")
include("integrator_accessors.jl")

# This is only used for oop stiff solvers
default_factorize(A) = lu(A; check = false)

if isdefined(SciMLBase, :AbstractParameterizedFunction)
    import SciMLBase: AbstractParameterizedFunction
else
    """
    $(TYPEDEF)
    """
    abstract type AbstractParameterizedFunction{iip} <: AbstractODEFunction{iip} end
end

"""
$(TYPEDEF)
"""
struct ConvergenceSetup{P, C}
    probs::P
    convergence_axis::C
end

export initialize!, finalize!

export SensitivityADPassThrough

include("precompilation.jl")

end # module
