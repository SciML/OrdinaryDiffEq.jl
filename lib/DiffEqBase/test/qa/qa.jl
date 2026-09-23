using SciMLTesting, DiffEqBase, Aqua, Test

run_qa(
    DiffEqBase;
    # `@reexport using SciMLBase` deliberately reexports SciMLBase's public API.
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    # `piracies` and the two sub-checks with hand-rolled bounds below are run
    # manually: the ambiguity check keeps the existing <=4 bound on ambiguities
    # located in DiffEqBase, and the unbound-args check keeps the direct
    # `detect_unbound_args_recursively` assertion.
    aqua_kwargs = (; piracies = false, ambiguities = false, unbound_args = false),
    ei_kwargs = (;
        # Most of these names are not used by DiffEqBase itself, but are
        # imported into its namespace and reached through `DiffEqBase.X` by its
        # extensions, the InternalEuler submodule, downstream sublibraries and
        # packages (e.g. DiffEqNoiseProcess calls `DiffEqBase.has_reinit`,
        # Sundials calls `DiffEqBase.update_coefficients!`), and test files.
        # ExplicitImports cannot see that cross-module usage, so it reports
        # them as stale; they must remain part of this package's namespace
        # contract.
        no_stale_explicit_imports = (;
            ignore = (
                Symbol("@add_kwonly"), Symbol("@def"),
                :AbstractAnalyticalSolution, :AbstractDAEFunction,
                :AbstractDAEIntegrator, :AbstractDAEProblem,
                :AbstractDAESolution, :AbstractDDEFunction,
                :AbstractDDEIntegrator, :AbstractDDEProblem,
                :AbstractDDESolution, :AbstractDEOptions,
                :AbstractDiffEqFunction, :AbstractDiffEqInterpolation,
                :AbstractDiscreteProblem, :AbstractDiscretization,
                :AbstractDynamicalODEProblem, :AbstractEnsembleSolution,
                :AbstractHistoryFunction, :AbstractNoTimeSolution,
                :AbstractNoiseProcess, :AbstractNonlinearFunction,
                :AbstractODEProblem, :AbstractODESolution,
                :AbstractOptimizationProblem, :AbstractRODEAlgorithm,
                :AbstractRODEFunction, :AbstractRODEIntegrator,
                :AbstractRODEProblem, :AbstractRODESolution,
                :AbstractSDDEAlgorithm, :AbstractSDDEFunction,
                :AbstractSDDEIntegrator, :AbstractSDDEProblem,
                :AbstractSDEFunction, :AbstractSDEIntegrator,
                :AbstractSDEProblem, :AbstractSciMLScalarOperator,
                :AbstractSensitivityAlgorithm, :AbstractTimeseriesSolution,
                :COMPLEX_SUPPORT_ERROR_MESSAGE, :COMPLEX_TSPAN_ERROR_MESSAGE,
                :CommonKwargError, :ConstantInterpolation, :DECache,
                :DEFAULT_REDUCTION, :DEFAULT_UPDATE_FUNC,
                :DIRECT_AUTODIFF_INCOMPATIBILITY_MESSAGE,
                :DISCRETE_INPLACE_DEFAULT, :DISCRETE_OUTOFPLACE_DEFAULT,
                :DualEltypeChecker, :EnsembleAlgorithm,
                :GENERIC_NUMBER_TYPE_ERROR_MESSAGE, :HermiteInterpolation,
                :IncompatibleInitialConditionError,
                :IncompatibleMassMatrixError, :JacobianWrapper,
                :KWARGERROR_MESSAGE, :KWARGWARN_MESSAGE, :KeywordArgError,
                :KeywordArgSilent, :KeywordArgWarn,
                :LATE_BINDING_TSTOPS_ERROR_MESSAGE, :LinearInterpolation,
                :MASS_MATRIX_ERROR_MESSAGE, :NAN_TSPAN_MESSAGE,
                :NOISE_SIZE_MESSAGE, :NONCONCRETE_ELTYPE_MESSAGE,
                :NONNUMBER_ELTYPE_MESSAGE, :NON_SOLVER_MESSAGE,
                :NO_DEFAULT_ALGORITHM_MESSAGE, :NO_TSPAN_MESSAGE, :NoAD,
                :NonlinearAliasSpecifier, :NullParameters,
                :PROBSOLVER_PAIRING_MESSAGE, :ParamJacobianWrapper,
                :SensitivityInterpolation, :StandardODEProblem,
                :TUPLE_STATE_ERROR_MESSAGE, :TimeDerivativeWrapper,
                :TimeGradientWrapper, :TupleStateError, :UDerivativeWrapper,
                :UJacobianWrapper, :__sum, :_reshape, :_vec, :allowedkeywords,
                :anyeltypedual, :calculate_ensemble_errors,
                :calculate_solution_errors!, :check_error!,
                :compatible_problem_types, :has_Wfact, :has_Wfact_t,
                :has_analytic, :has_colorvec, :has_jac, :has_paramjac,
                :has_reinit, :has_stats, :has_syms, :has_tgrad,
                :initialize_dae!, :interp_summary, :is_diagonal_noise,
                :isautodifferentiable, :numargs, :parameterless_type,
                :plot_indices, :sensitivity_solution, :set_ut!,
                :solution_new_tslocation, :sse, :totallength,
                :undefined_exports, :unitfulvalue, :unwrap_cache,
                # Referenced unqualified inside the `@verbosity_specifier`
                # macro expansion in src/verbosity.jl, which ExplicitImports
                # cannot see.
                :AbstractVerbosityPreset, :AbstractVerbositySpecifier,
                :MessageLevel, :Standard,
            ),
        ),
        # Internal (non-`public`) names of upstream packages that DiffEqBase
        # genuinely needs and that have no public replacement yet.
        all_qualified_accesses_are_public = (;
            ignore = (
                # Base / Base.Experimental / Base.FastMath / Base.Iterators internals
                :Experimental, Symbol("@max_methods"), :_nt_names, :diff_names,
                :promote_op, :structdiff, :FastMath, :sqrt_fast, :Zip,
                # StaticArraysCore — owner-internal, no public alternative
                :StaticArray,
                # SciMLBase internals with no public replacement yet
                :diagnose_symbolic_instability, :has_mtk_sys,
                :log_numerical_instability, :report_integrator_failure,
            ),
        ),
        # Internal (non-`public`) names imported from upstream packages.
        all_explicit_imports_are_public = (;
            ignore = (
                # SciMLBase internals genuinely needed by solve/remake
                # machinery, plus the error types/messages it throws; many are
                # additionally part of the downstream `DiffEqBase.X` namespace
                # contract (see above).
                :__sum, :_reshape, :_vec, :allowedkeywords, :anyeltypedual,
                :checkkwargs, :compatible_problem_types, :eltypedual,
                :extract_alg, :get_concrete_du0, :has_Wfact, :has_Wfact_t,
                :has_colorvec, :has_kwargs, :has_syms, :isconcretedu0,
                :plot_indices, :solution_new_tslocation, :sse, :totallength,
                :undefined_exports, :unwrap_cache, :DualEltypeChecker,
                :COMPLEX_SUPPORT_ERROR_MESSAGE, :COMPLEX_TSPAN_ERROR_MESSAGE,
                :CommonKwargError, :ComplexSupportError, :ComplexTspanError,
                :DIRECT_AUTODIFF_INCOMPATIBILITY_MESSAGE,
                :DISCRETE_INPLACE_DEFAULT, :DISCRETE_OUTOFPLACE_DEFAULT,
                :DirectAutodiffError, :GENERIC_NUMBER_TYPE_ERROR_MESSAGE,
                :GenericNumberTypeError, :IncompatibleInitialConditionError,
                :IncompatibleMassMatrixError, :KWARGERROR_MESSAGE,
                :KWARGWARN_MESSAGE, :KeywordArgSilent, :KeywordArgWarn,
                :LATE_BINDING_TSTOPS_ERROR_MESSAGE,
                :LateBindingTstopsNotSupportedError, :MASS_MATRIX_ERROR_MESSAGE,
                :NAN_TSPAN_MESSAGE, :NOISE_SIZE_MESSAGE,
                :NONCONCRETE_ELTYPE_MESSAGE, :NONNUMBER_ELTYPE_MESSAGE,
                :NON_SOLVER_MESSAGE, :NO_DEFAULT_ALGORITHM_MESSAGE,
                :NO_TSPAN_MESSAGE, :NaNTspanError, :NoDefaultAlgorithmError,
                :NoTspanError, :NoiseSizeIncompatibilityError,
                :NonConcreteEltypeError, :NonNumberEltypeError, :NonSolverError,
                :PROBSOLVER_PAIRING_MESSAGE, :ProblemSolverPairingError,
                :TUPLE_STATE_ERROR_MESSAGE, :TupleStateError,
                # SciMLOperators — owner-internal, no public alternative
                :DEFAULT_UPDATE_FUNC,
                # FunctionWrappers — owner-internal, no public alternative
                :FunctionWrapper,
            ),
        ),
    ),
)

@testset "Aqua tests (performance)" begin
    # This tests that we don't accidentally run into
    # https://github.com/JuliaLang/julia/issues/29393
    # Aqua.test_unbound_args(DiffEqBase) # fails
    ua = Aqua.detect_unbound_args_recursively(DiffEqBase)
    @test length(ua) == 0
    # Uncomment for debugging:
    # @show ua

    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
    # Test that we're not introducing method ambiguities across deps
    ambs = Aqua.detect_ambiguities(DiffEqBase; recursive = true)
    pkg_match(pkgname, pkdir::Nothing) = false
    pkg_match(pkgname, pkdir::AbstractString) = occursin(pkgname, pkdir)
    filter!(x -> pkg_match("DiffEqBase", pkgdir(last(x).module)), ambs)

    # Uncomment for debugging:
    # for method_ambiguity in ambs
    #     @show method_ambiguity
    # end
    @warn "Number of method ambiguities: $(length(ambs))"
    @test length(ambs) ≤ 4
end
