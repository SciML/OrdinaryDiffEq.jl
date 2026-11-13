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
        # These names are not used by DiffEqBase itself, but every one is
        # reached through `DiffEqBase.X` or `import DiffEqBase: X` by its
        # extensions, downstream sublibraries and packages (e.g.
        # DiffEqNoiseProcess calls `DiffEqBase.has_reinit`, Sundials calls
        # `DiffEqBase.update_coefficients!`, OrdinaryDiffEq imports
        # `DiffEqBase: DEFAULT_UPDATE_FUNC`), or test files — verified by
        # grepping every lib/*/src, lib/*/ext, src/, ext/ and all of
        # ~/.julia/packages. ExplicitImports cannot see that cross-module
        # usage, so it reports them as stale; they must remain part of this
        # package's namespace contract.
        no_stale_explicit_imports = (;
            ignore = (
                Symbol("@add_kwonly"), Symbol("@def"),
                :AbstractDAEProblem, :AbstractDAESolution,
                :AbstractDDEFunction, :AbstractDDEIntegrator,
                :AbstractDDEProblem,
                :AbstractDiffEqFunction, :AbstractDiffEqInterpolation,
                :AbstractDiscreteProblem,
                :AbstractDynamicalODEProblem, :AbstractEnsembleSolution,
                :AbstractHistoryFunction, :AbstractNoTimeSolution,
                :AbstractNoiseProcess,
                :AbstractODEProblem, :AbstractODESolution,
                :AbstractRODEAlgorithm, :AbstractRODEIntegrator,
                :AbstractRODEProblem, :AbstractRODESolution,
                :AbstractSDDEIntegrator, :AbstractSDDEProblem,
                :AbstractSDEIntegrator,
                :AbstractSDEProblem,
                :AbstractSensitivityAlgorithm, :AbstractTimeseriesSolution,
                :ConstantInterpolation, :DECache,
                :DEFAULT_UPDATE_FUNC,
                :DISCRETE_INPLACE_DEFAULT, :DISCRETE_OUTOFPLACE_DEFAULT,
                :DualEltypeChecker,
                :HermiteInterpolation,
                :LinearInterpolation,
                :NullParameters,
                :ParamJacobianWrapper,
                :SensitivityInterpolation, :StandardODEProblem,
                :TimeDerivativeWrapper,
                :TimeGradientWrapper, :UDerivativeWrapper,
                :UJacobianWrapper, :__sum, :_reshape, :_vec, :allowedkeywords,
                :anyeltypedual, :calculate_ensemble_errors,
                :calculate_solution_errors!, :check_error!,
                :has_Wfact, :has_Wfact_t,
                :has_analytic, :has_colorvec, :has_jac, :has_paramjac,
                :has_reinit, :has_stats, :has_tgrad,
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
                # machinery, plus the error types it throws; many are
                # additionally part of the downstream `DiffEqBase.X` namespace
                # contract (see above).
                :__sum, :_reshape, :_vec, :allowedkeywords, :anyeltypedual,
                :checkkwargs, :eltypedual,
                :extract_alg, :get_concrete_du0, :has_Wfact, :has_Wfact_t,
                :has_colorvec, :has_kwargs, :isconcretedu0,
                :plot_indices, :solution_new_tslocation, :sse, :totallength,
                :undefined_exports, :unwrap_cache, :DualEltypeChecker,
                :ComplexSupportError, :ComplexTspanError,
                :DISCRETE_INPLACE_DEFAULT, :DISCRETE_OUTOFPLACE_DEFAULT,
                :DirectAutodiffError,
                :GenericNumberTypeError,
                :LateBindingTstopsNotSupportedError,
                :NaNTspanError, :NoDefaultAlgorithmError,
                :NoTspanError, :NoiseSizeIncompatibilityError,
                :NonConcreteEltypeError, :NonNumberEltypeError, :NonSolverError,
                :ProblemSolverPairingError,
                # SciMLOperators — owner-internal, imported downstream via
                # `using DiffEqBase: DEFAULT_UPDATE_FUNC`
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
