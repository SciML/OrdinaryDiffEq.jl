using SciMLTesting, OrdinaryDiffEqSSPRK, Test

run_qa(
    OrdinaryDiffEqSSPRK;
    explicit_imports = true,
    ei_kwargs = (
        # Solver-internal names with no public counterpart in their owner.
        # OrdinaryDiffEqCore internals: the RK perform_step/cache/limiter API
        # (perform_step!, alg_cache, @cache, isfsal, constvalue, trivial_limiter!,
        #  explicit_rk_docstring, get_fsalfirstlast, ssp_coefficient, the
        #  _ode_interpolant/_ode_addsteps interpolation hooks, the abstract
        #  algorithm/cache supertypes). DiffEqBase residual hooks
        # (calculate_residuals[!]). SciMLBase macros (@def, @cache reexport).
        # See SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                Symbol("@cache"), Symbol("@def"),
                :OrdinaryDiffEqAlgorithm, :OrdinaryDiffEqMutableCache,
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqAdaptiveAlgorithm,
                :perform_step!, :alg_cache, :isfsal, :constvalue,
                :trivial_limiter!, :explicit_rk_docstring, :get_fsalfirstlast,
                :ssp_coefficient, :_ode_interpolant, :_ode_interpolant!,
                :_ode_addsteps!,
                :calculate_residuals, :calculate_residuals!,
            ),
        ),
        # Qualified accesses to solver-internal names: OrdinaryDiffEqCore
        # step bookkeeping/test problems (increment_nf!, set_EEst!, lorenz[_oop])
        # and SciMLBase specialization tags used only inside the precompile
        # workload (FunctionWrapperSpecialize, NoSpecialize).
        all_qualified_accesses_are_public = (;
            ignore = (
                :increment_nf!, :set_EEst!, :lorenz, :lorenz_oop,
                :FunctionWrapperSpecialize, :NoSpecialize,
            ),
        ),
    ),
)
