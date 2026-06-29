using SciMLTesting, OrdinaryDiffEqHighOrderRK, Test

run_qa(
    OrdinaryDiffEqHighOrderRK;
    explicit_imports = true,
    ei_kwargs = (
        # `@reexport using SciMLBase` is the intended public-API reexport; it
        # necessarily brings the `SciMLBase`/`Reexport`/`@reexport` bindings in
        # implicitly, which there is no way to make explicit.
        no_implicit_imports = (ignore = (:SciMLBase, :Reexport, Symbol("@reexport")),),
        # `@tight_loop_macros` is used throughout the perform_step/addsteps loops
        # but ExplicitImports cannot see macro invocations nested inside the
        # `@muladd`-expanded function bodies, so it is a false stale positive.
        no_stale_explicit_imports = (ignore = (Symbol("@tight_loop_macros"),),),
        # OrdinaryDiffEqCore-internal solver-statistics accessors; not yet
        # declared public upstream (SciML/OrdinaryDiffEq.jl#3776).
        all_qualified_accesses_are_public = (
            ignore = (:get_EEst, :increment_nf!, :set_EEst!),
        ),
        # Internal solver-extension API of OrdinaryDiffEqCore (and the
        # `calculate_residuals`/`@tight_loop_macros` helpers of DiffEqBase) that
        # is not yet declared public upstream (SciML/OrdinaryDiffEq.jl#3776).
        all_explicit_imports_are_public = (
            ignore = (
                Symbol("@cache"), Symbol("@tight_loop_macros"), :CompiledFloats,
                :DerivativeOrderNotPossibleError, :OrdinaryDiffEqAdaptiveAlgorithm,
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
                :_ode_addsteps!, :_ode_interpolant, :_ode_interpolant!, :alg_cache,
                :beta1_default, :beta2_default, :calculate_residuals,
                :calculate_residuals!, :constvalue, :explicit_rk_docstring,
                :get_fsalfirstlast, :isdp8, :perform_step!, :qmax_default,
                :qmin_default, :trivial_limiter!, :unwrap_alg,
            ),
        ),
    ),
)
