using SciMLTesting, OrdinaryDiffEqFeagin, Test

run_qa(
    OrdinaryDiffEqFeagin;
    explicit_imports = true,
    ei_kwargs = (;
        # `@reexport using SciMLBase` deliberately re-exports SciMLBase's API
        # as part of this solver package's public surface.
        no_implicit_imports = (; ignore = (:SciMLBase,)),
        # OrdinaryDiffEqCore solver-interface internals not declared `public`.
        all_qualified_accesses_are_public = (;
            ignore = (:increment_nf!, :set_EEst!)),
        # Non-public solver-interface internals (OrdinaryDiffEqCore) and
        # DiffEqBase residual/loop internals; not declared `public` upstream.
        all_explicit_imports_are_public = (;
            ignore = (
                Symbol("@cache"), Symbol("@tight_loop_macros"),
                :CompiledFloats, :OrdinaryDiffEqAdaptiveAlgorithm,
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
                :alg_cache, :calculate_residuals, :calculate_residuals!,
                :constvalue, :generic_solver_docstring, :get_fsalfirstlast,
                :perform_step!, :trivial_limiter!,
            )),
    ),
)
