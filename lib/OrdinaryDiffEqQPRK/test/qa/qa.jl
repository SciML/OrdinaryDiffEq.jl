using SciMLTesting, OrdinaryDiffEqQPRK, Test

run_qa(
    OrdinaryDiffEqQPRK;
    explicit_imports = true,
    ei_kwargs = (
        # OrdinaryDiffEqCore solver-interface internals (no public owner elsewhere) plus
        # DiffEqBase residual helpers that QPRK extends but are not yet marked public.
        all_explicit_imports_are_public = (
            ignore = (
                :OrdinaryDiffEqAdaptiveAlgorithm, :OrdinaryDiffEqConstantCache,
                :OrdinaryDiffEqMutableCache, :explicit_rk_docstring, :alg_cache,
                :trivial_limiter!, :perform_step!, :get_fsalfirstlast, :constvalue,
                Symbol("@cache"), Symbol("@fold"), Symbol("@OnDemandTableauExtract"),
                :calculate_residuals, :calculate_residuals!,
            ),
        ),
        # OrdinaryDiffEqCore step-statistics setters accessed by qualified name.
        all_qualified_accesses_are_public = (
            ignore = (:increment_nf!, :set_EEst!),
        ),
    ),
)
