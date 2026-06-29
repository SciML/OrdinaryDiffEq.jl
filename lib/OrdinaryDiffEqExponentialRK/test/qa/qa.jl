using SciMLTesting, OrdinaryDiffEqExponentialRK, Test

run_qa(
    OrdinaryDiffEqExponentialRK;
    explicit_imports = true,
    ei_kwargs = (
        # Internal (non-public) names from sibling solver packages and SciMLBase
        # that have no public replacement yet. Tracked in SciML/OrdinaryDiffEq.jl#3776
        # for upstream `public` declarations.
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internals
                Symbol("@cache"), :ExponentialAlgorithm,
                :OrdinaryDiffEqAdaptiveExponentialAlgorithm,
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqExponentialAlgorithm,
                :OrdinaryDiffEqMutableCache, :alg_cache, :alg_adaptive_order,
                :fsal_typeof, :generic_solver_docstring, :get_fsalfirstlast,
                :isdtchangeable, :ismultistep, :perform_step!, :unwrap_alg,
                :_fixup_ad, :full_cache,
                # SciMLBase internals
                :_unwrap_val, :UDerivativeWrapper, :UJacobianWrapper,
                # DiffEqBase internals
                :calculate_residuals, :calculate_residuals!,
                # OrdinaryDiffEqDifferentiation internals
                :build_jac_config, :calc_J, :calc_J!,
            )),
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internals
                :increment_nf!, :set_EEst!,
                # DiffEqBase internal
                :prepare_alg,
            )),
    ),
)
