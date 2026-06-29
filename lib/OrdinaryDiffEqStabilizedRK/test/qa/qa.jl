using SciMLTesting, OrdinaryDiffEqStabilizedRK, Test

run_qa(
    OrdinaryDiffEqStabilizedRK;
    explicit_imports = true,
    ei_kwargs = (
        # OrdinaryDiffEqCore solver-interface internals (no public owner elsewhere)
        # plus SciMLBase/DiffEqBase helpers that StabilizedRK extends/uses but which
        # are not yet marked public in their owner modules.
        all_explicit_imports_are_public = (
            ignore = (
                :OrdinaryDiffEqAdaptiveAlgorithm, :OrdinaryDiffEqConstantCache,
                :OrdinaryDiffEqMutableCache, :PredictiveController,
                :alg_adaptive_order, :alg_cache, :alg_extrapolates, :constvalue,
                :default_controller, :fac_default_gamma, :gamma_default,
                :generic_solver_docstring, :get_fsalfirstlast,
                :has_dtnew_modification, :perform_step!, :qmax_default, :unwrap_alg,
                Symbol("@cache"),
                :_vec, :value, :calculate_residuals, :calculate_residuals!,
            ),
        ),
        # OrdinaryDiffEqCore step-statistics setters accessed by qualified name.
        all_qualified_accesses_are_public = (
            ignore = (:increment_nf!, :set_EEst!),
        ),
    ),
)
