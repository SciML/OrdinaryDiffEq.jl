using SciMLTesting, OrdinaryDiffEqLinear, Test

run_qa(
    OrdinaryDiffEqLinear;
    explicit_imports = true,
    ei_kwargs = (;
        # All residual names below are non-public internals of upstream packages
        # that have no public-API replacement yet; they are tightly listed (not
        # blanket ei_broken) and feed the next make-public round (see SciML/OrdinaryDiffEq.jl#3776).
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internal solver-interface API (no public surface yet)
                Symbol("@cache"), :OrdinaryDiffEqAdaptiveAlgorithm,
                :OrdinaryDiffEqAlgorithm, :OrdinaryDiffEqConstantCache,
                :OrdinaryDiffEqLinearExponentialAlgorithm, :OrdinaryDiffEqMutableCache,
                :alg_cache, :alg_extrapolates, :dt_required, :generic_solver_docstring,
                :get_fsalfirstlast, :isdtchangeable, :perform_step!, :unwrap_alg,
                # DiffEqBase internal (owner of calculate_residuals!), non-public
                :calculate_residuals!,
                # SciMLBase internal (owner of _vec), non-public
                :_vec,
                # SciMLOperators internal abstract type, non-public
                :AbstractSciMLOperator,
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                :alloc_mem,          # ExponentialUtilities internal cache allocator
                :increment_nf!,      # OrdinaryDiffEqCore internal stats mutator
                :set_EEst!,          # OrdinaryDiffEqCore internal error-estimate setter
                :prepare_alg,        # DiffEqBase internal algorithm-preparation hook
            ),
        ),
    ),
)
