using SciMLTesting, OrdinaryDiffEqExplicitRK, Test

run_qa(
    OrdinaryDiffEqExplicitRK;
    explicit_imports = true,
    ei_kwargs = (;
        # Names imported from OrdinaryDiffEqCore/DiffEqBase/TruncatedStacktraces that are
        # genuinely internal (not exported, not `public`) in their owner module. These are
        # the solver-internal extension points OrdinaryDiffEqExplicitRK must hook into;
        # there is no public alternative. Tracked for make-public in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                :alg_adaptive_order, :alg_stability_size, :alg_cache, :unwrap_alg,
                :isfsal, :perform_step!, :get_fsalfirstlast,
                :OrdinaryDiffEqAdaptiveAlgorithm, :OrdinaryDiffEqConstantCache,
                :OrdinaryDiffEqMutableCache, :CompositeAlgorithm,
                :DerivativeOrderNotPossibleError,
                :_ode_interpolant, :_ode_interpolant!, :_ode_addsteps!,
                Symbol("@cache"),
                :calculate_residuals, :calculate_residuals!,
                Symbol("@truncate_stacktrace"),
            ),
        ),
        # Qualified accesses of solver-internal OrdinaryDiffEqCore names plus the Base
        # cartesian-loop macros. None are public; no public alternative exists.
        all_qualified_accesses_are_public = (;
            ignore = (
                :increment_nf!, :set_EEst!,
                :hermite_interpolant, :hermite_interpolant!,
                :interpolation_differential_vars,
                Symbol("@nexprs"), Symbol("@nif"),
            ),
        ),
    ),
)
