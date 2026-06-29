using SciMLTesting, OrdinaryDiffEqAdamsBashforthMoulton, Test

run_qa(
    OrdinaryDiffEqAdamsBashforthMoulton;
    explicit_imports = true,
    ei_kwargs = (
        # All residual names below are non-public internals of their owning
        # package (OrdinaryDiffEqCore / DiffEqBase / OrdinaryDiffEqLowOrderRK).
        # They are imported from their true owner; the only thing left is that
        # those owners have not yet declared them `public`. See SciML/OrdinaryDiffEq.jl#3776.
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internal stats/error-estimate accessors
                :get_EEst, :set_EEst!, :increment_nf!,
            )),
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore solver-interface internals
                Symbol("@cache"), :alg_cache, :perform_step!, :constvalue,
                :trivial_limiter!, :get_fsalfirstlast, :generic_solver_docstring,
                :default_controller, :IController,
                :OrdinaryDiffEqAlgorithm, :OrdinaryDiffEqAdaptiveAlgorithm,
                :OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm,
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
                # DiffEqBase residual-calculation internals
                :calculate_residuals, :calculate_residuals!,
                # OrdinaryDiffEqLowOrderRK starter-cache types
                :BS3Cache, :BS3ConstantCache, :RK4Cache, :RK4ConstantCache,
            )),
    ),
)
