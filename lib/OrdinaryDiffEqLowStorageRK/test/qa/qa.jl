using SciMLTesting, OrdinaryDiffEqLowStorageRK, Test

run_qa(
    OrdinaryDiffEqLowStorageRK;
    explicit_imports = true,
    ei_kwargs = (
        # All residual names below are non-public internals of their owning
        # package. They are imported/accessed from their true owner; the only
        # thing left is that those owners have not yet declared them `public`.
        # See SciML/OrdinaryDiffEq.jl#3776.
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore internal stats accessors / precompile probes
                :increment_nf!, :set_EEst!, :lorenz, :lorenz_oop,
                # SciMLBase specialization markers (used in @compile_workload)
                :FunctionWrapperSpecialize, :NoSpecialize,
                # Base.Broadcast internals used in ArrayFuse copyto!/materialize!
                :Broadcasted, :materialize!,
            )),
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore solver-interface internals
                Symbol("@cache"), :alg_cache, :perform_step!, :constvalue,
                :trivial_limiter!, :get_fsalfirstlast, :explicit_rk_docstring,
                :default_controller, :PIDController, :isfsal, :uses_uprev,
                :full_cache,
                :OrdinaryDiffEqAlgorithm, :OrdinaryDiffEqAdaptiveAlgorithm,
                :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
                # DiffEqBase residual-calculation internals
                :calculate_residuals, :calculate_residuals!,
            )),
    ),
)
