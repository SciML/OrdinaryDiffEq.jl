using SciMLTesting, StochasticDiffEqLeaping, Test
using JET

run_qa(
    StochasticDiffEqLeaping;
    reexports_allow = union(public_api_names(StochasticDiffEqCore), (:StochasticDiffEqCore,)),
    # Scope JET to this package in `:typo` mode (the OrdinaryDiffEq* solver-sublibrary
    # convention); `target_defined_modules = true` is deprecated in JET.
    jet_kwargs = (; target_modules = (StochasticDiffEqLeaping,), mode = :typo),
    explicit_imports = true,
    ei_kwargs = (;
        # `@..` is the SciML fused-broadcast macro, owned by FastBroadcast and
        # re-exported through DiffEqBase (FastBroadcast is not a direct dep here).
        all_explicit_imports_via_owners = (; ignore = (Symbol("@.."),)),
        # `pois_rand` is owned by PoissonRandom and re-exported by JumpProcesses,
        # the natural dependency that owns the jump-process interface here.
        all_qualified_accesses_via_owners = (; ignore = (:pois_rand,)),
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqNonlinearSolve — owner-internal, no public alternative
                :build_nlsolver, :markfirststage!, :nlsolve!, :nlsolvefail,
                # FastBroadcast fused-broadcast macro, re-exported via DiffEqBase.
                Symbol("@.."),
                # StochasticDiffEqCore internal cache-alloc macro (non-public).
                Symbol("@cache"),
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # OrdinaryDiffEqNonlinearSolve — owner-internal, no public alternative
                :build_nlsolver, :markfirststage!, :nlsolve!, :nlsolvefail,
                # OrdinaryDiffEqCore internal a-posteriori error-estimate predicate
                # (non-public).
                :isaposteriori,
                # PoissonRandom sampler, re-exported via JumpProcesses (non-public).
                :pois_rand,
            ),
        ),
    ),
)
