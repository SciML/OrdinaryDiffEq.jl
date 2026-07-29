using SciMLTesting, OrdinaryDiffEqStabilizedIRK, Test

run_qa(
    OrdinaryDiffEqStabilizedIRK;
    # No docs/ tree here; the umbrella manual renders this package's API.
    api_docs_kwargs = (; rendered = false),
    reexports_allow = union(public_api_names(SciMLBase), (:SciMLBase,)),
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # Owner-internal cross-sublibrary hooks with no public wrapper.
                # OrdinaryDiffEqDifferentiation:
                :dolinsolve, :update_W!,
                # OrdinaryDiffEqNonlinearSolve:
                :build_nlsolver, :du_alias_or_new, :markfirststage!, :nlsolve!,
                # SciMLBase internals, not public on registered releases.
                :_unwrap_val, :_reshape, :_vec,
            ),
        ),
    ),
)
