using SciMLTesting, OrdinaryDiffEqPRK, Test

run_qa(
    OrdinaryDiffEqPRK;
    explicit_imports = true,
    ei_kwargs = (;
        # OrdinaryDiffEqCore solver-interface names imported from their owner
        # module (OrdinaryDiffEqCore) but not yet declared `public`. These are
        # the OrdinaryDiffEq internal extension/dispatch API; no public
        # alternative exists on the registered releases. Tracked for
        # make-public in SciML/OrdinaryDiffEq.jl#3776.
        all_explicit_imports_are_public = (;
            ignore = (
                Symbol("@cache"),
                Symbol("@threaded"),
                :alg_cache,
                :constvalue,
                :generic_solver_docstring,
                :get_fsalfirstlast,
                :isthreaded,
                :OrdinaryDiffEqAlgorithm,
                :OrdinaryDiffEqConstantCache,
                :OrdinaryDiffEqMutableCache,
                :perform_step!,
                :unwrap_alg,
            ),
        ),
    ),
)
