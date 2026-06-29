using SciMLTesting, DiffEqDevTools, Test

run_qa(
    DiffEqDevTools;
    aqua_kwargs = (;
        ambiguities = false, piracies = false, unbound_args = false, stale_deps = false,
        deps_compat = (; check_extras = false),
    ),
    explicit_imports = true,
    ei_kwargs = (;
        # SciMLBase-owned solver-interface predicates that are accessed via SciMLBase
        # (their owner) but not yet declared `public` there.
        all_qualified_accesses_are_public = (;
            ignore = (:allowedkeywords, :calculate_ensemble_errors),
        ),
        # Abstract problem/solution/algorithm types + `@def` are owned by SciMLBase but
        # not yet `public` there; `ConvergenceSetup`/`ODERKTableau` live only in DiffEqBase
        # (not re-exported by SciMLBase) and are likewise not `public`.
        all_explicit_imports_are_public = (;
            ignore = (
                :AbstractDDEAlgorithm, :AbstractODESolution, :AbstractRODEProblem,
                :AbstractSDDEProblem, :AbstractBVProblem, Symbol("@def"),
                :ConvergenceSetup, :ODERKTableau,
            ),
        ),
    ),
)
