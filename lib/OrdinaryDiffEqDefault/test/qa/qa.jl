using SciMLTesting, OrdinaryDiffEqDefault, SciMLBase, Test

# `DefaultSolverChoice` is an `EnumX.@enumx`-generated submodule that
# ExplicitImports cannot statically analyze. Its members (`Tsit5`, `Vern7`,
# `Rosenbrock23`, `Rodas5P`, `FBDF`, ...) share names with the solver types
# imported into `OrdinaryDiffEqDefault`, so the stale-import analysis
# mis-attributes those still-used imports as unused.
const ENUM_SUBMODULE = (OrdinaryDiffEqDefault.DefaultSolverChoice,)

run_qa(
    OrdinaryDiffEqDefault;
    # Approve the SciMLBase names this package re-exports. The list itself and the rule
    # behind it are checked repo-wide by test/qa/qa_tests.jl against docs/src/api/reexports.md.
    reexports_allow = intersect(names(SciMLBase), names(OrdinaryDiffEqDefault)),
    aqua_kwargs = (; piracies = false),
    explicit_imports = true,
    ei_kwargs = (;
        no_implicit_imports = (; allow_unanalyzable = ENUM_SUBMODULE),
        no_stale_explicit_imports = (;
            allow_unanalyzable = ENUM_SUBMODULE,
            # Used by `default_alg.jl` (solver constructors) but mis-flagged
            # stale because of the unanalyzable enum submodule above.
            ignore = (:Tsit5, :Vern7, :Rosenbrock23, :Rodas5P, :FBDF),
        ),
    ),
)
