using SciMLTesting, OrdinaryDiffEqDefault, Test

# `DefaultSolverChoice` is an `EnumX.@enumx`-generated submodule that
# ExplicitImports cannot statically analyze. Its members (`Tsit5`, `Vern7`,
# `Rosenbrock23`, `Rodas5P`, `FBDF`, ...) share names with the solver types
# imported into `OrdinaryDiffEqDefault`, so the stale-import analysis
# mis-attributes those still-used imports as unused.
const ENUM_SUBMODULE = (OrdinaryDiffEqDefault.DefaultSolverChoice,)

run_qa(
    OrdinaryDiffEqDefault;
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
        all_qualified_accesses_are_public = (;
            # `lorenz`/`lorenz_oop` are OrdinaryDiffEqCore precompile-workload
            # test problems (defined in `precompilation_setup.jl`), deliberately
            # not part of its public extension surface.
            ignore = (:lorenz, :lorenz_oop),
        ),
    ),
)
