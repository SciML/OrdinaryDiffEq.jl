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
        # OrdinaryDiffEqCore internal API the default-alg machinery depends on
        # (no public replacement; see SciML/OrdinaryDiffEq.jl#3776).
        all_explicit_imports_are_public = (;
            ignore = (
                :AutoAlgSwitch, :AutoSwitchCache, :CompositeAlgorithm,
                :alg_stability_size, :beta1_default, :beta2_default,
                :default_autoswitch, :is_mass_matrix_alg, :isdefaultalg,
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                # SciMLBase specialization sentinels used only inside opt-in
                # (default-off) precompile-workload branches.
                :FunctionWrapperSpecialize, :NoSpecialize,
                # OrdinaryDiffEqCore precompile test problems.
                :lorenz, :lorenz_oop,
            ),
        ),
    ),
)
