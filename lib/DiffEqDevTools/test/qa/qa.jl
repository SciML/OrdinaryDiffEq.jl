using SciMLTesting, DiffEqDevTools, Test

run_qa(
    DiffEqDevTools;
    aqua_kwargs = (;
        ambiguities = false, piracies = false, unbound_args = false, stale_deps = false,
        deps_compat = (; check_extras = false),
    ),
    explicit_imports = true,
    ei_broken = (:no_implicit_imports, :no_stale_explicit_imports, :all_explicit_imports_via_owners, :all_qualified_accesses_via_owners, :all_qualified_accesses_are_public, :all_explicit_imports_are_public),  # known-broken; see SciML/OrdinaryDiffEq.jl#3776
)
