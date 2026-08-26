using SciMLTesting, OrdinaryDiffEqPRK, SciMLBase, Test

# `public` on a name another package owns counts as a public reexport to
# SciMLTesting, so the threading options need approving here too.
const THREADING_PUBLIC = (:Sequential, :BaseThreads, :PolyesterThreads)

run_qa(
    OrdinaryDiffEqPRK;
    reexports_allow = vcat(intersect(names(SciMLBase), names(OrdinaryDiffEqPRK)), collect(THREADING_PUBLIC)),
    explicit_imports = true,
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            ignore = (
                # OrdinaryDiffEqCore: deliberately non-public codegen macro
                # (kept internal alongside @fold/@OnDemandTableauExtract/@swap!).
                Symbol("@threaded"),
            ),
        ),
    ),
)
