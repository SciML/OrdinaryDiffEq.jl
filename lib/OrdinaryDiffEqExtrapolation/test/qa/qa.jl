using SciMLTesting, OrdinaryDiffEqExtrapolation, Test
# Load Polyester so the extension exists and ExplicitImports analyzes it.
using Polyester

# `public` on a name another package owns counts as a public reexport to
# SciMLTesting, so the threading options need approving here too.
const THREADING_PUBLIC = (:Sequential, :BaseThreads, :PolyesterThreads)

# SciMLBase names re-exported for ordinary ODE usage; everything else stays behind `SciMLBase.`.
const SCIMLBASE_REEXPORTS = (:ODEProblem, :solve)

run_qa(
    OrdinaryDiffEqExtrapolation;
    reexports_allow = (SCIMLBASE_REEXPORTS..., THREADING_PUBLIC...),
    ei_kwargs = (;
        all_explicit_imports_are_public = (;
            # Package-internal hook the Polyester extension implements; deliberately not public.
            ignore = (:_polyester_foreach,),
        ),
    ),
)
