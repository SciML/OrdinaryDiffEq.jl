using SciMLTesting, OrdinaryDiffEqExtrapolation, Test

# `public` on a name another package owns counts as a public reexport to
# SciMLTesting, so the threading options need approving here too.
const THREADING_PUBLIC = (:Sequential, :BaseThreads, :PolyesterThreads)

run_qa(
    OrdinaryDiffEqExtrapolation;
    reexports_allow = THREADING_PUBLIC,
)
