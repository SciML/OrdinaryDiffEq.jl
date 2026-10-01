using SciMLTesting, StochasticDiffEq, Test
using JET
using DiffEqBase, DiffEqNoiseProcess, OrdinaryDiffEqNonlinearSolve, SciMLBase, CommonSolve
using StochasticDiffEqCore, StochasticDiffEqLowOrder, StochasticDiffEqRODE,
    StochasticDiffEqHighOrder, StochasticDiffEqMilstein, StochasticDiffEqROCK,
    StochasticDiffEqImplicit, StochasticDiffEqWeak, StochasticDiffEqIIF,
    StochasticDiffEqLeaping

# Umbrella package: `@reexport using` of DiffEqBase, the StochasticDiffEq*
# solver sublibraries, and DiffEqNoiseProcess is the public surface. Mirror the
# OrdinaryDiffEq umbrella QA pattern (`test/qa/qa_tests.jl`): allow only names
# that are both public on this package and public on a reexported dependency.
const STOCHASTICDIFFEQ_REEXPORTS = intersect(
    public_api_names(StochasticDiffEq),
    union(
        public_api_names(DiffEqBase),
        public_api_names(DiffEqNoiseProcess),
        public_api_names(OrdinaryDiffEqNonlinearSolve),
        public_api_names(StochasticDiffEqCore),
        public_api_names(StochasticDiffEqLowOrder),
        public_api_names(StochasticDiffEqRODE),
        public_api_names(StochasticDiffEqHighOrder),
        public_api_names(StochasticDiffEqMilstein),
        public_api_names(StochasticDiffEqROCK),
        public_api_names(StochasticDiffEqImplicit),
        public_api_names(StochasticDiffEqWeak),
        public_api_names(StochasticDiffEqIIF),
        public_api_names(StochasticDiffEqLeaping),
        (
            :DiffEqBase, :DiffEqNoiseProcess, :OrdinaryDiffEqNonlinearSolve,
            :StochasticDiffEqCore, :StochasticDiffEqLowOrder, :StochasticDiffEqRODE,
            :StochasticDiffEqHighOrder, :StochasticDiffEqMilstein,
            :StochasticDiffEqROCK, :StochasticDiffEqImplicit,
            :StochasticDiffEqWeak, :StochasticDiffEqIIF,
            :StochasticDiffEqLeaping,
        ),
    ),
)

# Modules whose exports are intentionally brought in by `@reexport using X`
# (umbrella surface). Same `skip` idiom as OrdinaryDiffEq's umbrella QA for
# SciMLBase; CommonSolve owns `solve`/`init`/`solve!`/`step!` reached via
# DiffEqBase's reexport.
const REEXPORT_SKIP = (
    Base, Core,
    DiffEqBase, DiffEqNoiseProcess, OrdinaryDiffEqNonlinearSolve,
    StochasticDiffEqCore, StochasticDiffEqLowOrder, StochasticDiffEqRODE,
    StochasticDiffEqHighOrder, StochasticDiffEqMilstein, StochasticDiffEqROCK,
    StochasticDiffEqImplicit, StochasticDiffEqWeak, StochasticDiffEqIIF,
    StochasticDiffEqLeaping,
    SciMLBase, CommonSolve,
)

run_qa(
    StochasticDiffEq;
    reexports_allow = STOCHASTICDIFFEQ_REEXPORTS,
    # `supports_solve_rng(::AbstractSDEProblem, ::Nothing)` matches
    # OrdinaryDiffEqDefault's Nothing-algorithm default-solver hook; that sibling
    # disables Aqua piracy for the same pattern.
    # `check_extras = false` matches DiffEqDevTools: many test extras lack
    # [compat] entries and fixing them is out of scope for wiring QA.
    aqua_kwargs = (;
        piracies = false,
        deps_compat = (; check_extras = false),
    ),
    # Scope JET to this package in `:typo` mode, matching StochasticDiffEqLeaping /
    # StochasticDiffEqImplicit (and SciMLTesting's default for solver packages).
    jet_kwargs = (; target_modules = (StochasticDiffEq,), mode = :typo),
    explicit_imports = true,
    ei_kwargs = (;
        no_implicit_imports = (;
            skip = REEXPORT_SKIP,
            # Reexport machinery itself (OrdinaryDiffEqNordsieck / HighOrderRK).
            ignore = (:Reexport, Symbol("@reexport")),
        ),
        # `SciMLBase.__init` / `__solve` are the owner-path solve entry points
        # (same ignore as GlobalDiffEq for `__solve`); they are not declared
        # `public` on SciMLBase.
        all_qualified_accesses_are_public = (;
            ignore = (:__init, :__solve),
        ),
    ),
)
