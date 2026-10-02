using SciMLTesting, StochasticDiffEq, Test
using JET
using DiffEqBase, DiffEqNoiseProcess, OrdinaryDiffEqNonlinearSolve, SciMLBase, CommonSolve
using StochasticDiffEqCore, StochasticDiffEqLowOrder, StochasticDiffEqRODE,
    StochasticDiffEqHighOrder, StochasticDiffEqMilstein, StochasticDiffEqROCK,
    StochasticDiffEqImplicit, StochasticDiffEqWeak, StochasticDiffEqIIF,
    StochasticDiffEqLeaping

# OrdinaryDiffEq umbrella pattern: allow public names that are also public on a
# reexported dependency.
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

# Umbrella `@reexport using X` modules (+ CommonSolve via DiffEqBase).
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
    # Narrow piracy allowlist for the Nothing-algorithm default-solver hook
    # (StochasticDiffEqCore uses treat_as_own similarly). Missing test-extras
    # [compat] entries are a follow-up (DiffEqDevTools uses check_extras=false).
    aqua_kwargs = (;
        piracies = (;
            treat_as_own = (SciMLBase.supports_solve_rng, SciMLBase.AbstractSDEProblem),
        ),
        deps_compat = (; check_extras = false),
    ),
    jet_kwargs = (; target_modules = (StochasticDiffEq,), mode = :typo),
    explicit_imports = true,
    ei_kwargs = (;
        no_implicit_imports = (;
            skip = REEXPORT_SKIP,
            ignore = (:Reexport, Symbol("@reexport")),
        ),
        # GlobalDiffEq: SciMLBase.__init/__solve are not declared public.
        all_qualified_accesses_are_public = (;
            ignore = (:__init, :__solve),
        ),
    ),
)
