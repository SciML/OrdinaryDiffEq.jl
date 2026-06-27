using SciMLTesting, DelayDiffEq, Test
import SciMLBase
using OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, StochasticDiffEqAlgorithm,
    StochasticDiffEqRODEAlgorithm,
    StochasticDiffEqConstantCache, StochasticDiffEqMutableCache

run_qa(
    DelayDiffEq;
    aqua_kwargs = (;
        ambiguities = (; recursive = false),
        # piracy is allowed for the default solver methods and SDE integration
        piracies = (;
            treat_as_own = [
                SciMLBase.DDEProblem,
                OrdinaryDiffEqAlgorithm,
                StochasticDiffEqAlgorithm,
                StochasticDiffEqRODEAlgorithm,
                StochasticDiffEqConstantCache,
                StochasticDiffEqMutableCache,
            ],
        ),
    ),
    explicit_imports = true,
    ei_kwargs = (;
        no_implicit_imports = (; skip = (Base, Core), ignore = (Symbol("@reexport"),)),
        no_stale_explicit_imports = (; ignore = (:AbstractVerbositySpecifier, :Standard)),
    ),
    # known-broken; see SciML/OrdinaryDiffEq.jl#3776
    ei_broken = (
        :all_explicit_imports_via_owners,
        :all_qualified_accesses_are_public,
        :all_explicit_imports_are_public,
    ),
)
