using SciMLTesting, DelayDiffEq, Test
import SciMLBase
using OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, StochasticDiffEqAlgorithm,
    StochasticDiffEqRODEAlgorithm,
    StochasticDiffEqConstantCache, StochasticDiffEqMutableCache

# Names that are non-public in their owning package on the registered releases
# resolved here. The OrdinaryDiffEqCore/Differentiation/NonlinearSolve/DiffEqBase
# solver-author extension API is now declared `public` on this branch, so those
# names are dropped from the ignore lists below. What remains is the genuine
# residual: DDE-specific integrator hooks that OrdinaryDiffEqCore has not made
# public, the SciMLBase problem/function/solution abstract types still pending a
# `public` declaration (SciML/SciMLBase.jl#1412), and true external internals.
const QUALIFIED_INTERNAL = (
    # SciMLBase abstract problem/function types, pending SciMLBase#1412.
    :AbstractDDEFunction, :AbstractHistoryFunction, :AbstractSDDEFunction,
    :AbstractSDEFunction, :DAEInitializationAlgorithm,
    # SciMLBase init/save/step plumbing not yet public (make-public candidates).
    :get_save_idxs_and_saved_subsystem, :has_initializeprob, :last_step_failed,
    :postamble!, :save_discretes_if_enabled!, :unwrap_cache,
    # DiffEqBase internal stats type.
    :Stats,
    # OrdinaryDiffEqCore DDE-specific integrator stepping / discontinuity /
    # tstop / saveat hooks — reached by the DDE integrator but not part of the
    # OrdinaryDiffEqCore public extension surface declared on this branch.
    :_change_t_via_interpolation!, :current_interpolant!,
    :determine_controller_datatype, :first_discontinuity, :get_chunksize_int,
    :_get_next_step_tstop, :_get_tstop_target, :handle_discontinuities!,
    :handle_dt!, :handle_starting_time_discontinuity!, :handle_tstop!,
    :has_discontinuity, :initialize_callbacks!, :initialize_d_discontinuities,
    :initialize_saveat, :initialize_tstops, :init_ith_default_cache,
    :_loopfooter!, :loopfooter!, :loopheader!, :modify_dt_for_tstops!,
    :pop_discontinuity!, :_postamble!, :_savevalues!, :_set_tstop_flag!,
    :shift_past_discontinuity!, :update_uprev!,
    # ForwardDiff internal type used in discontinuity tracking.
    :Dual,
)

# Explicit `using A: name` of names non-public in `A`.
const EXPLICIT_INTERNAL = (
    # SciMLBase abstract algorithm/problem types + pairing error, pending
    # SciMLBase#1412.
    :AbstractDDEAlgorithm, :AbstractSDDEProblem, :ProblemSolverPairingError,
)

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
        no_implicit_imports = (; skip = (Base, Core)),
        all_qualified_accesses_are_public = (; ignore = QUALIFIED_INTERNAL),
        all_explicit_imports_are_public = (; ignore = EXPLICIT_INTERNAL),
    ),
)
