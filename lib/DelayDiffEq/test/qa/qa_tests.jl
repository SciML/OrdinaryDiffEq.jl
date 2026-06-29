using SciMLTesting, DelayDiffEq, Test
import SciMLBase
using OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, StochasticDiffEqAlgorithm,
    StochasticDiffEqRODEAlgorithm,
    StochasticDiffEqConstantCache, StochasticDiffEqMutableCache

# Names that are non-public in their owning package on the registered releases
# resolved here (SciMLBase 3.30, OrdinaryDiffEqCore 4.4, DiffEqBase 7, etc.).
# DelayDiffEq builds a DDE integrator on top of the OrdinaryDiffEq solver
# internals, so it legitimately reaches into these internal stepping/cache/
# discontinuity hooks. Each is a candidate for `public` declaration upstream
# (tracked in SciML/OrdinaryDiffEq.jl#3776); until then they are listed
# explicitly rather than hidden behind a blanket `ei_broken`.
const QUALIFIED_INTERNAL = (
    # SciMLBase internals (problem/function/solution interface, init plumbing)
    :AbstractDDEFunction, :AbstractHistoryFunction, :AbstractSDDEFunction,
    :AbstractSDEFunction, :DAEInitializationAlgorithm,
    :get_save_idxs_and_saved_subsystem, :has_initializeprob, :last_step_failed,
    :postamble!, :save_discretes_if_enabled!, :unwrap_cache,
    # DiffEqBase internals
    :calculate_residuals, :calculate_residuals!, :check_prob_alg_pairing,
    :ODE_DEFAULT_UNSTABLE_CHECK, :prepare_alg, :prob2dtmin, :Stats,
    # OrdinaryDiffEqCore internals (integrator stepping/cache/controller/discontinuity)
    :AbstractNLSolver, :alg_adaptive_order, :alg_autodiff, :alg_cache,
    :alg_difftype, :alg_extrapolates, :alg_maximum_order, :apply_step!,
    :AutoSwitchCache, :_change_t_via_interpolation!, :concrete_jac,
    :current_interpolant, :current_interpolant!, :DefaultCache,
    :default_controller, :DEOptions, :determine_controller_datatype,
    :first_discontinuity, :get_chunksize, :get_chunksize_int,
    :get_differential_vars, :get_EEst, :get_fsalfirstlast,
    :_get_next_step_tstop, :_get_tstop_target, :handle_callback_modifiers!,
    :handle_discontinuities!, :handle_dt!, :handle_starting_time_discontinuity!,
    :handle_tstop!, :has_discontinuity, :initialize_callbacks!,
    :_initialize_dae!, :initialize_d_discontinuities, :initialize_saveat,
    :initialize_tstops, :init_ith_default_cache, :InterpolationData,
    :isautoswitch, :isdtchangeable, :isimplicit, :ismultistep, :_loopfooter!,
    :loopfooter!, :loopheader!, :modify_dt_for_tstops!, :nlsolve_f, :NLStatus,
    :_ode_addsteps!, :ode_addsteps!, :ode_determine_initdt,
    :OrdinaryDiffEqAlgorithm, :OrdinaryDiffEqCompositeAlgorithm, :perform_step!,
    :pop_discontinuity!, :_postamble!, :reinit_controller!, :resize_J_W!,
    :resize_nlsolver!, :_savevalues!, :set_EEst!, :_set_tstop_flag!,
    :setup_controller_cache, :shift_past_discontinuity!, :standardtag,
    :update_uprev!, :uses_uprev,
    # OrdinaryDiffEqNonlinearSolve internals (fixed-point solver hooks)
    :anderson, :anderson!, :compute_step!, :initial_η, :nlsolve!, :nlsolvefail,
    # OrdinaryDiffEqDifferentiation internals
    :resize_grad_config!, :resize_jac_config!,
    # ForwardDiff internal type used in discontinuity tracking
    :Dual,
)

# Explicit `using A: name` of names that are non-public in `A`. Same upstream
# make-public story as above (subset of the qualified names plus the abstract
# algorithm/cache supertypes DelayDiffEq subtypes and piracy-dispatches on).
const EXPLICIT_INTERNAL = (
    # SciMLBase
    :AbstractDDEAlgorithm, :AbstractSDDEProblem, :ProblemSolverPairingError,
    # OrdinaryDiffEqCore
    :AbstractNLSolverCache, :alg_extrapolates, :alg_maximum_order, :AutoSwitch,
    :CompositeAlgorithm, :SlowConvergence, :StochasticDiffEqAlgorithm,
    :StochasticDiffEqRODEAlgorithm,
    # OrdinaryDiffEqNonlinearSolve
    :NLAnderson, :NLFunctional,
    # OrdinaryDiffEqRosenbrock
    :RosenbrockMutableCache,
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
