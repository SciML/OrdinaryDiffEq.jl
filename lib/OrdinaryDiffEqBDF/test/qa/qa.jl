using SciMLTesting, OrdinaryDiffEqBDF, Test

# The names below are internal (non-public) API of sibling solver sublibraries
# (OrdinaryDiffEqCore / OrdinaryDiffEqSDIRK / OrdinaryDiffEqNonlinearSolve) plus a
# handful of SciMLBase / DiffEqBase / TruncatedStacktraces internals that BDF must
# use to plug into the integrator framework. They are owned by those packages and
# have no public replacement, so they are imported/accessed from their owner and
# ignored here. Each becomes a candidate for a `public` declaration upstream; see
# SciML/OrdinaryDiffEq.jl#3776.
const QUALIFIED_NONPUBLIC = (
    :CompositeControllerCache, :accept_step_controller, :get_EEst, :increment_nf!,
    :post_newton_controller!, :set_EEst!, :step_accept_controller!,
    :step_reject_controller!, :stepsize_controller!,
    :lorenz, :lorenz_oop,                       # precompile-workload test problems
    :FunctionWrapperSpecialize, :NoSpecialize,  # SciMLBase specialization markers
)

const EXPLICIT_IMPORT_NONPUBLIC = (
    Symbol("@cache"), Symbol("@truncate_stacktrace"),
    :AbstractController, :AbstractControllerCache, :CommonControllerOptions,
    :DAEAlgorithm, :DerivativeOrderNotPossibleError, :IController,
    :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
    :OrdinaryDiffEqNewtonAdaptiveAlgorithm, :OrdinaryDiffEqNewtonAlgorithm,
    :COEFFICIENT_MULTISTEP, :DIRK,
    :ESDIRKIMEXCache, :ESDIRKIMEXConstantCache, :ImplicitEulerESDIRKIMEXTableau,
    :NLNewton, :build_nlsolver, :du_alias_or_new, :markfirststage!,
    :nlsolve!, :nlsolvefail, :isnewton, :set_new_W!,
    :alg_cache, :alg_extrapolates, :constvalue, :default_controller, :error_constant,
    :_fixup_ad, :gamma_default, :generic_solver_docstring,
    :get_current_adaptive_order, :get_current_alg_order, :get_current_qmax,
    :get_EEst, :get_failfactor, :get_fsalfirstlast, :get_gamma, :get_qmax,
    :get_qsteady_max, :get_qsteady_min, :has_special_newton_error,
    :has_stiff_interpolation, :isfsal, :issplit,
    :_ode_addsteps!, :_ode_interpolant, :_ode_interpolant!,
    :perform_step!, :post_newton_controller!, :qmax_default,
    :qsteady_max_default, :qsteady_min_default, :resolve_basic, :_resolved_QT,
    :set_discontinuity, :setup_controller_cache,
    :step_accept_controller!, :step_reject_controller!, :stepsize_controller!,
    :trivial_limiter!, :unwrap_alg, :_unwrap_val,
    :calculate_residuals, :calculate_residuals!,  # DiffEqBase-owned, non-public
)

run_qa(
    OrdinaryDiffEqBDF;
    explicit_imports = true,
    ei_kwargs = (
        all_qualified_accesses_are_public = (; ignore = QUALIFIED_NONPUBLIC),
        all_explicit_imports_are_public = (; ignore = EXPLICIT_IMPORT_NONPUBLIC),
    ),
)
