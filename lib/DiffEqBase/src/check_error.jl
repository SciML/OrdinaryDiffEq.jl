# `check_error`/`check_error!` themselves stay in SciMLBase so that a DiffEqBase older than
# this one still has working implementations. Only the wording of the diagnostics and the
# `DEVerbosity` toggles that gate them are differential-equation specific, so those live
# here, as methods of the `report_integrator_failure` hook SciMLBase calls.

#cap diagnostic output to avoid OOM errors when printing large symbolic systems or user data
const DIAGNOSTIC_OBJECT_CHARS = 160
const DIAGNOSTIC_REPORT_CHARS = 4000

@noinline function truncate_str(x, limit::Int = DIAGNOSTIC_OBJECT_CHARS)::String
    buf = IOBuffer(maxsize = 4limit)
    print(IOContext(buf, :limit => true, :displaysize => (10, limit)), x)
    s = String(take!(buf))
    length(s) <= limit && return s
    return first(s, limit) * "… (truncated)"
end

if isdefined(SciMLBase, :report_integrator_failure)
    """
        SciMLBase.report_integrator_failure(integrator::DEIntegrator, ::Val{reason})

    Emit the diagnostic for a failure mode detected by `SciMLBase.check_error`, gated by
    the `DEVerbosity` toggle named `reason`. `reason` is one of `:dt_NaN`, `:max_iters`,
    `:dt_min_unstable`, `:dt_epsilon`, `:instability` or `:newton_convergence`.

    Dispatch is on the integrator type as well, so a solver package can replace the wording
    for one failure and inherit the rest. Implementations must not affect control flow and
    their return value is ignored.

    `reason` is a static parameter rather than a runtime `Symbol` so the toggle lookup
    constant-folds and `check_error` stays allocation-free when nothing is emitted.
    """
    @inline function SciMLBase.report_integrator_failure(
            integ::DEIntegrator, ::Val{reason}
        ) where {reason}
        # emit_message never calls the closure for silent toggle
        @SciMLMessage(
            () -> "$(failure_message(integ, Val(reason)))$(instability_diagnostic(integ))",
            integ.opts.verbose, reason
        )
        return nothing
    end
end

# Trailing symbolic/numeric detail appended to every failure message. Getting here already
# means the toggle for this failure is active, so only the symbolic half needs its own
# gate; it additionally needs an MTK system. When it says something the numeric half drops
# its Jacobian report rather than repeating it.
function instability_diagnostic(integ)
    verbose = integ.opts.verbose
    symbolic = if SciMLBase.has_mtk_sys(integ) && verbosity_to_bool(verbose.symbolic_diagnostic)
        # the symbolic report embeds the model's equations, so it is the one unbounded half
        truncate_str(
            SciMLBase.diagnose_symbolic_instability(integ.f.sys, integ.u, integ.uprev),
            DIAGNOSTIC_REPORT_CHARS
        )
    else
        ""
    end
    numeric = SciMLBase.log_numerical_instability(integ; jacobian_logging = symbolic == "")
    return numeric * symbolic
end

eest_suffix(integ) = isdefined(integ, :EEst) ?
    lazy", step error estimate = $(integ.EEst)" : ""

failure_message(integ, ::Val{:dt_NaN}) =
    "NaN dt detected. Likely a NaN value in the state, parameters, or derivative value caused this outcome."

failure_message(integ, ::Val{:max_iters}) =
    "Interrupted. Larger maxiters is needed. If you are using an integrator for non-stiff ODEs or an automatic switching algorithm (the default), you may want to consider using a method for stiff equations. See the solver pages for more details (e.g. https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Stiff-Problems)."

failure_message(integ, ::Val{:dt_min_unstable}) =
    lazy"dt($(integ.dt)) <= dtmin($(integ.opts.dtmin)) at t=$(integ.t)$(eest_suffix(integ)). Aborting. There is either an error in your model specification or the true solution is unstable."

failure_message(integ, ::Val{:dt_epsilon}) =
    lazy"At t=$(integ.t), dt was forced below floating point epsilon $(integ.dt)$(eest_suffix(integ)). Aborting. There is either an error in your model specification or the true solution is unstable (or it cannot be represented in $(eltype(integ.u)) precision)."

failure_message(integ, ::Val{:instability}) = "Instability detected. Aborting."

failure_message(integ, ::Val{:newton_convergence}) =
    "Newton steps could not converge and algorithm is not adaptive. Use a lower dt."

# `a && b()` and `a || b()` for the conditions below. A non-`Bool` operand, such as a boolean traced
# by a compiler, cannot short-circuit, so both sides are evaluated.
@inline check_error_and(b, a::Bool) = a && b()
@inline check_error_and(b, a) = a & b()
@inline check_error_or(b, a::Bool) = a || b()
@inline check_error_or(b, a) = a | b()

check_error_step_accepted(integrator) =
    !hasproperty(integrator, :accept_step) || integrator.accept_step

"""
    check_error_failed_retcode(integrator) -> Bool

Whether `integrator.sol.retcode` already records a failure, i.e. is neither
`ReturnCode.Success` nor `ReturnCode.Default`. `de_check_error` returns such a code
unchanged before checking anything else.
"""
check_error_failed_retcode(integrator) =
    integrator.sol.retcode ∉ (ReturnCode.Success, ReturnCode.Default)

"""
    check_error_dt_below_time_eps(integrator) -> Bool

Whether `abs(integrator.dt) <= abs(eps(integrator.t))`, the comparison of the `dt_epsilon`
failure check of `de_check_error`; `false` when the time is not an `AbstractFloat`.
Integrators whose time can be a non-`AbstractFloat` wrapper of a float, such as a traced
number, extend this with the exact comparison.
"""
check_error_dt_below_time_eps(integrator) =
    integrator.t isa AbstractFloat && abs(integrator.dt) <= abs(eps(integrator.t))

check_error_dt_nan(integrator, step_accepted) = isnan(integrator.dt)

check_error_maxiters(integrator, step_accepted) = integrator.iter > integrator.opts.maxiters

# Bail out if we take a step with dt less than the minimum value (which may be time
# dependent), except when such a small timestep is successfully hitting a tstop exactly.
function check_error_dtmin(integrator, step_accepted)
    opts = integrator.opts
    (!opts.force_dtmin && opts.adaptive) || return false
    return check_error_and(abs(integrator.dt) <= abs(opts.dtmin)) do
        check_error_or(!step_accepted) do
            hasproperty(opts, :tstops) ?
                integrator.t + integrator.dt < integrator.tdir * first(opts.tstops) : true
        end
    end
end

function check_error_dt_epsilon(integrator, step_accepted)
    opts = integrator.opts
    (!opts.force_dtmin && opts.adaptive) || return false
    return check_error_and(() -> check_error_dt_below_time_eps(integrator), !step_accepted)
end

# Only judge accepted steps as unstable, to avoid bailing out as unstable when we just took
# way too big a step.
function check_error_instability(integrator, step_accepted)
    return check_error_and(step_accepted) do
        integrator.opts.unstable_check(integrator.dt, integrator.u, integrator.p, integrator.t)
    end
end

check_error_newton_convergence(integrator, step_accepted) = last_step_failed(integrator)

# The failure conditions after an existing failure code, in priority order: the first that
# holds determines the code and the diagnostic.
const CHECK_ERROR_CONDITIONS = (
    (check_error_dt_nan, ReturnCode.DtNaN, Val(:dt_NaN)),
    (check_error_maxiters, ReturnCode.MaxIters, Val(:max_iters)),
    (check_error_dtmin, ReturnCode.DtLessThanMin, Val(:dt_min_unstable)),
    (check_error_dt_epsilon, ReturnCode.Unstable, Val(:dt_epsilon)),
    (check_error_instability, ReturnCode.Unstable, Val(:instability)),
    (check_error_newton_convergence, ReturnCode.ConvergenceFailure, Val(:newton_convergence)),
)

@inline report_first_check_error_failure(integrator, step_accepted, ::Tuple{}) =
    ReturnCode.Success
@inline function report_first_check_error_failure(integrator, step_accepted, conditions)
    condition, code, reason = first(conditions)
    if condition(integrator, step_accepted)
        SciMLBase.report_integrator_failure(integrator, reason)
        return code
    end
    return report_first_check_error_failure(integrator, step_accepted, Base.tail(conditions))
end

"""
    de_check_error(integrator::DEIntegrator)

The `DEIntegrator` implementation of [`SciMLBase.check_error`](@ref): inspect `integrator`
and return the `ReturnCode` describing whether integration may continue, reporting any
failure through `SciMLBase.report_integrator_failure`. Does not mutate the solution.
Intended for ODE and SDE integrators. `staged_check_error` evaluates the same
conditions without branching.
"""
function de_check_error(integrator::DEIntegrator)
    check_error_failed_retcode(integrator) && return integrator.sol.retcode
    return report_first_check_error_failure(
        integrator, check_error_step_accepted(integrator), CHECK_ERROR_CONDITIONS
    )
end

@inline staged_first_check_error_failure(integrator, step_accepted, ::Tuple{}) =
    (false, ReturnCode.Success)
@inline function staged_first_check_error_failure(integrator, step_accepted, conditions)
    condition, code, _ = first(conditions)
    holds = condition(integrator, step_accepted)
    later_failed, later_code = staged_first_check_error_failure(
        integrator, step_accepted, Base.tail(conditions)
    )
    return holds | later_failed, ifelse(holds, code, later_code)
end

"""
    staged_check_error(integrator::DEIntegrator) -> (failed, code)

Evaluate every condition of `de_check_error` and select the code of the first that
holds with `ifelse` instead of returning early, so that a tracing compiler can stage the
check when the integrator state is traced. `code` equals `de_check_error(integrator)` and
`failed` is `code != ReturnCode.Success`; no diagnostic is reported. Conditions are
evaluated even when an earlier one holds, so the integrator's `unstable_check` must be free
of side effects.
"""
function staged_check_error(integrator::DEIntegrator)
    failed, code = staged_first_check_error_failure(
        integrator, check_error_step_accepted(integrator), CHECK_ERROR_CONDITIONS
    )
    stored = check_error_failed_retcode(integrator)
    return stored | failed, ifelse(stored, integrator.sol.retcode, code)
end

"""
    de_check_error!(integrator::DEIntegrator)

Run `SciMLBase.check_error`, store the code in `integrator.sol.retcode`, and return it,
calling `postamble!` when the code is not `ReturnCode.Success`. Dispatching through
`SciMLBase.check_error` rather than [`de_check_error`](@ref) keeps solver-specific
`check_error` overrides in effect.
"""
function de_check_error!(integrator::DEIntegrator)
    code = SciMLBase.check_error(integrator)
    integrator.sol = solution_new_retcode(integrator.sol, code)
    if code != ReturnCode.Success
        postamble!(integrator)
    end
    return code
end

# Wire the implementations in only if SciMLBase stopped shipping its own, so DiffEqBase can
# stand alone without overwriting SciMLBase's methods while they are still there. Defining
# them unconditionally is a hard `Method overwriting is not permitted during Module
# precompilation` error, not a warning, so these guards cannot simply be dropped -- they come
# out only once SciMLBase stops shipping its `DEIntegrator` methods.
#below will be removed once scimlbase stops shipping, but we keep it for now
if !hasmethod(SciMLBase.check_error, Tuple{DEIntegrator})
    SciMLBase.check_error(integrator::DEIntegrator) = de_check_error(integrator)
end
if !hasmethod(SciMLBase.check_error!, Tuple{DEIntegrator})
    SciMLBase.check_error!(integrator::DEIntegrator) = de_check_error!(integrator)
end
