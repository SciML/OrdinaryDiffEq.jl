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

"""
    de_check_error(integrator::DEIntegrator)

The `DEIntegrator` implementation of [`SciMLBase.check_error`](@ref): inspect `integrator`
and return the `ReturnCode` describing whether integration may continue, reporting any
failure through `SciMLBase.report_integrator_failure`. Does not mutate the solution.
"""
function de_check_error(integrator::DEIntegrator)
    if integrator.sol.retcode ∉ (ReturnCode.Success, ReturnCode.Default)
        return integrator.sol.retcode
    end
    opts = integrator.opts
    # This implementation is intended to be used for ODEIntegrator and SDEIntegrator.

    if isnan(integrator.dt)
        SciMLBase.report_integrator_failure(integrator, Val(:dt_NaN))
        return ReturnCode.DtNaN
    end
    if integrator.iter > opts.maxiters
        SciMLBase.report_integrator_failure(integrator, Val(:max_iters))
        return ReturnCode.MaxIters
    end

    # Bail out if we take a step with dt less than the minimum value (which may be time
    # dependent), except when such a small timestep is successfully hitting a tstop exactly.
    # We also exit if the ODE is unstable according to a user chosen callback, but only if we
    # accepted the step, to avoid bailing out as unstable when we just took way too big a step.
    step_accepted = !hasproperty(integrator, :accept_step) || integrator.accept_step
    if !opts.force_dtmin && opts.adaptive
        if abs(integrator.dt) <= abs(opts.dtmin) &&
                (
                !step_accepted || (
                    hasproperty(opts, :tstops) ?
                        integrator.t + integrator.dt < integrator.tdir * first(opts.tstops) :
                        true
                )
            )
            SciMLBase.report_integrator_failure(integrator, Val(:dt_min_unstable))
            return ReturnCode.DtLessThanMin
        elseif !step_accepted && integrator.t isa AbstractFloat &&
                abs(integrator.dt) <= abs(eps(integrator.t))
            SciMLBase.report_integrator_failure(integrator, Val(:dt_epsilon))
            return ReturnCode.Unstable
        end
    end
    if step_accepted &&
            opts.unstable_check(integrator.dt, integrator.u, integrator.p, integrator.t)
        SciMLBase.report_integrator_failure(integrator, Val(:instability))
        return ReturnCode.Unstable
    end
    if last_step_failed(integrator)
        SciMLBase.report_integrator_failure(integrator, Val(:newton_convergence))
        return ReturnCode.ConvergenceFailure
    end
    return ReturnCode.Success
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
