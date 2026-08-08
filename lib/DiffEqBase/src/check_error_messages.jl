# `SciMLBase.check_error` only detects failures; it reports each one through the
# `report_integrator_failure` hook. The wording and the `DEVerbosity` toggles that
# gate it are differential-equation specific, so they live here.

"""
    SciMLBase.report_integrator_failure(integrator, ::Val{reason})

Emit the diagnostic for a failure mode detected by `SciMLBase.check_error`, gated by
the `DEVerbosity` toggle named `reason`.

`reason` is a static parameter rather than a runtime `Symbol` so the toggle lookup
constant-folds and `check_error` stays allocation-free when nothing is emitted.
"""
@inline function SciMLBase.report_integrator_failure(integ::DEIntegrator, ::Val{reason}) where {reason}
    @SciMLMessage(
        lazy"$(failure_message(integ, Val(reason)))$(instability_diagnostic(integ))",
        integ.opts.verbose, reason
    )
    return nothing
end

# Trailing symbolic/numeric detail appended to every failure message. Getting here already
# means the toggle for this failure is active, so only the symbolic half needs its own
# gate; it additionally needs an MTK system. When it says something the numeric half drops
# its Jacobian report rather than repeating it.
function instability_diagnostic(integ)
    verbose = integ.opts.verbose
    symbolic = if SciMLBase.has_mtk_sys(integ) && verbosity_to_bool(verbose.symbolic_diagnostic)
        SciMLBase.diagnose_symbolic_instability(integ.f.sys, integ.u, integ.uprev)
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
