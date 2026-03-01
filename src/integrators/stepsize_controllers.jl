function stepsize_controller!(integrator::SDEIntegrator, controller::PIController, alg)
    integrator.q11 = DiffEqBase.value(FastPower.fastpower(integrator.EEst, controller.beta1))
    integrator.q = DiffEqBase.value(integrator.q11 / FastPower.fastpower(integrator.qold, controller.beta2))
    return @fastmath integrator.q = DiffEqBase.value(
        max(
            inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), integrator.q / integrator.opts.gamma)
        )
    )
end

# Bridge: ODE's _loopfooter! calls step_accept_controller!(integrator, alg, q).
# SDE uses 2-arg dispatch (ignores q, uses integrator.q instead).
@inline function step_accept_controller!(integrator::SDEIntegrator, alg, q)
    return step_accept_controller!(integrator, alg)
end

@inline function step_accept_controller!(integrator::SDEIntegrator, alg)
    return step_accept_controller!(integrator, integrator.opts.controller, alg)
end

function step_accept_controller!(integrator::SDEIntegrator, controller::PIController, alg)
    integrator.qold = max(integrator.EEst, integrator.opts.qoldinit)
    return DiffEqBase.value(integrator.dt / integrator.q) * oneunit(integrator.dt)
end

function step_reject_controller!(integrator::SDEIntegrator, controller::PIController, alg)
    return integrator.dt = integrator.dt / min(inv(integrator.opts.qmin), integrator.q11 / integrator.opts.gamma)
end

function stepsize_controller!(integrator::SDEIntegrator, alg::TauLeaping)
    return nothing
end

function step_accept_controller!(integrator::SDEIntegrator, alg::TauLeaping)
    integrator.q = min(integrator.opts.gamma / integrator.EEst, integrator.opts.qmax)
    return integrator.dt * integrator.q
end

function step_reject_controller!(integrator::SDEIntegrator, alg::TauLeaping)
    return integrator.dt = integrator.opts.gamma * integrator.dt / integrator.EEst
end

function stepsize_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping)
    return nothing
end

function step_accept_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping)
    return integrator.EEst # use EEst for the τ
end

function step_reject_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping)
    error("CaoTauLeaping should never reject steps")
end
