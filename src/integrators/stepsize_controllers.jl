# PIController stepsize methods now provided by ODE's generic versions
# (stepsize_controller!, step_accept_controller!, step_reject_controller! for PIController)
# which are functionally equivalent for SDE (with better iszero(EEst) handling).

# TauLeaping and CaoTauLeaping use algorithm-specific stepsize logic.
# Bridge: ODE calls step_accept_controller!(integrator, alg, q) with 3 args,
# but TauLeaping's versions take 2 args. These bridges handle the conversion.
@inline function step_accept_controller!(integrator::SDEIntegrator, alg::TauLeaping, q)
    return step_accept_controller!(integrator, alg)
end
@inline function step_accept_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping, q)
    return step_accept_controller!(integrator, alg)
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
