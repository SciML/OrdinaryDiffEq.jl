@inline function step_accept_controller!(integrator, alg::TauLeaping, q)
    return step_accept_controller!(integrator, alg)
end
@inline function step_accept_controller!(integrator, alg::CaoTauLeaping, q)
    return step_accept_controller!(integrator, alg)
end

function stepsize_controller!(integrator, alg::TauLeaping)
    return nothing
end

function step_accept_controller!(integrator, alg::TauLeaping)
    q = min(integrator.alg.gamma / integrator.EEst, integrator.alg.qmax)
    return integrator.dt * q
end

function step_reject_controller!(integrator, alg::TauLeaping)
    return integrator.dt = integrator.alg.gamma * integrator.dt / integrator.EEst
end

function stepsize_controller!(integrator, alg::CaoTauLeaping)
    return nothing
end

function step_accept_controller!(integrator, alg::CaoTauLeaping)
    return integrator.EEst # use EEst for the τ
end

function step_reject_controller!(integrator, alg::CaoTauLeaping)
    error("CaoTauLeaping should never reject steps")
end
