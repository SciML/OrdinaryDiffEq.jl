# JVODE
function stepsize_controller!(integrator, alg::JVODE)
    if iszero(integrator.EEst)
        η = alg.qmax
    else
        η = integrator.cache.η
        integrator.cache.ηold = η
    end
    return η
end

function step_accept_controller!(integrator, alg::JVODE, η)
    q = inv(η)
    if q <= alg.qsteady_max && q >= alg.qsteady_min
        q = one(q)
    end
    return integrator.dt / q  # dtnew
end

function step_reject_controller!(integrator, alg::JVODE)
    return integrator.dt *= integrator.cache.η
end
