struct DummyController <: AbstractController
end

# JVODE
function stepsize_controller!(integrator, alg::JVODE)
    if iszero(integrator.EEst)
        η = integrator.opts.qmax
    else
        η = integrator.cache.η
        integrator.qold = η
    end
    return η
end

function step_accept_controller!(integrator, alg::JVODE, η)
    q = inv(η)
    if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
        q = one(q)
    end
    return integrator.dt / q  # dtnew
end

function step_reject_controller!(integrator, alg::JVODE)
    return integrator.dt *= integrator.qold
end
