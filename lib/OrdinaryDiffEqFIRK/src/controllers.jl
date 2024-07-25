@inline function stepsize_controller!(integrator, controller::PredictiveController, alg)
    @unpack qmin, qmax, gamma = integrator.opts
    EEst = DiffEqBase.value(integrator.EEst)

    if iszero(EEst)
        q = inv(qmax)
    else
        if fac_default_gamma(alg)
            fac = gamma
        else
            if isfirk(alg::Union{RadauIIA3, RadauIIA5, RadauIIA7})
                @unpack iter = integrator.cache
                @unpack maxiters = alg
            else
                @unpack iter, maxiters = integrator.cache.nlsolver
            end
            fac = min(gamma, (1 + 2 * maxiters) * gamma / (iter + 2 * maxiters))
        end
        expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
        qtmp = DiffEqBase.fastpow(EEst, expo) / fac
        @fastmath q = DiffEqBase.value(max(inv(qmax), min(inv(qmin), qtmp)))
        integrator.qold = q
    end
    q
end

function step_accept_controller!(integrator, controller::PredictiveController, alg, q)
    @unpack qmin, qmax, gamma, qsteady_min, qsteady_max = integrator.opts
    EEst = DiffEqBase.value(integrator.EEst)

    if integrator.success_iter > 0
        expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
        qgus = (integrator.dtacc / integrator.dt) *
               DiffEqBase.fastpow((EEst^2) / integrator.erracc, expo)
        qgus = max(inv(qmax), min(inv(qmin), qgus / gamma))
        qacc = max(q, qgus)
    else
        qacc = q
    end
    if qsteady_min <= qacc <= qsteady_max
        qacc = one(qacc)
    end
    integrator.dtacc = integrator.dt
    integrator.erracc = max(1e-2, EEst)
    return integrator.dt / qacc
end

function step_reject_controller!(integrator, controller::PredictiveController, alg)
    @unpack dt, success_iter, qold = integrator
    integrator.dt = success_iter == 0 ? 0.1 * dt : dt / qold
end
