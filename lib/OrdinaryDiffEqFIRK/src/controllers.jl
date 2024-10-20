@inline function stepsize_controller!(integrator, controller::PredictiveController, alg)
    @unpack qmin, qmax, gamma = integrator.opts
    EEst = DiffEqBase.value(integrator.EEst)
    if iszero(EEst)
        q = inv(qmax)
    else
        if fac_default_gamma(alg)
            fac = gamma
        else
            if isfirk(alg)
                @unpack iter = integrator.cache
                @unpack maxiters = alg
            else
                @unpack iter, maxiters = integrator.cache.nlsolver
            end
            fac = min(gamma, (1 + 2 * maxiters) * gamma / (iter + 2 * maxiters))
        end
        expo = 1 / (get_current_adaptive_order(alg, integrator.cache) + 1)
        qtmp = FastPower.fastpower(EEst, expo) / fac
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
               FastPower.fastpower((EEst^2) / integrator.erracc, expo)
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


function step_accept_controller!(integrator, controller::PredictiveController, alg::AdaptiveRadau, q)
    @unpack qmin, qmax, gamma, qsteady_min, qsteady_max = integrator.opts
    @unpack cache = integrator
    @unpack num_stages, step, iter, hist_iter = cache
 
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
    cache.step = step + 1
    hist_iter = hist_iter * 0.8 + iter * 0.2
    cache.hist_iter = hist_iter
    if (step > 10)
        if (hist_iter < 2.6 && num_stages < alg.max_stages)
            cache.num_stages += 2
            cache.step = 1
            cache.hist_iter = iter
        elseif ((hist_iter > 8 || cache.status == VerySlowConvergence || cache.status == Divergence) && num_stages > alg.min_stages)
            cache.num_stages -= 2 
            cache.step = 1
            cache.hist_iter = iter
        end
    end
    return integrator.dt / qacc
end

function step_reject_controller!(integrator, controller::PredictiveController, alg)
    @unpack dt, success_iter, qold = integrator
    integrator.dt = success_iter == 0 ? 0.1 * dt : dt / qold
end

function step_reject_controller!(integrator, controller::PredictiveController, alg::AdaptiveRadau)
    @unpack dt, success_iter, qold = integrator
    @unpack cache = integrator 
    @unpack num_stages, step, iter, hist_iter = cache
    integrator.dt = success_iter == 0 ? 0.1 * dt : dt / qold
    cache.step = step + 1
    hist_iter = hist_iter * 0.8 + iter * 0.2
    cache.hist_iter = hist_iter
    if (step > 10)
        if ((hist_iter > 8 || cache.status == VerySlowConvergence || cache.status == Divergence) && num_stages > alg.min_stages)
            cache.num_stages -= 2 
            cache.step = 1
            cache.hist_iter = iter
        end
    end
end

