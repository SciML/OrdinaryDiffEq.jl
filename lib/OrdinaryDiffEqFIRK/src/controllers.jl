function step_accept_controller!(
        integrator, controller::PredictiveController, alg::AdaptiveRadau, q)
    @unpack qmin, qmax, gamma, qsteady_min, qsteady_max = integrator.opts
    @unpack cache = integrator
    @unpack num_stages, step, iter, hist_iter, index = cache

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
    max_stages = (alg.max_order - 1) ÷ 4 * 2 + 1
    min_stages = (alg.min_order - 1) ÷ 4 * 2 + 1
    if (step > 10)
        if (hist_iter < 2.6 && num_stages < max_stages)
            cache.num_stages += 2
            cache.index += 1
            cache.step = 1
            cache.hist_iter = iter
        elseif ((hist_iter > 8 || cache.status == VerySlowConvergence ||
                 cache.status == Divergence) && num_stages > min_stages)
            cache.num_stages -= 2
            cache.index -= 1
            cache.step = 1
            cache.hist_iter = iter
        end
    end
    return integrator.dt / qacc
end

function step_reject_controller!(
        integrator, controller::PredictiveController, alg::AdaptiveRadau)
    @unpack dt, success_iter, qold = integrator
    @unpack cache = integrator
    @unpack num_stages, step, iter, hist_iter = cache
    integrator.dt = success_iter == 0 ? 0.1 * dt : dt / qold
    cache.step = step + 1
    hist_iter = hist_iter * 0.8 + iter * 0.2
    cache.hist_iter = hist_iter
    min_stages = (alg.min_order - 1) ÷ 4 * 2 + 1
    if (step > 10)
        if ((hist_iter > 8 || cache.status == VerySlowConvergence ||
             cache.status == Divergence) && num_stages > min_stages)
            cache.num_stages -= 2
            cache.index -= 1
            cache.step = 1
            cache.hist_iter = iter
        end
    end
end
