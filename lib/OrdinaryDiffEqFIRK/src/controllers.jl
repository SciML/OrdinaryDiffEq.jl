@inline function stepsize_controller!(integrator, controller::PredictiveController, alg)
    @unpack qmin, qmax, gamma = integrator.opts
    EEst = DiffEqBase.value(integrator.EEst)

    if iszero(EEst)
        q = inv(qmax)
    else
        if fac_default_gamma(alg)
            fac = gamma
        else
            if alg isa Union{RadauIIA3, RadauIIA5}
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