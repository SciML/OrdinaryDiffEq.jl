function step_reject_controller!(integrator, ::DFBDF)
    bdf_step_reject_controller!(integrator, integrator.cache.terkm1)
end

function post_newton_controller!(integrator, alg)
    integrator.dt = integrator.dt / integrator.opts.failfactor
    nothing
end

function post_newton_controller!(integrator, alg::DFBDF)
    @unpack cache = integrator
    if cache.order > 1 && cache.nlsolver.nfails >= 3
        cache.order -= 1
    end
    integrator.dt = integrator.dt / integrator.opts.failfactor
    integrator.cache.consfailcnt += 1
    integrator.cache.nconsteps = 0
    nothing
end

function choose_order!(alg::DFBDF, integrator,
        cache::OrdinaryDiffEqMutableCache,
        ::Val{max_order}) where {max_order}
    @unpack t, dt, u, cache, uprev = integrator
    @unpack atmp, ts_tmp, terkm2, terkm1, terk, terkp1, terk_tmp, u_history = cache
    k = cache.order
    # only when the order of amount of terk follows the order of step size, and achieve enough constant step size, the order could be increased.
    if k < max_order && integrator.cache.nconsteps >= integrator.cache.order + 2 &&
       ((k == 1 && terk > terkp1) ||
        (k == 2 && terkm1 > terk > terkp1) ||
        (k > 2 && terkm2 > terkm1 > terk > terkp1))
        k += 1
        terk = terkp1
    else
        while !(terkm2 > terkm1 > terk > terkp1) && k > 2
            terkp1 = terk
            terk = terkm1
            terkm1 = terkm2
            fd_weights = calc_finite_difference_weights(ts_tmp, t + dt, k - 2,
                Val(max_order))
            terk_tmp = @.. broadcast=false fd_weights[k - 2, 1]*u
            vc = _vec(terk_tmp)
            for i in 2:(k - 2)
                @.. broadcast=false @views vc += fd_weights[i, k - 2] * u_history[:, i - 1]
            end
            @.. broadcast=false terk_tmp*=abs(dt^(k - 2))
            calculate_residuals!(atmp, _vec(terk_tmp), _vec(uprev), _vec(u),
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            terkm2 = integrator.opts.internalnorm(atmp, t)
            k -= 1
        end
    end
    return k, terk
end

function choose_order!(alg::DFBDF, integrator,
        cache::OrdinaryDiffEqConstantCache,
        ::Val{max_order}) where {max_order}
    @unpack t, dt, u, cache, uprev = integrator
    @unpack ts_tmp, terkm2, terkm1, terk, terkp1, u_history = cache
    k = cache.order
    if k < max_order && integrator.cache.nconsteps >= integrator.cache.order + 2 &&
       ((k == 1 && terk > terkp1) ||
        (k == 2 && terkm1 > terk > terkp1) ||
        (k > 2 && terkm2 > terkm1 > terk > terkp1))
        k += 1
        terk = terkp1
    else
        while !(terkm2 > terkm1 > terk > terkp1) && k > 2
            terkp1 = terk
            terk = terkm1
            terkm1 = terkm2
            fd_weights = calc_finite_difference_weights(ts_tmp, t + dt, k - 2,
                Val(max_order))
            terk_tmp = @.. broadcast=false fd_weights[k - 2, 1]*u
            if u isa Number
                for i in 2:(k - 2)
                    terk_tmp += fd_weights[i, k - 2] * u_history[i - 1]
                end
                terk_tmp *= abs(dt^(k - 2))
            else
                vc = _vec(terk_tmp)
                for i in 2:(k - 2)
                    @.. broadcast=false @views vc += fd_weights[i, k - 2] *
                                                     u_history[:, i - 1]
                end
                terk_tmp = reshape(vc, size(terk_tmp))
                terk_tmp *= @.. broadcast=false abs(dt^(k - 2))
            end
            atmp = calculate_residuals(_vec(terk_tmp), _vec(uprev), _vec(u),
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            terkm2 = integrator.opts.internalnorm(atmp, t)
            k -= 1
        end
    end
    return k, terk
end

function stepsize_controller!(integrator,
        alg::DFBDF{max_order}) where {
        max_order,
}
    @unpack cache = integrator
    cache.prev_order = cache.order
    k, terk = choose_order!(alg, integrator, cache, Val(max_order))
    if k != cache.order
        integrator.cache.nconsteps = 0
        cache.order = k
    end
    if iszero(terk)
        q = inv(integrator.opts.qmax)
    else
        q = ((2 * terk / (k + 1))^(1 / (k + 1)))
    end
    integrator.qold = q
    q
end

function step_accept_controller!(integrator, alg::DFBDF{max_order},
        q) where {max_order}
    integrator.cache.consfailcnt = 0
    if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
        q = one(q)
    end
    integrator.cache.nconsteps += 1
    integrator.cache.iters_from_event += 1
    return integrator.dt / q
end