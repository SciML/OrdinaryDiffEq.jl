@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.4"
    @eval begin
        function legacy_default_controller(alg::Union{QNDF, FBDF}, args...)
            return DummyController()
        end

        function default_controller_v7(QT, alg::Union{QNDF, FBDF}, args...)
            return DummyController()
        end
    end
else
    @eval begin
        function default_controller(alg::Union{QNDF, FBDF}, args...)
            return DummyController()
        end
    end
end

# QNBDF
stepsize_controller!(integrator, alg::QNDF) = nothing

# this stepsize and order controller is taken from
# Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s paper
# E. Alberdi Celaya, J. J. Anza Aguirrezabala, and P. Chatzipantelidis

function step_accept_controller!(integrator, alg::QNDF, q)
    return step_accept_controller!(integrator, integrator.cache, alg, q)
end

function step_accept_controller!(integrator, cache::Union{QNDFCache, QNDFConstantCache}, alg::QNDF{max_order}, q) where {max_order}
    #step is accepted, reset count of consecutive failed steps
    integrator.cache.consfailcnt = 0
    integrator.cache.nconsteps += 1
    if iszero(integrator.EEst)
        return integrator.dt * get_current_qmax(integrator, integrator.opts.qmax)
    else
        est = integrator.EEst
        estₖ₋₁ = integrator.cache.EEst1
        estₖ₊₁ = integrator.cache.EEst2
        h = integrator.dt
        k = integrator.cache.order
        cache = integrator.cache
        prefer_const_step = integrator.cache.nconsteps < integrator.cache.order + 2
        zₛ = 1.2 # equivalent to integrator.opts.gamma
        zᵤ = 0.1
        Fᵤ = 10
        expo = 1 / (k + 1)
        z = zₛ * ((est)^expo)
        F = inv(z)
        hₙ = h
        kₙ = k
        if z <= zᵤ
            hₖ = Fᵤ * h
        else
            hₖ = F * h
        end
        hₖ₋₁ = 0.0
        hₖ₊₁ = 0.0

        if k > 1
            expo = 1 / k
            zₖ₋₁ = 1.3 * ((estₖ₋₁)^expo)
            Fₖ₋₁ = inv(zₖ₋₁)
            if zₖ₋₁ <= 0.1
                hₖ₋₁ = 10 * h
            elseif zₖ₋₁ <= 1.3
                hₖ₋₁ = Fₖ₋₁ * h
            end
            if hₖ₋₁ > hₖ
                hₙ = hₖ₋₁
                kₙ = k - 1
            else
                hₙ = hₖ
                kₙ = k
            end
        else
            hₙ = hₖ
            kₙ = k
        end

        if k < max_order
            expo = 1 / (k + 2)
            zₖ₊₁ = 1.4 * ((estₖ₊₁)^expo)
            Fₖ₊₁ = inv(zₖ₊₁)

            if zₖ₊₁ <= 0.1
                hₖ₊₁ = 10 * h
            elseif 0.1 < zₖ₊₁ <= 1.4
                hₖ₊₁ = Fₖ₊₁ * h
            end
            if hₖ₊₁ > hₙ
                hₙ = hₖ₊₁
                kₙ = k + 1
            end
        end
        cache.order = kₙ
        q = integrator.dt / hₙ
    end
    if prefer_const_step
        if q < 1.2 && q > 0.6
            return integrator.dt
        end
    end
    if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
        return integrator.dt
    end
    return integrator.dt / q
end

function bdf_step_reject_controller!(integrator, EEst1)
    k = integrator.cache.order
    h = integrator.dt
    integrator.cache.consfailcnt += 1
    integrator.cache.nconsteps = 0
    if integrator.cache.consfailcnt > 1
        h = h / 2
    end
    zₛ = 1.2  # equivalent to integrator.opts.gamma
    expo = 1 / (k + 1)
    z = zₛ * ((integrator.EEst)^expo)
    F = inv(z)
    if z <= 10
        hₖ = F * h
    else # z > 10
        hₖ = 0.1 * h
    end
    hₙ = hₖ
    kₙ = k
    if k > 1
        expo = 1 / k
        zₖ₋₁ = 1.3 * (EEst1^expo)
        Fₖ₋₁ = inv(zₖ₋₁)
        if zₖ₋₁ <= 10
            hₖ₋₁ = Fₖ₋₁ * h
        else # zₖ₋₁ > 10
            hₖ₋₁ = 0.1 * h
        end
        if integrator.cache.consfailcnt > 2 || hₖ₋₁ > hₖ
            hₙ = min(h, hₖ₋₁)
            kₙ = k - 1
        end
    end
    # Restart BDf (clear history) when we failed repeatedly
    if kₙ == 1 && integrator.cache.consfailcnt > 3
        u_modified!(integrator, true)
    end
    integrator.dt = hₙ
    return integrator.cache.order = kₙ
end

function step_reject_controller!(integrator, ::QNDF)
    return bdf_step_reject_controller!(integrator, integrator.cache.EEst1)
end

function step_reject_controller!(integrator, ::FBDF)
    return bdf_step_reject_controller!(integrator, integrator.cache.terkm1)
end

function post_newton_controller!(integrator, alg::FBDF)
    (; cache) = integrator
    if cache.order > 1 && cache.nlsolver.nfails >= 3
        cache.order -= 1
    end
    integrator.dt = integrator.dt / integrator.opts.failfactor
    integrator.cache.consfailcnt += 1
    integrator.cache.nconsteps = 0
    return nothing
end

function choose_order!(
        alg::FBDF, integrator,
        cache::OrdinaryDiffEqMutableCache,
        ::Val{max_order}
    ) where {max_order}
    (; t, dt, u, cache, uprev) = integrator
    (; atmp, ts_tmp, terkm2, terkm1, terk, terkp1, terk_tmp, u_history, fd_weights) = cache
    k = cache.order
    # Use CVODE-style qwait countdown: only consider order increase when qwait reaches 0
    if k < max_order && integrator.cache.qwait == 0 &&
            (
            (k == 1 && terk > terkp1) ||
                (k == 2 && terkm1 > terk > terkp1) ||
                (k > 2 && terkm2 > terkm1 > terk > terkp1)
        )
        k += 1
        terk = terkp1
    else
        while !(terkm2 > terkm1 > terk > terkp1) && k > 2
            terkp1 = terk
            terk = terkm1
            terkm1 = terkm2
            calc_finite_difference_weights!(
                fd_weights, ts_tmp, t + dt, k - 2
            )
            @.. broadcast = false terk_tmp = fd_weights[k - 2, 1] * u
            for i in 2:(k - 2)
                @.. broadcast = false terk_tmp += fd_weights[i, k - 2] * u_history[i - 1]
            end
            @.. broadcast = false terk_tmp *= abs(dt^(k - 2))
            calculate_residuals!(
                atmp, terk_tmp, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            terkm2 = integrator.opts.internalnorm(atmp, t)
            k -= 1
        end
    end
    return k, terk
end

function choose_order!(
        alg::FBDF, integrator,
        cache::OrdinaryDiffEqConstantCache,
        ::Val{max_order}
    ) where {max_order}
    (; t, dt, u, cache, uprev) = integrator
    (; ts_tmp, terkm2, terkm1, terk, terkp1, u_history) = cache
    k = cache.order
    if k < max_order && integrator.cache.qwait == 0 &&
            (
            (k == 1 && terk > terkp1) ||
                (k == 2 && terkm1 > terk > terkp1) ||
                (k > 2 && terkm2 > terkm1 > terk > terkp1)
        )
        k += 1
        terk = terkp1
    else
        while !(terkm2 > terkm1 > terk > terkp1) && k > 2
            terkp1 = terk
            terk = terkm1
            terkm1 = terkm2
            fd_weights = calc_finite_difference_weights(
                ts_tmp, t + dt, k - 2,
                Val(max_order)
            )
            local terk_tmp
            if u isa Number
                terk_tmp = fd_weights[k - 2, 1] * u
                for i in 2:(k - 2)
                    terk_tmp += fd_weights[i, k - 2] * u_history[i - 1]
                end
                terk_tmp *= abs(dt^(k - 2))
            else
                # we need terk_tmp to be mutable.
                # so it can be updated
                terk_tmp = fd_weights[k - 2, 1] * u
                for i in 2:(k - 2)
                    terk_tmp = @.. terk_tmp + fd_weights[i, k - 2] * u_history[i - 1]
                end
                terk_tmp = @.. terk_tmp * abs(dt^(k - 2))
            end
            atmp = calculate_residuals(
                terk_tmp, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            terkm2 = integrator.opts.internalnorm(atmp, t)
            k -= 1
        end
    end
    return k, terk
end

function stepsize_controller!(
        integrator,
        alg::FBDF
    )
    return stepsize_controller!(integrator, integrator.cache, alg)
end

function stepsize_controller!(
        integrator,
        cache::Union{FBDFCache, FBDFConstantCache},
        alg::FBDF{max_order}
    ) where {
        max_order,
    }
    cache.prev_order = cache.order
    k, terk = choose_order!(alg, integrator, cache, Val(max_order))
    if k != cache.order
        cache.nconsteps = 0
        cache.order = k
    end
    if iszero(terk)
        q = inv(get_current_qmax(integrator, integrator.opts.qmax))
    else
        q = ((2 * terk / (k + 1))^(1 / (k + 1)))
    end
    integrator.qold = q
    return q
end

function step_accept_controller!(integrator, alg::FBDF, q)
    return step_accept_controller!(integrator, integrator.cache, alg, q)
end

function step_accept_controller!(
        integrator, cache::Union{FBDFCache, FBDFConstantCache}, alg::FBDF{max_order},
        q
    ) where {max_order}
    cache.consfailcnt = 0
    if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
        q = one(q)
    end
    cache.nconsteps += 1
    cache.iters_from_event += 1
    # CVODE-style qwait countdown for order change gating
    if cache.order != cache.prev_order
        cache.qwait = cache.order + 2 # reset after order change, matching nconsteps >= order + 2
    elseif cache.qwait > 0
        cache.qwait -= 1 # countdown
    end
    return integrator.dt / q
end

function step_reject_controller!(integrator, ::DFBDF)
    return bdf_step_reject_controller!(integrator, integrator.cache.terkm1)
end

function post_newton_controller!(integrator, alg::DFBDF)
    (; cache) = integrator
    if cache.order > 1 && cache.nlsolver.nfails >= 3
        cache.order -= 1
    end
    integrator.dt = integrator.dt / integrator.opts.failfactor
    integrator.cache.consfailcnt += 1
    integrator.cache.nconsteps = 0
    return nothing
end

function choose_order!(
        alg::DFBDF, integrator,
        cache::OrdinaryDiffEqMutableCache,
        ::Val{max_order}
    ) where {max_order}
    (; t, dt, u, cache, uprev) = integrator
    (; atmp, ts_tmp, terkm2, terkm1, terk, terkp1, terk_tmp, u_history, fd_weights) = cache
    k = cache.order
    # Use CVODE-style qwait countdown: only consider order increase when qwait reaches 0
    if k < max_order && integrator.cache.qwait == 0 &&
            (
            (k == 1 && terk > terkp1) ||
                (k == 2 && terkm1 > terk > terkp1) ||
                (k > 2 && terkm2 > terkm1 > terk > terkp1)
        )
        k += 1
        terk = terkp1
    else
        while !(terkm2 > terkm1 > terk > terkp1) && k > 2
            terkp1 = terk
            terk = terkm1
            terkm1 = terkm2
            calc_finite_difference_weights!(
                fd_weights, ts_tmp, t + dt, k - 2
            )
            @.. broadcast = false terk_tmp = fd_weights[k - 2, 1] * u
            for i in 2:(k - 2)
                @.. broadcast = false terk_tmp += fd_weights[i, k - 2] * u_history[i - 1]
            end
            @.. broadcast = false terk_tmp *= abs(dt^(k - 2))
            calculate_residuals!(
                atmp, terk_tmp, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            terkm2 = integrator.opts.internalnorm(atmp, t)
            k -= 1
        end
    end
    return k, terk
end

function choose_order!(
        alg::DFBDF, integrator,
        cache::OrdinaryDiffEqConstantCache,
        ::Val{max_order}
    ) where {max_order}
    (; t, dt, u, cache, uprev) = integrator
    (; ts_tmp, terkm2, terkm1, terk, terkp1, u_history) = cache
    k = cache.order
    if k < max_order && integrator.cache.qwait == 0 &&
            (
            (k == 1 && terk > terkp1) ||
                (k == 2 && terkm1 > terk > terkp1) ||
                (k > 2 && terkm2 > terkm1 > terk > terkp1)
        )
        k += 1
        terk = terkp1
    else
        while !(terkm2 > terkm1 > terk > terkp1) && k > 2
            terkp1 = terk
            terk = terkm1
            terkm1 = terkm2
            fd_weights = calc_finite_difference_weights(
                ts_tmp, t + dt, k - 2,
                Val(max_order)
            )
            terk_tmp = @.. broadcast = false fd_weights[k - 2, 1] * u
            if u isa Number
                for i in 2:(k - 2)
                    terk_tmp += fd_weights[i, k - 2] * u_history[i - 1]
                end
                terk_tmp *= abs(dt^(k - 2))
            else
                for i in 2:(k - 2)
                    terk_tmp = @.. terk_tmp + fd_weights[i, k - 2] * u_history[i - 1]
                end
                terk_tmp = @.. broadcast = false terk_tmp * abs(dt^(k - 2))
            end
            atmp = calculate_residuals(
                terk_tmp, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            terkm2 = integrator.opts.internalnorm(atmp, t)
            k -= 1
        end
    end
    return k, terk
end

function stepsize_controller!(
        integrator,
        alg::DFBDF
    )
    return stepsize_controller!(integrator, integrator.cache, alg)
end

function stepsize_controller!(
        integrator,
        cache::Union{DFBDFCache, DFBDFConstantCache},
        alg::DFBDF{max_order}
    ) where {
        max_order,
    }
    cache.prev_order = cache.order
    k, terk = choose_order!(alg, integrator, cache, Val(max_order))
    if k != cache.order
        cache.nconsteps = 0
        cache.order = k
    end
    if iszero(terk)
        q = inv(get_current_qmax(integrator, integrator.opts.qmax))
    else
        q = ((2 * terk / (k + 1))^(1 / (k + 1)))
    end
    integrator.qold = q
    return q
end

function step_accept_controller!(integrator, alg::DFBDF, q)
    return step_accept_controller!(integrator, integrator.cache, alg, q)
end

function step_accept_controller!(
        integrator, cache::Union{DFBDFCache, DFBDFConstantCache}, alg::DFBDF{max_order},
        q
    ) where {max_order}
    cache.consfailcnt = 0
    if q <= integrator.opts.qsteady_max && q >= integrator.opts.qsteady_min
        q = one(q)
    end
    cache.nconsteps += 1
    cache.iters_from_event += 1
    # CVODE-style qwait countdown for order change gating
    if cache.order != cache.prev_order
        cache.qwait = cache.order + 2 # reset after order change, matching nconsteps >= order + 2
    elseif cache.qwait > 0
        cache.qwait -= 1 # countdown
    end
    return integrator.dt / q
end
