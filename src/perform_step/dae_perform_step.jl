function initialize!(integrator, cache::DImplicitEulerConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.k[1] = integrator.du
end

function initialize!(integrator, cache::DImplicitEulerCache)
    integrator.kshortsize = 2
    @unpack k₁, k₂ = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k .= [k₁, k₂]
    integrator.k[1] .= integrator.du
    nothing
end

@muladd function perform_step!(integrator, cache::DImplicitEulerConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack nlsolver = cache

    nlsolver.z = zero(u)
    nlsolver.tmp = zero(u)
    nlsolver.γ = 1
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = uprev + z

    if integrator.opts.adaptive && integrator.success_iter > 0
        # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
        # use 2nd divided differences (DD) a la SPICE and Shampine

        # TODO: check numerical stability
        uprev2 = integrator.uprev2
        tprev = integrator.tprev

        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
        r = c * dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

        tmp = r *
              integrator.opts.internalnorm.((u - uprev) / dt1 - (uprev - uprev2) / dt2, t)
        atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    else
        integrator.EEst = 1
    end

    integrator.u = u
    integrator.du = (u - uprev) / dt

    if integrator.opts.calck
        integrator.k[2] = integrator.k[1]
        integrator.k[1] = integrator.du
    end
end

@muladd function perform_step!(integrator, cache::DImplicitEulerCache, repeat_step = false)
    @unpack t, dt, uprev, du, u, f, p = integrator
    @unpack atmp, nlsolver = cache
    @unpack tmp = nlsolver

    @. nlsolver.z = false
    @. nlsolver.tmp = false
    nlsolver.γ = 1
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast=false u=uprev + z
    @.. broadcast=false du=z * inv(dt)

    if integrator.opts.adaptive && integrator.success_iter > 0
        # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
        # use 2nd divided differences (DD) a la SPICE and Shampine

        # TODO: check numerical stability
        uprev2 = integrator.uprev2
        tprev = integrator.tprev

        dt1 = dt * (t + dt - tprev)
        dt2 = (t - tprev) * (t + dt - tprev)
        c = 7 / 12 # default correction factor in SPICE (LTE overestimated by DD)
        r = c * dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

        @.. broadcast=false tmp=r * integrator.opts.internalnorm((u - uprev) / dt1 -
                                                             (uprev - uprev2) / dt2, t)
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    else
        integrator.EEst = 1
    end

    if integrator.opts.calck
        integrator.k[2] .= integrator.k[1]
        integrator.k[1] .= du
    end
end

function initialize!(integrator, cache::DABDF2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.k[1] = integrator.du
end

@muladd function perform_step!(integrator, cache::DABDF2ConstantCache, repeat_step = false)
    @unpack t, f, p = integrator
    @unpack dtₙ₋₁, nlsolver = cache
    dtₙ, uₙ, uₙ₋₁, uₙ₋₂ = integrator.dt, integrator.u, integrator.uprev, integrator.uprev2

    if integrator.iter == 1 && !integrator.u_modified
        cache.dtₙ₋₁ = dtₙ
        perform_step!(integrator, cache.eulercache, repeat_step)
        integrator.fsalfirst = @.. broadcast=false (integrator.u - integrator.uprev)/dtₙ
        cache.fsalfirstprev = integrator.fsalfirst
        return
    end

    # precalculations
    ρ = dtₙ / dtₙ₋₁
    c1 = ρ^2 / (1 + 2ρ)

    nlsolver.γ = (1 + ρ) / (1 + 2ρ)
    nlsolver.α = 1 // 1

    nlsolver.z = zero(uₙ)

    nlsolver.tmp = -c1 * uₙ₋₁ + c1 * uₙ₋₂
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    uₙ = uₙ₋₁ + z
    integrator.fsallast = @.. broadcast=false z/dtₙ

    if integrator.opts.adaptive
        tmp = integrator.fsallast - (1 + dtₙ / dtₙ₋₁) * integrator.fsalfirst +
              (dtₙ / dtₙ₋₁) * cache.fsalfirstprev
        est = (dtₙ₋₁ + dtₙ) / 6 * tmp
        atmp = calculate_residuals(est, uₙ₋₁, uₙ, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    ################################### Finalize

    if integrator.EEst < one(integrator.EEst)
        cache.fsalfirstprev = integrator.fsalfirst
        cache.dtₙ₋₁ = dtₙ
    end

    integrator.u = uₙ
    integrator.du = du = (nlsolver.α * z + nlsolver.tmp) * inv(nlsolver.γ * dtₙ)

    if integrator.opts.calck
        integrator.k[2] = integrator.k[1]
        integrator.k[1] = integrator.du
    end
    return
end

function initialize!(integrator, cache::DABDF2Cache)
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)

    integrator.kshortsize = 2
    @unpack k₁, k₂ = cache.eulercache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k .= [k₁, k₂]
    integrator.k[1] .= integrator.du
    nothing
end

@muladd function perform_step!(integrator, cache::DABDF2Cache, repeat_step = false)
    @unpack t, dt, du, f, p = integrator
    @unpack atmp, dtₙ₋₁, nlsolver = cache
    @unpack z, tmp = nlsolver
    uₙ, uₙ₋₁, uₙ₋₂, dtₙ = integrator.u, integrator.uprev, integrator.uprev2, integrator.dt

    if integrator.iter == 1 && !integrator.u_modified
        cache.dtₙ₋₁ = dtₙ
        perform_step!(integrator, cache.eulercache, repeat_step)
        @.. broadcast=false integrator.fsalfirst=(uₙ - uₙ₋₁) / dt
        cache.fsalfirstprev .= integrator.fsalfirst
        return
    end

    # precalculations
    ρ = dtₙ / dtₙ₋₁
    c1 = ρ^2 / (1 + 2ρ)

    nlsolver.γ = (1 + ρ) / (1 + 2ρ)
    nlsolver.α = 1 // 1
    @.. broadcast=false nlsolver.tmp=-c1 * uₙ₋₁ + c1 * uₙ₋₂
    nlsolver.z .= zero(eltype(z))
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast=false uₙ=uₙ₋₁ + z
    @.. broadcast=false du=(nlsolver.α * z + nlsolver.tmp) * inv(nlsolver.γ * dt)

    @.. broadcast=false integrator.fsallast=du
    integrator.stats.nf += 1
    if integrator.opts.adaptive
        btilde0 = (dtₙ₋₁ + dtₙ) * 1 // 6
        btilde1 = 1 + ρ
        btilde2 = ρ
        @.. broadcast=false tmp=btilde0 *
                                (integrator.fsallast - btilde1 * integrator.fsalfirst +
                                 btilde2 * cache.fsalfirstprev)
        calculate_residuals!(atmp, tmp, uₙ₋₁, uₙ, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    ################################### Finalize

    if integrator.EEst < one(integrator.EEst)
        @.. broadcast=false cache.fsalfirstprev=integrator.fsalfirst
        cache.dtₙ₋₁ = dtₙ
    end
    return
end

function initialize!(integrator, cache::DFBDFConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.du, integrator.uprev, integrator.p,
        integrator.t) # Pre-start fsal
    integrator.stats.nf += 1
    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::DFBDFConstantCache{max_order},
        repeat_step = false) where {max_order}
    @unpack ts, u_history, order, u_corrector, bdf_coeffs, r, nlsolver, weights, ts_tmp, iters_from_event, nconsteps = cache
    @unpack t, dt, u, f, p, uprev = integrator

    k = order
    reinitFBDF!(integrator, cache)
    cache.u₀ = zero(u)
    if iters_from_event >= 1
        cache.u₀ = calc_Lagrange_interp(k, weights, t + dt, ts, u_history, cache.u₀)
    else
        cache.u₀ = u
    end
    markfirststage!(nlsolver)
    equi_ts = zeros(k - 1)
    for i in 1:(k - 1)
        equi_ts[i] = t - dt * i
    end

    fill!(u_corrector, zero(eltype(u)))
    if u isa Number
        for i in 1:(k - 1)
            u_corrector[i] = calc_Lagrange_interp(k, weights, equi_ts[i], ts, u_history,
                u_corrector[i])
        end
        tmp = uprev * bdf_coeffs[k, 2]
        for i in 1:(k - 1)
            tmp += u_corrector[i] * bdf_coeffs[k, i + 2]
        end
    else
        for i in 1:(k - 1)
            @.. broadcast=false @views u_corrector[:, i] = $calc_Lagrange_interp(k, weights,
                equi_ts[i],
                ts,
                u_history,
                u_corrector[:,
                    i])
        end
        tmp = uprev * bdf_coeffs[k, 2]
        vc = _vec(tmp)
        for i in 1:(k - 1)
            vc += @.. broadcast=false u_corrector[:, i]*bdf_coeffs[k, i + 2]
        end
        tmp = reshape(vc, size(tmp))
    end

    nlsolver.tmp = tmp + cache.u₀
    nlsolver.z = zero(nlsolver.z)
    nlsolver.γ = bdf_coeffs[k, 1]
    nlsolver.α = 1 // 1
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = z + cache.u₀

    for j in 2:k
        r[j] = (1 - j)
        for i in 2:(k + 1)
            r[j] *= ((t + dt - j * dt) - ts[i]) / (i * dt)
        end
    end

    terkp1 = z
    for j in 1:(k + 1)
        terkp1 *= j * dt / (t + dt - ts[j])
    end

    lte = -1 / (1 + k)
    for j in 2:k
        lte -= (bdf_coeffs[k, j] // bdf_coeffs[k, 1]) * r[j]
    end
    lte *= terkp1

    if integrator.opts.adaptive
        for i in 1:(k + 1)
            ts_tmp[i + 1] = ts[i]
        end
        ts_tmp[1] = t + dt
        atmp = calculate_residuals(_vec(lte), _vec(uprev), _vec(u), integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        terk = estimate_terk(integrator, cache, k + 1, Val(max_order), u)
        atmp = calculate_residuals(_vec(terk), _vec(uprev), _vec(u), integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        cache.terk = integrator.opts.internalnorm(atmp, t)

        if k > 1
            terkm1 = estimate_terk(integrator, cache, k, Val(max_order), u)
            atmp = calculate_residuals(_vec(terkm1), _vec(uprev), _vec(u),
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            cache.terkm1 = integrator.opts.internalnorm(atmp, t)
        end
        if k > 2
            terkm2 = estimate_terk(integrator, cache, k - 1, Val(max_order), u)
            atmp = calculate_residuals(_vec(terkm2), _vec(uprev), _vec(u),
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            cache.terkm2 = integrator.opts.internalnorm(atmp, t)
        end
        if nconsteps > k + 1 && k < max_order
            atmp = calculate_residuals(_vec(terkp1), _vec(uprev), _vec(u),
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            cache.terkp1 = integrator.opts.internalnorm(atmp, t)
        else
            cache.terkp1 = zero(cache.terk)
        end
    end
    integrator.u = u
    integrator.fsallast = integrator.du = (nlsolver.α * z + nlsolver.tmp) *
                                          inv(nlsolver.γ * dt)
end

function initialize!(integrator, cache::DFBDFCache)
    integrator.kshortsize = 2
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    #integrator.f(integrator.fsalfirst, integrator.du, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    #integrator.stats.nf += 1
end

function perform_step!(integrator, cache::DFBDFCache{max_order},
        repeat_step = false) where {max_order}
    @unpack ts, u_history, order, u_corrector, bdf_coeffs, r, nlsolver, weights, terk_tmp, terkp1_tmp, atmp, tmp, equi_ts, u₀, ts_tmp = cache
    @unpack t, dt, u, f, p, uprev = integrator

    reinitFBDF!(integrator, cache)
    k = order
    @.. broadcast=false u₀=zero(u)
    if cache.iters_from_event >= 1
        calc_Lagrange_interp!(k, weights, t + dt, ts, u_history, u₀)
    else
        @.. broadcast=false u₀=u
    end
    markfirststage!(nlsolver)

    for i in 1:(k - 1)
        equi_ts[i] = t - dt * i
    end

    fill!(u_corrector, zero(eltype(u)))
    for i in 1:(k - 1)
        @views calc_Lagrange_interp!(k, weights, equi_ts[i], ts, u_history,
            u_corrector[:, i])
    end

    @.. broadcast=false tmp=uprev * bdf_coeffs[k, 2]
    vc = _vec(tmp)
    for i in 1:(k - 1)
        @.. broadcast=false @views vc += u_corrector[:, i] * bdf_coeffs[k, i + 2]
    end

    @.. broadcast=false nlsolver.tmp=tmp + u₀
    @.. broadcast=false nlsolver.z=zero(eltype(nlsolver.z))
    nlsolver.γ = bdf_coeffs[k, 1]
    nlsolver.α = 1 // 1
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast=false u=z + u₀

    for j in 2:k
        r[j] = (1 - j)
        for i in 2:(k + 1)
            r[j] *= ((t + dt - j * dt) - ts[i]) / (i * dt)
        end
    end

    @.. broadcast=false terkp1_tmp=z
    for j in 1:(k + 1)
        @.. broadcast=false terkp1_tmp*=j * dt / (t + dt - ts[j])
    end

    lte = -1 / (1 + k)
    for j in 2:k
        lte -= (bdf_coeffs[k, j] // bdf_coeffs[k, 1]) * r[j]
    end
    @.. broadcast=false terk_tmp=lte * terkp1_tmp
    if integrator.opts.adaptive
        @unpack abstol, reltol, internalnorm = integrator.opts
        for i in 1:(k + 1)
            ts_tmp[i + 1] = ts[i]
        end
        ts_tmp[1] = t + dt
        calculate_residuals!(atmp, _vec(terk_tmp), _vec(uprev), _vec(u), abstol, reltol,
            internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
        estimate_terk!(integrator, cache, k + 1, Val(max_order))
        calculate_residuals!(atmp, _vec(terk_tmp), _vec(uprev), _vec(u), abstol, reltol,
            internalnorm, t)
        cache.terk = integrator.opts.internalnorm(atmp, t)

        if k > 1
            estimate_terk!(integrator, cache, k, Val(max_order))
            calculate_residuals!(atmp, _vec(terk_tmp), _vec(uprev), _vec(u),
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            cache.terkm1 = integrator.opts.internalnorm(atmp, t)
        end
        if k > 2
            estimate_terk!(integrator, cache, k - 1, Val(max_order))
            calculate_residuals!(atmp, _vec(terk_tmp), _vec(uprev), _vec(u),
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            cache.terkm2 = integrator.opts.internalnorm(atmp, t)
        end
        if cache.nconsteps > k + 1 && k < max_order
            calculate_residuals!(atmp, _vec(terkp1_tmp), _vec(uprev), _vec(u),
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            cache.terkp1 = integrator.opts.internalnorm(atmp, t)
        else
            cache.terkp1 = zero(cache.terkp1)
        end
    end
    @.. broadcast=false integrator.fsallast=integrator.du = (nlsolver.α * z +
                                                             nlsolver.tmp) *
                                                            inv(nlsolver.γ * dt) #TODO Lorenz plot seems not smooth
end
