function initialize!(integrator, cache::ROCK2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    alg = unwrap_alg(integrator, true)
    cache.max_stage = (alg.max_stages < 1 || alg.max_stages > 200) ? 200 : alg.max_stages
    cache.min_stage = (alg.min_stages > cache.max_stage) ? cache.max_stage : alg.min_stages
    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ROCK2ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; ms, fp1, fp2, recf) = cache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    # The the number of degree for Chebyshev polynomial
    mdeg = floor(Int, sqrt((1.5 + abs(dt) * integrator.eigen_est) / 0.811)) + 1
    mdeg = min(max(mdeg, cache.min_stage), cache.max_stage)
    cache.mdeg = max(mdeg, 3) - 2
    choosedeg!(cache)
    # recurrence
    # for the first stage
    tᵢ₋₁ = t + dt * recf[cache.start]
    tᵢ₋₂ = t + dt * recf[cache.start]
    tᵢ₋₃ = t
    uᵢ₋₂ = copy(uprev)
    uᵢ₋₁ = uprev + (dt * recf[cache.start]) * fsalfirst
    cache.mdeg < 2 && (u = uᵢ₋₁)
    # for the second to the ms[cache.mdeg] th stages
    for i in 2:(cache.mdeg)
        μ, κ = recf[cache.start + (i - 2) * 2 + 1], recf[cache.start + (i - 2) * 2 + 2]
        ν = -1 - κ
        u = f(uᵢ₋₁, p, tᵢ₋₁)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tᵢ₋₁ = dt * μ - ν * tᵢ₋₂ - κ * tᵢ₋₃
        u = (dt * μ) * u - ν * uᵢ₋₁ - κ * uᵢ₋₂
        if i < cache.mdeg
            uᵢ₋₂ = uᵢ₋₁
            uᵢ₋₁ = u
        end
        tᵢ₋₃ = tᵢ₋₂
        tᵢ₋₂ = tᵢ₋₁
    end # end if
    # two-stage finishing procedure.
    δt₁ = dt * fp1[cache.deg_index]
    δt₂ = dt * fp2[cache.deg_index]
    uᵢ₋₂ = f(u, p, tᵢ₋₁)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    uᵢ₋₁ = u + δt₁ * uᵢ₋₂
    tᵢ₋₁ += δt₁
    u = f(uᵢ₋₁, p, tᵢ₋₁)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    tmp = u  # Initialize for JET
    if integrator.opts.adaptive
        tmp = δt₂ * (u - uᵢ₋₂)
        u = uᵢ₋₁ + δt₁ * u + tmp
    else
        u = uᵢ₋₁ + δt₁ * u + δt₂ * (u - uᵢ₋₂)
    end
    # error estimate
    if integrator.opts.adaptive
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

function initialize!(integrator, cache::ROCK2Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    alg = unwrap_alg(integrator, true)
    cache.constantcache.max_stage = (alg.max_stages < 1 || alg.max_stages > 200) ? 200 :
        alg.max_stages
    cache.constantcache.min_stage = (alg.min_stages > cache.constantcache.max_stage) ?
        cache.constantcache.max_stage : alg.min_stages

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::ROCK2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; k, tmp, uᵢ₋₂, uᵢ₋₁, atmp) = cache
    (; ms, fp1, fp2, recf) = cache.constantcache
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    # The the number of degree for Chebyshev polynomial
    mdeg = floor(Int, sqrt((1.5 + abs(dt) * integrator.eigen_est) / 0.811)) + 1
    mdeg = min(max(mdeg, ccache.min_stage), ccache.max_stage)
    ccache.mdeg = max(mdeg, 3) - 2
    choosedeg!(cache)
    # recurrence
    # for the first stage
    tᵢ₋₁ = t + dt * recf[ccache.start]
    tᵢ₋₂ = t + dt * recf[ccache.start]
    tᵢ₋₃ = t
    @.. broadcast = false uᵢ₋₂ = uprev
    @.. broadcast = false uᵢ₋₁ = uprev + (dt * recf[ccache.start]) * fsalfirst
    ccache.mdeg < 2 && (@.. broadcast = false u = uᵢ₋₁)
    # for the second to the ms[ccache.mdeg] th stages
    for i in 2:(ccache.mdeg)
        μ, κ = recf[ccache.start + (i - 2) * 2 + 1], recf[ccache.start + (i - 2) * 2 + 2]
        ν = -1 - κ
        f(k, uᵢ₋₁, p, tᵢ₋₁)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tᵢ₋₁ = dt * μ - ν * tᵢ₋₂ - κ * tᵢ₋₃
        @.. broadcast = false u = (dt * μ) * k - ν * uᵢ₋₁ - κ * uᵢ₋₂
        if i < ccache.mdeg
            @.. broadcast = false uᵢ₋₂ = uᵢ₋₁
            @.. broadcast = false uᵢ₋₁ = u
        end
        tᵢ₋₃ = tᵢ₋₂
        tᵢ₋₂ = tᵢ₋₁
    end # end if
    # two-stage finishing procedure.
    δt₁ = dt * fp1[ccache.deg_index]
    δt₂ = dt * fp2[ccache.deg_index]
    f(k, u, p, tᵢ₋₁)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false uᵢ₋₁ = u + δt₁ * k
    if integrator.opts.adaptive
        @.. broadcast = false tmp = -δt₂ * k
    else
        @.. broadcast = false u = -δt₂ * k
    end
    c = DiffEqBase.value(sign(δt₁)) * integrator.opts.internalnorm(δt₁, t)
    tᵢ₋₁ += c
    f(k, uᵢ₋₁, p, tᵢ₋₁)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        @.. broadcast = false tmp += δt₂ * k
        @.. broadcast = false u = uᵢ₋₁ + δt₁ * k + tmp
    else
        @.. broadcast = false u += uᵢ₋₁ + (δt₁ + δt₂) * k
    end

    # error estimate
    if integrator.opts.adaptive
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::ROCK4ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    alg = unwrap_alg(integrator, true)
    cache.max_stage = (alg.max_stages < 1 || alg.max_stages > 152) ? 152 : alg.max_stages
    cache.min_stage = (alg.min_stages > cache.max_stage) ? cache.max_stage : alg.min_stages
    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ROCK4ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; ms, fpa, fpb, fpbe, recf) = cache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    # The the number of degree for Chebyshev polynomial
    mdeg = floor(Int, sqrt((3 + abs(dt) * integrator.eigen_est) / 0.353)) + 1
    mdeg = min(max(mdeg, cache.min_stage), cache.max_stage)
    cache.mdeg = max(mdeg, 5) - 4
    choosedeg!(cache)
    # recurrence
    # for the first stage
    tᵢ₋₁ = t + dt * recf[cache.start]
    tᵢ₋₂ = t + dt * recf[cache.start]
    tᵢ₋₃ = t
    uᵢ₋₂ = copy(uprev)
    uᵢ₋₁ = uprev + (dt * recf[cache.start]) * fsalfirst
    if cache.mdeg < 2
        u = uᵢ₋₁
    end
    # for the second to the cache.mdeg th stages
    for i in 2:(cache.mdeg)
        μ, κ = recf[cache.start + (i - 2) * 2 + 1], recf[cache.start + (i - 2) * 2 + 2]
        ν = -1 - κ
        u = f(uᵢ₋₁, p, tᵢ₋₁)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tᵢ₋₁ = dt * μ - ν * tᵢ₋₂ - κ * tᵢ₋₃
        u = (dt * μ) * u - ν * uᵢ₋₁ - κ * uᵢ₋₂
        if i < cache.mdeg
            uᵢ₋₂ = uᵢ₋₁
            uᵢ₋₁ = u
        end
        tᵢ₋₃ = tᵢ₋₂
        tᵢ₋₂ = tᵢ₋₁
    end

    # These constants correspond to the Buther Tableau coefficients of explicit RK methods
    a₂₁ = dt * fpa[cache.deg_index][1]
    a₃₁ = dt * fpa[cache.deg_index][2]
    a₃₂ = dt * fpa[cache.deg_index][3]
    a₄₁ = dt * fpa[cache.deg_index][4]
    a₄₂ = dt * fpa[cache.deg_index][5]
    a₄₃ = dt * fpa[cache.deg_index][6]
    B₁ = dt * fpb[cache.deg_index][1]
    B₂ = dt * fpb[cache.deg_index][2]
    B₃ = dt * fpb[cache.deg_index][3]
    B₄ = dt * fpb[cache.deg_index][4]
    # coefficients of embedded method for error estimation
    B̂₁ = dt * (fpbe[cache.deg_index][1] - fpb[cache.deg_index][1])
    B̂₂ = dt * (fpbe[cache.deg_index][2] - fpb[cache.deg_index][2])
    B̂₃ = dt * (fpbe[cache.deg_index][3] - fpb[cache.deg_index][3])
    B̂₄ = dt * (fpbe[cache.deg_index][4] - fpb[cache.deg_index][4])
    B̂₅ = dt * fpbe[cache.deg_index][5]

    # 4-stage finishing procedure.
    # Stage-1
    uᵢ₋₁ = f(u, p, tᵢ₋₁)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    uᵢ₋₂ = u + a₃₁ * uᵢ₋₁
    uᵢ₋₃ = u + a₄₁ * uᵢ₋₁
    u += B₁ * uᵢ₋₁
    tmp = u  # Initialize for JET
    if integrator.opts.adaptive
        tmp = B̂₁ * uᵢ₋₁
    end
    uᵢ₋₁ = u + (a₂₁ - B₁) * uᵢ₋₁

    # Stage-2
    c₂ = a₂₁
    _c₂ = DiffEqBase.value(sign(c₂)) * integrator.opts.internalnorm(c₂, t)
    tᵢ₋₂ = tᵢ₋₁ + _c₂
    uᵢ₋₁ = f(uᵢ₋₁, p, tᵢ₋₂)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    uᵢ₋₂ += a₃₂ * uᵢ₋₁
    uᵢ₋₃ += a₄₂ * uᵢ₋₁
    u += B₂ * uᵢ₋₁
    if integrator.opts.adaptive
        tmp += B̂₂ * uᵢ₋₁
    end

    # Stage-3
    c₃ = a₃₁ + a₃₂
    _c₃ = DiffEqBase.value(sign(c₃)) * integrator.opts.internalnorm(c₃, t)
    tᵢ₋₂ = tᵢ₋₁ + _c₃
    uᵢ₋₂ = f(uᵢ₋₂, p, tᵢ₋₂)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    uᵢ₋₃ += a₄₃ * uᵢ₋₂
    u += B₃ * uᵢ₋₂
    if integrator.opts.adaptive
        tmp += B̂₃ * uᵢ₋₂
    end

    #Stage-4
    c₄ = a₄₁ + a₄₂ + a₄₃
    _c₄ = DiffEqBase.value(sign(c₄)) * integrator.opts.internalnorm(c₄, t)
    tᵢ₋₂ = tᵢ₋₁ + _c₄
    uᵢ₋₃ = f(uᵢ₋₃, p, tᵢ₋₂)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    u += B₄ * uᵢ₋₃
    if integrator.opts.adaptive
        tmp += B̂₄ * uᵢ₋₃
    end

    uᵢ₋₁ = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    #Error estimate (embedded method of order 3)
    if integrator.opts.adaptive
        tmp += B̂₅ * uᵢ₋₁
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast = uᵢ₋₁
    integrator.u = u
end

function initialize!(integrator, cache::ROCK4Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)

    alg = unwrap_alg(integrator, true)
    cache.constantcache.max_stage = (alg.max_stages < 1 || alg.max_stages > 152) ? 152 :
        alg.max_stages
    cache.constantcache.min_stage = (alg.min_stages > cache.constantcache.max_stage) ?
        cache.constantcache.max_stage : alg.min_stages

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::ROCK4Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; uᵢ₋₁, uᵢ₋₂, uᵢ₋₃, tmp, atmp, k) = cache
    (; ms, fpa, fpb, fpbe, recf) = cache.constantcache
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    # The the number of degree for Chebyshev polynomial
    mdeg = floor(Int, sqrt((3 + abs(dt) * integrator.eigen_est) / 0.353)) + 1
    mdeg = min(max(mdeg, ccache.min_stage), ccache.max_stage)
    ccache.mdeg = max(mdeg, 5) - 4
    choosedeg!(cache)
    # recurrence
    # for the first stage
    tᵢ₋₁ = t + dt * recf[ccache.start]
    tᵢ₋₂ = t + dt * recf[ccache.start]
    tᵢ₋₃ = t
    @.. broadcast = false uᵢ₋₂ = uprev
    @.. broadcast = false uᵢ₋₁ = uprev + (dt * recf[ccache.start]) * fsalfirst
    if ccache.mdeg < 2
        @.. broadcast = false u = uᵢ₋₁
    end
    # for the second to the ccache.mdeg th stages
    for i in 2:(ccache.mdeg)
        μ, κ = recf[ccache.start + (i - 2) * 2 + 1], recf[ccache.start + (i - 2) * 2 + 2]
        ν = -1 - κ
        f(k, uᵢ₋₁, p, tᵢ₋₁)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tᵢ₋₁ = (dt * μ) - ν * tᵢ₋₂ - κ * tᵢ₋₃
        @.. broadcast = false u = (dt * μ) * k - ν * uᵢ₋₁ - κ * uᵢ₋₂
        if i < ccache.mdeg
            @.. broadcast = false uᵢ₋₂ = uᵢ₋₁
            @.. broadcast = false uᵢ₋₁ = u
        end
        tᵢ₋₃ = tᵢ₋₂
        tᵢ₋₂ = tᵢ₋₁
    end

    # These constants correspond to the Buther Tableau coefficients of explicit RK methods
    a₂₁ = dt * fpa[ccache.deg_index][1]
    a₃₁ = dt * fpa[ccache.deg_index][2]
    a₃₂ = dt * fpa[ccache.deg_index][3]
    a₄₁ = dt * fpa[ccache.deg_index][4]
    a₄₂ = dt * fpa[ccache.deg_index][5]
    a₄₃ = dt * fpa[ccache.deg_index][6]
    B₁ = dt * fpb[ccache.deg_index][1]
    B₂ = dt * fpb[ccache.deg_index][2]
    B₃ = dt * fpb[ccache.deg_index][3]
    B₄ = dt * fpb[ccache.deg_index][4]
    # coefficients of embedded method for error estimation
    B̂₁ = dt * (fpbe[ccache.deg_index][1] - fpb[ccache.deg_index][1])
    B̂₂ = dt * (fpbe[ccache.deg_index][2] - fpb[ccache.deg_index][2])
    B̂₃ = dt * (fpbe[ccache.deg_index][3] - fpb[ccache.deg_index][3])
    B̂₄ = dt * (fpbe[ccache.deg_index][4] - fpb[ccache.deg_index][4])
    B̂₅ = dt * fpbe[ccache.deg_index][5]

    # 4-stage finishing procedure.
    # Stage-1
    f(k, u, p, tᵢ₋₁)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false uᵢ₋₂ = u + a₃₁ * k
    @.. broadcast = false uᵢ₋₃ = u + a₄₁ * k
    @.. broadcast = false uᵢ₋₁ = u + a₂₁ * k
    @.. broadcast = false u += B₁ * k
    if integrator.opts.adaptive
        @.. broadcast = false tmp = B̂₁ * k
    end

    # Stage-2
    c₂ = a₂₁
    _c₂ = value(sign(c₂)) * integrator.opts.internalnorm(c₂, t)
    tᵢ₋₂ = tᵢ₋₁ + _c₂
    f(k, uᵢ₋₁, p, tᵢ₋₂)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false uᵢ₋₂ += a₃₂ * k
    @.. broadcast = false uᵢ₋₃ += a₄₂ * k
    @.. broadcast = false u += B₂ * k
    if integrator.opts.adaptive
        @.. broadcast = false tmp += B̂₂ * k
    end

    # Stage-3
    c₃ = a₃₁ + a₃₂
    _c₃ = DiffEqBase.value(sign(c₃)) * integrator.opts.internalnorm(c₃, t)
    tᵢ₋₂ = tᵢ₋₁ + _c₃
    f(k, uᵢ₋₂, p, tᵢ₋₂)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false uᵢ₋₃ += a₄₃ * k
    @.. broadcast = false u += B₃ * k
    if integrator.opts.adaptive
        @.. broadcast = false tmp += B̂₃ * k
    end

    #Stage-4
    c₄ = a₄₁ + a₄₂ + a₄₃
    _c₄ = DiffEqBase.value(sign(c₄)) * integrator.opts.internalnorm(c₄, t)
    tᵢ₋₂ = tᵢ₋₁ + _c₄
    f(k, uᵢ₋₃, p, tᵢ₋₂)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false u += B₄ * k
    if integrator.opts.adaptive
        @.. broadcast = false tmp += B̂₄ * k
    end

    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    #Error estimate (embedded method of order 3)
    if integrator.opts.adaptive
        @.. broadcast = false tmp += B̂₅ * k
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    @.. broadcast = false integrator.fsallast = k
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::RKCConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::RKCConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    # The number of degree for Chebyshev polynomial
    mdeg = floor(Int, sqrt(1.54 * abs(dt) * integrator.eigen_est + 1)) + 1

    T = typeof(t)

    w0 = one(T) + T(2) / (13 * (mdeg^2))
    temp1 = w0^2 - one(T)
    temp2 = sqrt(temp1)
    arg = mdeg * log(w0 + temp2)
    w1 = (sinh(arg) * temp1) / (cosh(arg) * mdeg * temp2 - w0 * sinh(arg))
    b1 = one(T) / ((T(2) * w0)^2)
    b2 = b1

    # stage-1
    gprev2 = copy(uprev)
    μs = w1 * b1
    gprev = uprev + dt * μs * fsalfirst
    th2 = zero(eltype(u))
    th1 = μs
    z1 = w0
    z2 = one(eltype(u))
    dz1 = one(eltype(u))
    dz2 = zero(eltype(u))
    d2z1 = zero(eltype(u))
    d2z2 = zero(eltype(u))

    # stage 2 - mdeg
    for iter in 2:mdeg
        z = 2 * w0 * z1 - z2
        dz = 2 * w0 * dz1 - dz2 + 2 * z1
        d2z = 2 * w0 * d2z1 - d2z2 + 4 * dz1
        b = d2z / (dz^2)
        νs = 1 - z1 * b1
        μ = (2 * w0 * b) / b1
        ν = -b / b2
        μs = μ * w1 / w0
        #using u as temporary storage
        u = f(gprev, p, t + dt * th1)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        u = μ * gprev + ν * gprev2 + (1 - μ - ν) * uprev + dt * μs * (u - νs * fsalfirst)
        th = μ * th1 + ν * th2 + μs * (1 - νs)
        if (iter < mdeg)
            gprev2 = gprev
            gprev = u
            th2 = th1
            th1 = th
            b2 = b1
            b1 = b
            z2 = z1
            z1 = z
            dz2 = dz1
            dz1 = dz
            d2z2 = d2z1
            d2z1 = d2z
        end
    end
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # error estimate
    if integrator.opts.adaptive
        tmp = 0.8 * (uprev - u) + 0.4 * dt * (fsalfirst + integrator.fsallast)
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::RKCCache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::RKCCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; k, tmp, gprev2, gprev, atmp) = cache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    # The number of degree for Chebyshev polynomial
    mdeg = floor(Int, sqrt(1.54 * abs(dt) * integrator.eigen_est + 1)) + 1

    T = typeof(t)

    w0 = one(T) + T(2) / (13 * (mdeg^2))
    temp1 = w0^2 - one(T)
    temp2 = sqrt(temp1)
    arg = mdeg * log(w0 + temp2)
    w1 = (sinh(arg) * temp1) / (cosh(arg) * mdeg * temp2 - w0 * sinh(arg))
    b1 = one(T) / ((T(2) * w0)^2)
    b2 = b1

    # stage-1
    @.. broadcast = false gprev2 = uprev
    μs = w1 * b1
    @.. broadcast = false gprev = uprev + dt * μs * fsalfirst
    th2 = zero(eltype(u))
    th1 = μs
    z1 = w0
    z2 = one(eltype(u))
    dz1 = one(eltype(u))
    dz2 = zero(eltype(u))
    d2z1 = zero(eltype(u))
    d2z2 = zero(eltype(u))

    # stage 2 - mdeg
    for iter in 2:mdeg
        z = 2 * w0 * z1 - z2
        dz = 2 * w0 * dz1 - dz2 + 2 * z1
        d2z = 2 * w0 * d2z1 - d2z2 + 4 * dz1
        b = d2z / (dz^2)
        νs = 1 - z1 * b1
        μ = (2 * w0 * b) / b1
        ν = -b / b2
        μs = μ * w1 / w0
        f(k, gprev, p, t + dt * th1)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false u = μ * gprev + ν * gprev2 + (1 - μ - ν) * uprev +
            dt * μs * (k - νs * fsalfirst)
        th = μ * th1 + ν * th2 + μs * (1 - νs)
        if (iter < mdeg)
            @.. broadcast = false gprev2 = gprev
            @.. broadcast = false gprev = u
            th2 = th1
            th1 = th
            b2 = b1
            b1 = b
            z2 = z1
            z1 = z
            dz2 = dz1
            dz1 = dz
            d2z2 = d2z1
            d2z1 = d2z
        end
    end
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # error estimate
    if integrator.opts.adaptive
        @.. broadcast = false tmp = 0.8 * (uprev - u) + 0.4 * dt * (fsalfirst + integrator.fsallast)
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::ESERK4ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ESERK4ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; ms, Cᵤ, Cₑ) = cache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    mdeg = Int(floor(sqrt(abs(dt) * integrator.eigen_est)) + 1)
    mdeg = (mdeg > 4000) ? 4000 : mdeg
    cache.mdeg = mdeg
    choosedeg_SERK!(integrator, cache)
    mdeg = cache.mdeg
    start = cache.start
    internal_deg = cache.internal_deg
    α = 2.0 / (mdeg^2)

    u = zero(uprev)
    tmp = zero(uprev)

    for i in 1:4
        hᵢ = dt / i
        tᵢ = t
        Sᵢ = zero(u)
        uᵢ₋₁ = uprev
        uᵢ₋₂ = zero(u)
        for j in 1:i
            r = tᵢ
            Sᵢ = (cache.Bᵢ[start]) * uᵢ₋₁
            for st in 1:mdeg
                k = f(uᵢ₋₁, p, r)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

                if st % internal_deg == 1
                    uᵢ = uᵢ₋₁ + α * hᵢ * k
                else
                    uᵢ = 2 * uᵢ₋₁ - uᵢ₋₂ + 2 * α * hᵢ * k
                end
                q = convert(Int, floor(st / internal_deg))
                r = tᵢ + α * (st^2 + q * internal_deg^2) * hᵢ
                Sᵢ = Sᵢ + (cache.Bᵢ[start + st]) * uᵢ
                if st < mdeg
                    uᵢ₋₂ = uᵢ₋₁
                    uᵢ₋₁ = uᵢ
                end
            end

            if j < i
                tᵢ = tᵢ + hᵢ
                uᵢ₋₁ = Sᵢ
            end
        end

        u = u + Cᵤ[i] * Sᵢ
        integrator.opts.adaptive && (tmp = tmp + Cₑ[i] * Sᵢ)
    end

    u = u / 6
    if integrator.opts.adaptive
        tmp = tmp / 6
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

function initialize!(integrator, cache::ESERK4Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::ESERK4Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; uᵢ, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, k) = cache
    (; ms, Cᵤ, Cₑ) = cache.constantcache
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    mdeg = Int(floor(sqrt(abs(dt) * integrator.eigen_est)) + 1)
    mdeg = (mdeg > 4000) ? 4000 : mdeg
    ccache.mdeg = mdeg
    choosedeg_SERK!(integrator, cache)
    mdeg = ccache.mdeg
    start = ccache.start
    internal_deg = ccache.internal_deg
    α = 2.0 / (mdeg^2)

    @.. broadcast = false u = zero(uprev)
    @.. broadcast = false tmp = zero(uprev)
    for i in 1:4
        hᵢ = dt / i
        tᵢ = t
        @.. broadcast = false Sᵢ = zero(u)
        @.. broadcast = false uᵢ₋₁ = uprev
        @.. broadcast = false uᵢ₋₂ = zero(u)
        for j in 1:i
            r = tᵢ
            @.. broadcast = false Sᵢ = (cache.constantcache.Bᵢ[start]) * uᵢ₋₁
            for st in 1:mdeg
                f(k, uᵢ₋₁, p, r)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

                if st % internal_deg == 1
                    @.. broadcast = false uᵢ = uᵢ₋₁ + α * hᵢ * k
                else
                    @.. broadcast = false uᵢ = 2 * uᵢ₋₁ - uᵢ₋₂ + 2 * α * hᵢ * k
                end
                q = convert(Int, floor(st / internal_deg))
                r = tᵢ + α * (st^2 + q * internal_deg^2) * hᵢ
                @.. broadcast = false Sᵢ = Sᵢ + (cache.constantcache.Bᵢ[start + st]) * uᵢ
                if st < mdeg
                    @.. broadcast = false uᵢ₋₂ = uᵢ₋₁
                    @.. broadcast = false uᵢ₋₁ = uᵢ
                end
            end

            if j < i
                tᵢ = tᵢ + hᵢ
                @.. broadcast = false uᵢ₋₁ = Sᵢ
            end
        end

        @.. broadcast = false u = u + Cᵤ[i] * Sᵢ
        integrator.opts.adaptive && (@.. broadcast = false tmp = tmp + Cₑ[i] * Sᵢ)
    end

    @.. broadcast = false u = u / 6

    if integrator.opts.adaptive
        @.. broadcast = false tmp = tmp / 6
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::ESERK5ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ESERK5ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; ms, Cᵤ, Cₑ, Bᵢ) = cache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    mdeg = Int(floor(sqrt(abs(dt) * integrator.eigen_est / 0.98)) + 1)
    mdeg = (mdeg > 2000) ? 2000 : mdeg
    cache.mdeg = mdeg
    choosedeg_SERK!(integrator, cache)
    mdeg = cache.mdeg
    start = cache.start
    internal_deg = cache.internal_deg
    α = 100.0 / (49.0 * mdeg^2)

    u = zero(uprev)
    tmp = zero(uprev)
    for i in 1:5
        hᵢ = dt / i
        tᵢ = t
        Sᵢ = zero(u)
        uᵢ₋₁ = uprev
        uᵢ₋₂ = zero(u)
        for j in 1:i
            r = tᵢ
            Sᵢ = (Bᵢ[start]) * uᵢ₋₁
            for st in 1:mdeg
                k = f(uᵢ₋₁, p, r)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

                if st % internal_deg == 1
                    uᵢ = uᵢ₋₁ + α * hᵢ * k
                else
                    uᵢ = 2 * uᵢ₋₁ - uᵢ₋₂ + 2 * α * hᵢ * k
                end
                q = convert(Int, floor(st / internal_deg))
                r = tᵢ + α * (st^2 + q * internal_deg^2) * hᵢ
                Sᵢ = Sᵢ + (Bᵢ[start + st]) * uᵢ
                if st < mdeg
                    uᵢ₋₂ = uᵢ₋₁
                    uᵢ₋₁ = uᵢ
                end
            end

            if j < i
                tᵢ = tᵢ + hᵢ
                uᵢ₋₁ = Sᵢ
            end
        end

        u = u + Cᵤ[i] * Sᵢ
        integrator.opts.adaptive && (tmp = tmp + Cₑ[i] * Sᵢ)
    end

    u = u / 24
    if integrator.opts.adaptive
        tmp = tmp / 24
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

function initialize!(integrator, cache::ESERK5Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::ESERK5Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; uᵢ, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, k) = cache
    (; ms, Cᵤ, Cₑ, Bᵢ) = cache.constantcache
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    mdeg = Int(floor(sqrt(abs(dt) * integrator.eigen_est / 0.98)) + 1)
    mdeg = (mdeg > 2000) ? 2000 : mdeg
    ccache.mdeg = mdeg
    choosedeg_SERK!(integrator, cache)
    mdeg = ccache.mdeg
    start = ccache.start
    internal_deg = ccache.internal_deg
    α = 100.0 / (49.0 * mdeg^2)

    @.. broadcast = false u = zero(uprev)
    @.. broadcast = false tmp = zero(uprev)
    for i in 1:5
        hᵢ = dt / i
        tᵢ = t
        @.. broadcast = false Sᵢ = zero(u)
        @.. broadcast = false uᵢ₋₁ = uprev
        @.. broadcast = false uᵢ₋₂ = zero(u)
        for j in 1:i
            r = tᵢ
            @.. broadcast = false Sᵢ = (Bᵢ[start]) * uᵢ₋₁
            for st in 1:mdeg
                f(k, uᵢ₋₁, p, r)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

                if st % internal_deg == 1
                    @.. broadcast = false uᵢ = uᵢ₋₁ + α * hᵢ * k
                else
                    @.. broadcast = false uᵢ = 2 * uᵢ₋₁ - uᵢ₋₂ + 2 * α * hᵢ * k
                end
                q = convert(Int, floor(st / internal_deg))
                r = tᵢ + α * (st^2 + q * internal_deg^2) * hᵢ
                @.. broadcast = false Sᵢ = Sᵢ + (Bᵢ[start + st]) * uᵢ
                if st < mdeg
                    @.. broadcast = false uᵢ₋₂ = uᵢ₋₁
                    @.. broadcast = false uᵢ₋₁ = uᵢ
                end
            end

            if j < i
                tᵢ = tᵢ + hᵢ
                @.. broadcast = false uᵢ₋₁ = Sᵢ
            end
        end

        @.. broadcast = false u = u + Cᵤ[i] * Sᵢ
        integrator.opts.adaptive && (@.. broadcast = false tmp = tmp + Cₑ[i] * Sᵢ)
    end

    @.. broadcast = false u = u / 24

    if integrator.opts.adaptive
        @.. broadcast = false tmp = tmp / 24
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::SERK2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::SERK2ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; ms, Bᵢ) = cache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    mdeg = Int(floor(sqrt(abs(dt) * integrator.eigen_est / 0.8)) + 1)
    mdeg = (mdeg > 250) ? 250 : mdeg
    cache.mdeg = mdeg
    choosedeg_SERK!(integrator, cache)
    mdeg = cache.mdeg
    start = cache.start
    internal_deg = cache.internal_deg
    α = 1.0 / (0.4 * mdeg^2)

    uᵢ₋₁ = uprev
    uᵢ₋₂ = uprev
    Sᵢ = Bᵢ[start] * uprev
    for i in 1:10
        k = f(uᵢ₋₁, p, t + (1 + (i - 1) * internal_deg^2) * α * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        u = uᵢ₋₁ + α * dt * k
        Sᵢ = Sᵢ + Bᵢ[start + (i - 1) * internal_deg + 1] * u
        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = u
        for j in 2:internal_deg
            k = f(uᵢ₋₁, p, t + (j^2 + (i - 1) * internal_deg^2) * α * dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            u = 2 * uᵢ₋₁ - uᵢ₋₂ + 2 * α * dt * k
            Sᵢ = Sᵢ + Bᵢ[start + j + (i - 1) * internal_deg] * u
            if j * i < mdeg
                uᵢ₋₂ = uᵢ₋₁
                uᵢ₋₁ = u
            end
        end
    end
    u = Sᵢ
    k = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        tmp = u - uprev - dt * k
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast = k
    integrator.u = u
end

function initialize!(integrator, cache::SERK2Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::SERK2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, k) = cache
    (; ms, Bᵢ) = cache.constantcache
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    mdeg = Int(floor(sqrt(abs(dt) * integrator.eigen_est / 0.8)) + 1)
    mdeg = (mdeg > 250) ? 250 : mdeg
    ccache.mdeg = mdeg
    choosedeg_SERK!(integrator, cache)
    mdeg = ccache.mdeg
    start = ccache.start
    internal_deg = ccache.internal_deg
    α = 1.0 / (0.4 * mdeg^2)

    @.. broadcast = false uᵢ₋₁ = uprev
    @.. broadcast = false uᵢ₋₂ = uprev
    @.. broadcast = false Sᵢ = Bᵢ[start] * uprev
    for i in 1:10
        f(k, uᵢ₋₁, p, t + (1 + (i - 1) * internal_deg^2) * α * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false u = uᵢ₋₁ + α * dt * k
        @.. broadcast = false Sᵢ = Sᵢ + Bᵢ[start + (i - 1) * internal_deg + 1] * u
        @.. broadcast = false uᵢ₋₂ = uᵢ₋₁
        @.. broadcast = false uᵢ₋₁ = u
        for j in 2:internal_deg
            f(k, uᵢ₋₂, p, t + (j^2 + (i - 1) * internal_deg^2) * α * dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            @.. broadcast = false u = 2 * uᵢ₋₁ - uᵢ₋₂ + 2 * α * dt * k
            @.. broadcast = false Sᵢ = Sᵢ + Bᵢ[start + j + (i - 1) * internal_deg] * u
            if j < mdeg
                @.. broadcast = false uᵢ₋₂ = uᵢ₋₁
                @.. broadcast = false uᵢ₋₁ = u
            end
        end
    end
    @.. broadcast = false u = Sᵢ
    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        @.. broadcast = false tmp = u - uprev - dt * k
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast = k
    integrator.u = u
end


function initialize!(integrator, cache::TSRKC2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::TSRKC2ConstantCache, repeat_step = false)
    (; t, tprev, dt, uprev, u, f, p, fsalfirst, uprev2) = integrator

    q = (t - tprev) / dt
    onemq = 1 - q

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    dtsw0 = zero(typeof(cache.tsw0))
    d2tsw0 = zero(typeof(cache.tsw0))

    mdeg = floor(Int, sqrt(1 + 0.759782816506459 * abs(dt) * integrator.eigen_est * (onemq + sqrt(1 + q * (q - 0.598626091572911))))) + 1
    tsw0 = cache.tsw0
    acoshtsw0 = cache.acoshtsw0
    acoshtsw0dm = acoshtsw0 / mdeg
    w0 = cosh(acoshtsw0dm)
    w0sqm1 = w0 * w0 - 1
    dtsw0 = mdeg * sinh(acoshtsw0) / sqrt(w0sqm1)
    d2tsw0 = (mdeg * mdeg * tsw0 - w0 * dtsw0) / w0sqm1
    w1 = (onemq * dtsw0 + sqrt(onemq * onemq * dtsw0 * dtsw0 + 4 * q * tsw0 * d2tsw0)) / (2 * d2tsw0)

    # stage-1
    gprev2 = uprev
    μs = w1 / w0
    gprev = uprev + dt * μs * fsalfirst
    th2 = zero(eltype(u))
    th1 = μs
    z1 = w0
    z2 = one(eltype(u))
    z = 2 * w0 * z1 - z2

    # stage 2 - mdeg
    for iter in 2:mdeg
        μ = 2 * w0 * z1 / z
        ν = -z2 / z
        μs = μ * w1 / w0
        #using u as temporary storage
        u = f(gprev, p, t + dt * th1)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        u = μ * gprev + ν * gprev2 + dt * μs * u
        if (iter < mdeg)
            gprev2 = gprev
            gprev = u
            th = μ * th1 + ν * th2 + μs
            th2 = th1
            th1 = th
            z2 = z1
            z1 = z
            z = 2 * w0 * z1 - z2
        end
    end

    g = (1 + q) * tsw0 / (q * tsw0 + w1 * dtsw0)
    μ = 1 - g
    u = μ * uprev2 + g * u

    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # error estimate
    if integrator.opts.adaptive
        tmp = 0.3333333333333333 * (uprev - u + dt * integrator.fsallast)
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::TSRKC2Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::TSRKC2Cache, repeat_step = false)
    (; t, tprev, dt, uprev, u, f, p, fsalfirst, uprev2) = integrator
    (; k, tmp, gprev2, gprev, atmp, constantcache) = cache

    q = (t - tprev) / dt
    onemq = 1 - q

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    dtsw0 = zero(typeof(constantcache.tsw0))
    d2tsw0 = zero(typeof(constantcache.tsw0))

    mdeg = floor(Int, sqrt(1 + 0.759782816506459 * abs(dt) * integrator.eigen_est * (onemq + sqrt(1 + q * (q - 0.598626091572911))))) + 1
    tsw0 = constantcache.tsw0
    acoshtsw0 = constantcache.acoshtsw0
    acoshtsw0dm = acoshtsw0 / mdeg
    w0 = cosh(acoshtsw0dm)
    w0sqm1 = w0 * w0 - 1
    dtsw0 = mdeg * sinh(acoshtsw0) / sqrt(w0sqm1)
    d2tsw0 = (mdeg * mdeg * tsw0 - w0 * dtsw0) / w0sqm1
    w1 = (onemq * dtsw0 + sqrt(onemq * onemq * dtsw0 * dtsw0 + 4 * q * tsw0 * d2tsw0)) / (2 * d2tsw0)

    # stage-1
    @.. broadcast = false gprev2 = uprev
    μs = w1 / w0
    @.. broadcast = false gprev = uprev + dt * μs * fsalfirst
    th2 = zero(eltype(u))
    th1 = μs
    z1 = w0
    z2 = one(eltype(u))
    z = 2 * w0 * z1 - z2

    # stage 2 - mdeg
    for iter in 2:mdeg
        μ = 2 * w0 * z1 / z
        ν = -z2 / z
        μs = μ * w1 / w0
        f(k, gprev, p, t + dt * th1)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false u = μ * gprev + ν * gprev2 + dt * μs * k
        if (iter < mdeg)
            @.. broadcast = false gprev2 = gprev
            @.. broadcast = false gprev = u
            th = μ * th1 + ν * th2 + μs
            th2 = th1
            th1 = th
            z2 = z1
            z1 = z
            z = 2 * w0 * z1 - z2
        end
    end

    g = (1 + q) * tsw0 / (q * tsw0 + w1 * dtsw0)
    μ = 1 - g
    @.. broadcast = false u = μ * uprev2 + g * u

    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # error estimate
    if integrator.opts.adaptive
        @.. broadcast = false tmp = 0.3333333333333333 * (uprev - u + dt * integrator.fsallast)
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end


function initialize!(integrator, cache::TSRKC3ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::TSRKC3ConstantCache, repeat_step = false)
    (; t, tprev, dt, uprev, u, f, p, fsalfirst, uprev2) = integrator

    q = (t - tprev) / dt
    onemq = 1 - q
    onepq = 1 + q
    onepq2 = onepq * onepq

    # The first and possibly second steps are calculated via the one-step RKC method.
    rkcstep = q < 0.49

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    dtsw0 = zero(typeof(cache.tsw0))
    d2tsw0 = zero(typeof(cache.tsw0))
    if (rkcstep)
        mdeg = floor(Int, sqrt(1.54 * abs(dt) * integrator.eigen_est + 1)) + 1
        w0 = 1 + 2 / (13 * (mdeg^2))
        temp1 = w0^2 - 1
        temp2 = sqrt(temp1)
        arg = mdeg * log(w0 + temp2)
        w1 = (sinh(arg) * temp1) / (cosh(arg) * mdeg * temp2 - w0 * sinh(arg))
        b1 = 1 / ((2 * w0)^2)
        b = b1
        b2 = b1
    else
        mdeg = floor(Int, sqrt(4 + 1.267029788142009 * abs(dt) * integrator.eigen_est * (onemq + sqrt(1 + q * (0.44256220745562963 + q))))) + 1
        mdeg2 = mdeg^2
        tsw0 = cache.tsw0
        acoshtsw0 = cache.acoshtsw0
        acoshtsw0dm = acoshtsw0 / mdeg
        w0 = cosh(acoshtsw0dm)
        w0sq = w0^2
        w0sqm1 = w0sq - 1
        dtsw0 = mdeg * sinh(acoshtsw0) / sqrt(w0sqm1)
        d2tsw0 = (mdeg2 * tsw0 - w0 * dtsw0) / w0sqm1
        d3tsw0 = ((1 + 2 * w0sq + mdeg2 * w0sqm1) * dtsw0 - 3 * mdeg2 * w0 * tsw0) / (w0sqm1^2)
        w1 = (onemq * d2tsw0 + sqrt((onemq * d2tsw0)^2 + 4 * q * dtsw0 * d3tsw0)) / (2 * d3tsw0)

        b1 = 0.25 * sinh((mdeg - 2) * acoshtsw0dm) / sinh((mdeg - 1) * acoshtsw0dm)
        b = 15 / ((8 * w0)^2)
        b2 = b1
    end

    # stage-1
    gprev2 = uprev
    μs = w1 * b1
    gprev = uprev + dt * μs * fsalfirst
    th2 = zero(eltype(u))
    th1 = μs
    z1 = w0
    z2 = one(eltype(u))
    dz1 = one(eltype(u))
    dz2 = zero(eltype(u))
    d2z1 = zero(eltype(u))
    d2z2 = zero(eltype(u))
    z = 2 * w0 * z1 - z2
    dz = 2 * w0 * dz1 - dz2 + 2 * z1
    d2z = 2 * w0 * d2z1 - d2z2 + 4 * dz1

    # stage 2 - mdeg
    for iter in 2:mdeg
        νs = 1 - z1 * b1
        μ = (2 * w0 * b) / b1
        ν = -b / b2
        μs = μ * w1 / w0
        #using u as temporary storage
        u = f(gprev, p, t + dt * th1)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        u = μ * gprev + ν * gprev2 + (1 - μ - ν) * uprev + dt * μs * (u - νs * fsalfirst)
        if (iter < mdeg)
            gprev2 = gprev
            gprev = u
            th = μ * th1 + ν * th2 + μs * (1 - νs)
            th2 = th1
            th1 = th
            b2 = b1
            b1 = b
            z2 = z1
            z1 = z
            dz2 = dz1
            dz1 = dz
            d2z2 = d2z1
            d2z1 = d2z

            z = 2 * w0 * z1 - z2
            dz = 2 * w0 * dz1 - dz2 + 2 * z1
            d2z = 2 * w0 * d2z1 - d2z2 + 4 * dz1

            b = d2z / (dz^2)
        end
    end

    if (!rkcstep)
        a = onepq / (q * dtsw0 + w1 * d2tsw0)
        g = (a * dtsw0 - 1) / q
        b = a / (b * w1)
        μ = 1 - g - b
        u = μ * uprev + g * uprev2 + b * u
    end

    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # error estimate
    if integrator.opts.adaptive
        if (rkcstep)
            tmp = 1.2 * (uprev - u) + 0.6 * dt * (fsalfirst + integrator.fsallast)
        else
            tmp = 0.6 * (uprev / q - uprev2 / (q * onepq2) - u * (2 + q) / onepq2 + dt * integrator.fsallast / onepq)
        end
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::TSRKC3Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::TSRKC3Cache, repeat_step = false)
    (; t, tprev, dt, uprev, u, f, p, fsalfirst, uprev2) = integrator
    (; k, tmp, gprev2, gprev, atmp, constantcache) = cache

    q = (t - tprev) / dt
    onemq = 1 - q
    onepq = 1 + q
    onepq2 = onepq * onepq

    # The first and possibly second steps are calculated via the one-step RKC method.
    rkcstep = q < 0.49

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    dtsw0 = zero(typeof(constantcache.tsw0))
    d2tsw0 = zero(typeof(constantcache.tsw0))
    if (rkcstep)
        mdeg = floor(Int, sqrt(1.54 * abs(dt) * integrator.eigen_est + 1)) + 1
        w0 = 1 + 2 / (13 * (mdeg^2))
        temp1 = w0^2 - 1
        temp2 = sqrt(temp1)
        arg = mdeg * log(w0 + temp2)
        w1 = (sinh(arg) * temp1) / (cosh(arg) * mdeg * temp2 - w0 * sinh(arg))
        b1 = 1 / ((2 * w0)^2)
        b = b1
        b2 = b1
    else
        mdeg = floor(Int, sqrt(4 + 1.267029788142009 * abs(dt) * integrator.eigen_est * (onemq + sqrt(1 + q * (0.44256220745562963 + q))))) + 1
        mdeg2 = mdeg^2
        tsw0 = constantcache.tsw0
        acoshtsw0 = constantcache.acoshtsw0
        acoshtsw0dm = acoshtsw0 / mdeg
        w0 = cosh(acoshtsw0dm)
        w0sq = w0^2
        w0sqm1 = w0sq - 1
        dtsw0 = mdeg * sinh(acoshtsw0) / sqrt(w0sqm1)
        d2tsw0 = (mdeg2 * tsw0 - w0 * dtsw0) / w0sqm1
        d3tsw0 = ((1 + 2 * w0sq + mdeg2 * w0sqm1) * dtsw0 - 3 * mdeg2 * w0 * tsw0) / (w0sqm1^2)
        w1 = (onemq * d2tsw0 + sqrt((onemq * d2tsw0)^2 + 4 * q * dtsw0 * d3tsw0)) / (2 * d3tsw0)

        b1 = 0.25 * sinh((mdeg - 2) * acoshtsw0dm) / sinh((mdeg - 1) * acoshtsw0dm)
        b = 15 / ((8 * w0)^2)
        b2 = b1
    end

    # stage-1
    @.. broadcast = false gprev2 = uprev
    μs = w1 * b1
    @.. broadcast = false gprev = uprev + dt * μs * fsalfirst
    th2 = zero(eltype(u))
    th1 = μs
    z1 = w0
    z2 = one(eltype(u))
    dz1 = one(eltype(u))
    dz2 = zero(eltype(u))
    d2z1 = zero(eltype(u))
    d2z2 = zero(eltype(u))
    z = 2 * w0 * z1 - z2
    dz = 2 * w0 * dz1 - dz2 + 2 * z1
    d2z = 2 * w0 * d2z1 - d2z2 + 4 * dz1

    # stage 2 - mdeg
    for iter in 2:mdeg
        νs = 1 - z1 * b1
        μ = (2 * w0 * b) / b1
        ν = -b / b2
        μs = μ * w1 / w0
        f(k, gprev, p, t + dt * th1)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false u = μ * gprev + ν * gprev2 + (1 - μ - ν) * uprev +
            dt * μs * (k - νs * fsalfirst)
        if (iter < mdeg)
            @.. broadcast = false gprev2 = gprev
            @.. broadcast = false gprev = u
            th = μ * th1 + ν * th2 + μs * (1 - νs)
            th2 = th1
            th1 = th
            b2 = b1
            b1 = b
            z2 = z1
            z1 = z
            dz2 = dz1
            dz1 = dz
            d2z2 = d2z1
            d2z1 = d2z

            z = 2 * w0 * z1 - z2
            dz = 2 * w0 * dz1 - dz2 + 2 * z1
            d2z = 2 * w0 * d2z1 - d2z2 + 4 * dz1

            b = d2z / (dz^2)
        end
    end

    if (!rkcstep)
        a = onepq / (q * dtsw0 + w1 * d2tsw0)
        g = (a * dtsw0 - 1) / q
        b = a / (b * w1)
        μ = 1 - g - b
        @.. broadcast = false u = μ * uprev + g * uprev2 + b * u
    end

    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # error estimate
    if integrator.opts.adaptive
        if (rkcstep)
            @.. broadcast = false tmp = 1.2 * (uprev - u) + 0.6 * dt * (fsalfirst + integrator.fsallast)
        else
            @.. broadcast = false tmp = 0.6 * (uprev / q - uprev2 / (q * onepq2) - u * (2 + q) / onepq2 + dt * integrator.fsallast / onepq)
        end
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end


function initialize!(integrator, cache::RKL1ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::RKL1Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::RKL1ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    # choose s from the eigenvalue estimate; min/max_stage are odd from RKL1(; ...) constructor
    s = ceil(Int, 0.5 * (-1.0 + sqrt(1.0 + 4.0 * abs(dt) * integrator.eigen_est)))
    s = max(s, cache.min_stage)
    s = isodd(s) ? s : s + 1
    s = min(s, cache.max_stage)
    cache.mdeg = s

    w1 = 2.0 / (s^2 + s)

    # stage 1
    uᵢ₋₂ = uprev
    uᵢ₋₁ = uprev + (dt * w1) * fsalfirst

    # stages 2 to s
    for j in 2:s
        μⱼ = (2j - 1) / j
        νⱼ = -(j - 1) / j
        μ̃ⱼ = μⱼ * w1
        fYm1 = f(uᵢ₋₁, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        u = μⱼ * uᵢ₋₁ + νⱼ * uᵢ₋₂ + μ̃ⱼ * dt * fYm1
        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = u
    end

    # error estimate from a single forward-Euler step
    if integrator.opts.adaptive
        tmp = u - (uprev + dt * fsalfirst)
        atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
        OrdinaryDiffEqCore.set_EEst!(
            integrator,
            integrator.opts.internalnorm(atmp, t) / s
        )
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RKL1Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; uᵢ₋₁, uᵢ₋₂, tmp, k, atmp) = cache
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    # choose s from the eigenvalue estimate based on the RKL1 stability bound
    s = ceil(Int, 0.5 * (-1.0 + sqrt(1.0 + 4.0 * abs(dt) * integrator.eigen_est)))
    s = max(s, ccache.min_stage)
    s = isodd(s) ? s : s + 1
    s = min(s, ccache.max_stage)
    ccache.mdeg = s

    w1 = 2.0 / (s^2 + s)

    # stage 1
    @.. broadcast = false uᵢ₋₂ = uprev
    @.. broadcast = false uᵢ₋₁ = uprev + (dt * w1) * fsalfirst

    # stages 2 to s
    for j in 2:s
        μⱼ = (2j - 1) / j
        νⱼ = -(j - 1) / j
        μ̃ⱼ = μⱼ * w1
        f(k, uᵢ₋₁, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false u = μⱼ * uᵢ₋₁ + νⱼ * uᵢ₋₂ + (dt * μ̃ⱼ) * k
        @.. broadcast = false uᵢ₋₂ = uᵢ₋₁
        @.. broadcast = false uᵢ₋₁ = u
    end

    if integrator.opts.adaptive
        @.. broadcast = false tmp = u - (uprev + dt * fsalfirst)
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(
            integrator,
            integrator.opts.internalnorm(atmp, t) / s
        )
    end

    integrator.k[1] = integrator.fsalfirst
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[2] = integrator.fsallast
end


function initialize!(integrator, cache::RKL2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::RKL2Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end


@muladd function perform_step!(integrator, cache::RKL2ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    # choose s from the eigenvalue estimate from the RKL2 stability bound
    s = ceil(Int, 0.5 * (-1.0 + sqrt(9.0 + 8.0 * abs(dt) * integrator.eigen_est)))
    s = max(s, cache.min_stage)
    s = isodd(s) ? s : s + 1
    s = min(s, cache.max_stage)
    cache.mdeg = s
    w1 = 4.0 / (s^2 + s - 2)

    # stage 1
    μ̃₁ = (1.0 / 3.0) * w1
    uᵢ₋₂ = uprev
    uᵢ₋₁ = uprev + (dt * μ̃₁) * fsalfirst

    # stages 2 to s
    for j in 2:s
        bj = j <= 2 ? 1.0 / 3.0 : (j * j + j - 2) / (2 * j * (j + 1))
        jm1 = j - 1
        bjm1 = jm1 <= 2 ? 1.0 / 3.0 : (jm1 * jm1 + jm1 - 2) / (2 * jm1 * (jm1 + 1))
        jm2 = j - 2
        bjm2 = jm2 <= 2 ? 1.0 / 3.0 : (jm2 * jm2 + jm2 - 2) / (2 * jm2 * (jm2 + 1))
        aj = 1.0 - bj

        μⱼ = (2j - 1) / j * bj / bjm1
        νⱼ = -(j - 1) / j * bj / bjm2
        μ̃ⱼ = μⱼ * w1

        ajm1 = 1.0 - bjm1
        γ̃ⱼ = -ajm1 * μ̃ⱼ

        fYm1 = f(uᵢ₋₁, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        u = μⱼ * uᵢ₋₁ + νⱼ * uᵢ₋₂ + (1.0 - μⱼ - νⱼ) * uprev +
            dt * μ̃ⱼ * fYm1 + dt * γ̃ⱼ * fsalfirst

        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = u
    end

    # error estimate based on a one step prediction
    if integrator.opts.adaptive
        tmp = u - (uprev + dt * fsalfirst)
        atmp = calculate_residuals(
            tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(
            integrator,
            integrator.opts.internalnorm(atmp, t) / s
        )
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RKL2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; uᵢ₋₁, uᵢ₋₂, tmp, k, atmp) = cache
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    s = ceil(Int, 0.5 * (-1.0 + sqrt(9.0 + 8.0 * abs(dt) * integrator.eigen_est)))
    s = max(s, ccache.min_stage)
    s = isodd(s) ? s : s + 1
    s = min(s, ccache.max_stage)
    ccache.mdeg = s
    w1 = 4.0 / (s^2 + s - 2)

    # stage 1 uses previous stage final values for legendre polynomial
    μ̃₁ = (1.0 / 3.0) * w1
    @.. broadcast = false uᵢ₋₂ = uprev
    @.. broadcast = false uᵢ₋₁ = uprev + (dt * μ̃₁) * fsalfirst

    # stages 2 to s
    for j in 2:s
        bj = j <= 2 ? 1.0 / 3.0 : (j * j + j - 2) / (2 * j * (j + 1))
        jm1 = j - 1
        bjm1 = jm1 <= 2 ? 1.0 / 3.0 : (jm1 * jm1 + jm1 - 2) / (2 * jm1 * (jm1 + 1))
        jm2 = j - 2
        bjm2 = jm2 <= 2 ? 1.0 / 3.0 : (jm2 * jm2 + jm2 - 2) / (2 * jm2 * (jm2 + 1))
        ajm1 = 1.0 - bjm1

        μⱼ = (2j - 1) / j * bj / bjm1
        νⱼ = -(j - 1) / j * bj / bjm2
        μ̃ⱼ = μⱼ * w1
        γ̃ⱼ = -ajm1 * μ̃ⱼ

        f(k, uᵢ₋₁, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false u = μⱼ * uᵢ₋₁ + νⱼ * uᵢ₋₂ +
            (1.0 - μⱼ - νⱼ) * uprev +
            (dt * μ̃ⱼ) * k +
            (dt * γ̃ⱼ) * fsalfirst
        @.. broadcast = false uᵢ₋₂ = uᵢ₋₁
        @.. broadcast = false uᵢ₋₁ = u
    end

    if integrator.opts.adaptive
        @.. broadcast = false tmp = u - (uprev + dt * fsalfirst)
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t) / s)
    end

    integrator.k[1] = integrator.fsalfirst
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::RKG1ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::RKG1Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::RKG1ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    # stability bound for RKG1 is dt*λ ≤ s(s+3)/4
    # solving s² + 3s - 4*dt*λ ≥ 0 → s ≥ (-3 + sqrt(9 + 16*dt*λ))/2
    s = ceil(Int, 0.5 * (-3.0 + sqrt(9.0 + 16.0 * abs(dt) * integrator.eigen_est)))
    s = max(s, cache.min_stage)
    s = min(s, cache.max_stage)
    cache.mdeg = s
    w1 = 4.0 / (s * (s + 3))
    rkg1_b(j) = 2.0 / ((j + 1) * (j + 2))

    # first stage coefficients
    b0 = rkg1_b(0)
    b1 = rkg1_b(1)
    μ̃₁ = (3.0 * b1 / b0) * w1
    uᵢ₋₂ = uprev
    uᵢ₋₁ = uprev + (dt * μ̃₁) * fsalfirst

    # stages 2 to s
    for j in 2:s
        bj = rkg1_b(j)
        bjm1 = rkg1_b(j - 1)
        bjm2 = rkg1_b(j - 2)

        # μⱼ = (2j+1)/j * bⱼ/bⱼ₋₁
        # νⱼ = -(j+1)/j * bⱼ/bⱼ₋₂
        μⱼ = (2j + 1) / j * bj / bjm1
        νⱼ = -(j + 1) / j * bj / bjm2
        μ̃ⱼ = μⱼ * w1

        fYm1 = f(uᵢ₋₁, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        u = μⱼ * uᵢ₋₁ + νⱼ * uᵢ₋₂ + μ̃ⱼ * dt * fYm1
        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = u
    end

    if integrator.opts.adaptive
        tmp = u - (uprev + dt * fsalfirst)
        atmp = calculate_residuals(
            tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(
            integrator,
            integrator.opts.internalnorm(atmp, t) / s
        )
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RKG1Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; uᵢ₋₁, uᵢ₋₂, tmp, k, atmp) = cache
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    # choose s from the eigenvalue estimate from the RKG1 stability bound
    s = ceil(Int, 0.5 * (-3.0 + sqrt(9.0 + 16.0 * abs(dt) * integrator.eigen_est)))
    s = max(s, ccache.min_stage)
    s = min(s, ccache.max_stage)
    ccache.mdeg = s

    w1 = 4.0 / (s * (s + 3))
    rkg1_b(j) = 2.0 / ((j + 1) * (j + 2))

    # first stage coefficients
    b0 = rkg1_b(0); b1 = rkg1_b(1)
    μ̃₁ = (3.0 * b1 / b0) * w1

    @.. broadcast = false uᵢ₋₂ = uprev
    @.. broadcast = false uᵢ₋₁ = uprev + (dt * μ̃₁) * fsalfirst

    # stages 2 to s
    for j in 2:s
        bj = rkg1_b(j)
        bjm1 = rkg1_b(j - 1)
        bjm2 = rkg1_b(j - 2)
        μⱼ = (2j + 1) / j * bj / bjm1
        νⱼ = -(j + 1) / j * bj / bjm2
        μ̃ⱼ = μⱼ * w1

        f(k, uᵢ₋₁, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false u = μⱼ * uᵢ₋₁ + νⱼ * uᵢ₋₂ + (dt * μ̃ⱼ) * k
        @.. broadcast = false uᵢ₋₂ = uᵢ₋₁
        @.. broadcast = false uᵢ₋₁ = u
    end

    if integrator.opts.adaptive
        @.. broadcast = false tmp = u - (uprev + dt * fsalfirst)
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(
            integrator,
            integrator.opts.internalnorm(atmp, t) / s
        )
    end

    integrator.k[1] = integrator.fsalfirst
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::RKG2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::RKG2Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::RKG2ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    # stability bound for RKG2 is dt*λ ≤ (s+4)(s-1)/6
    # Solve: s² + 3s - (4 + 6*dt*λ) ≥ 0 → s ≥ (-3 + sqrt(25 + 24*dt*λ))/2
    s = ceil(Int, 0.5 * (-3.0 + sqrt(25.0 + 24.0 * abs(dt) * integrator.eigen_est)))
    s = max(s, cache.min_stage)
    s = min(s, cache.max_stage)
    cache.mdeg = s
    w1 = 6.0 / ((s + 4) * (s - 1))

    # first stage coefficients
    μ̃₁ = w1
    uᵢ₋₂ = uprev
    uᵢ₋₁ = uprev + (dt * μ̃₁) * fsalfirst

    # stages 2 to s: b_0=1, b_1=1/3, b_j=4(j-1)(j+4)/(3j(j+1)(j+2)(j+3)) for j>=2
    # a_{j-1} = 1 - j(j+1)/2 * b_{j-1}
    for j in 2:s
        bj = 4.0 * (j - 1) * (j + 4) / (3.0 * j * (j + 1) * (j + 2) * (j + 3))
        bjm1 = j == 2 ? (1.0 / 3.0) :
            4.0 * (j - 2) * (j + 3) / (3.0 * (j - 1) * j * (j + 1) * (j + 2))
        bjm2 = j == 2 ? 1.0 :
            j == 3 ? (1.0 / 3.0) :
            4.0 * (j - 3) * (j + 2) / (3.0 * (j - 2) * (j - 1) * j * (j + 1))
        ajm1 = 1.0 - j * (j + 1) / 2.0 * bjm1

        μⱼ = (2j + 1) / j * bj / bjm1
        νⱼ = -(j + 1) / j * bj / bjm2
        μ̃ⱼ = μⱼ * w1
        γ̃ⱼ = -μ̃ⱼ * ajm1

        fYm1 = f(uᵢ₋₁, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        u = μⱼ * uᵢ₋₁ + νⱼ * uᵢ₋₂ + (1.0 - μⱼ - νⱼ) * uprev +
            dt * μ̃ⱼ * fYm1 + dt * γ̃ⱼ * fsalfirst

        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = u
    end

    if integrator.opts.adaptive
        tmp = u - (uprev + dt * fsalfirst)
        atmp = calculate_residuals(
            tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(
            integrator,
            integrator.opts.internalnorm(atmp, t) / s
        )
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RKG2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; uᵢ₋₁, uᵢ₋₂, tmp, k, atmp) = cache
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    # choose s from the eigenvalue estimate from the RKG2 stability bound
    s = ceil(Int, 0.5 * (-3.0 + sqrt(25.0 + 24.0 * abs(dt) * integrator.eigen_est)))
    s = max(s, ccache.min_stage)
    s = min(s, ccache.max_stage)
    ccache.mdeg = s
    w1 = 6.0 / ((s + 4) * (s - 1))
    μ̃₁ = w1

    @.. broadcast = false uᵢ₋₂ = uprev
    @.. broadcast = false uᵢ₋₁ = uprev + (dt * μ̃₁) * fsalfirst

    # stages 2 to s
    for j in 2:s
        bj = 4.0 * (j - 1) * (j + 4) / (3.0 * j * (j + 1) * (j + 2) * (j + 3))
        bjm1 = j == 2 ? (1.0 / 3.0) :
            4.0 * (j - 2) * (j + 3) / (3.0 * (j - 1) * j * (j + 1) * (j + 2))
        bjm2 = j == 2 ? 1.0 :
            j == 3 ? (1.0 / 3.0) :
            4.0 * (j - 3) * (j + 2) / (3.0 * (j - 2) * (j - 1) * j * (j + 1))
        ajm1 = 1.0 - j * (j + 1) / 2.0 * bjm1

        μⱼ = (2j + 1) / j * bj / bjm1
        νⱼ = -(j + 1) / j * bj / bjm2
        μ̃ⱼ = μⱼ * w1
        γ̃ⱼ = -μ̃ⱼ * ajm1

        f(k, uᵢ₋₁, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        @.. broadcast = false u = μⱼ * uᵢ₋₁ + νⱼ * uᵢ₋₂ +
            (1.0 - μⱼ - νⱼ) * uprev +
            (dt * μ̃ⱼ) * k +
            (dt * γ̃ⱼ) * fsalfirst
        @.. broadcast = false uᵢ₋₂ = uᵢ₋₁
        @.. broadcast = false uᵢ₋₁ = u
    end

    if integrator.opts.adaptive
        @.. broadcast = false tmp = u - (uprev + dt * fsalfirst)
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(
            integrator, integrator.opts.internalnorm(atmp, t) / s
        )
    end

    integrator.k[1] = integrator.fsalfirst
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::RKMC2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::RKMC2Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::RKMC2ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    # selects stage count from eigenvalue estimate 
    mdeg = max(3, ceil(Int, -0.8306782178712795 + 1.8547887825836553 * (abs(dt) * integrator.eigen_est)^0.533871357807877))
    mdeg = min(max(mdeg, cache.min_stage), cache.max_stage)

    # solves the equaton for w0, the parameter used to shift the Chebyshev polynomial to be monotone, using bisection on α = acosh(w0)
    if mdeg != cache.mdeg
        cache.mdeg = mdeg
        s = mdeg
        sign_s = iseven(s) ? 1 : -1
        αlo, αhi = 1e-10, 2.0
        for _ in 1:50
            α = (αlo + αhi) / 2
            Ts    = cosh(s * α)
            Tsm1  = cosh((s - 1) * α)
            Tsm2  = cosh((s - 2) * α)
            dTsm1 = (s - 1) * sinh((s - 1) * α) / sinh(α)
            fval  = 1 + sign_s / (s * (s - 2)) + cosh(α) + Ts / (2s) -
                Tsm2 / (2 * (s - 2)) - (1 + Tsm1)^2 / dTsm1
            fval > 0 ? (αlo = α) : (αhi = α)
        end
        α     = (αlo + αhi) / 2
        Tsm1  = cosh((s - 1) * α)
        dTsm1 = (s - 1) * sinh((s - 1) * α) / sinh(α)
        cache.w0 = cosh(α)
        cache.w1 = (1 + Tsm1) / dTsm1
    end
    w0, w1 = cache.w0, cache.w1

    # initializes the b value tracking: bⱼ = 1/(1+Tⱼ(w0)) via Chebyshev recurrence on Tⱼ(w0)
    Tj₋₂ = one(w0)        
    Tj₋₁ = w0              
    bj₋₂ = 1 / (1 + Tj₋₂) 
    bj₋₁ = 1 / (1 + Tj₋₁) 

    # stage 1: Y₁ = y₀ + h·b₁·w₁·F₀
    gprev2 = copy(uprev)
    μs = bj₋₁ * w1
    gprev = uprev + dt * μs * fsalfirst
    th2 = zero(eltype(u))
    th1 = μs
    bs = bj₋₁

    # stages 2 to s using the shifted Chebyshev reccurance
    for j in 2:mdeg
        Tj = 2 * w0 * Tj₋₁ - Tj₋₂
        bj = 1 / (1 + Tj)
        bs = bj
        μ  = 2 * w0 * bj / bj₋₁
        ν  = -bj / bj₋₂
        μ̃  = 2 * w1 * bj / bj₋₁
        u  = f(gprev, p, t + dt * th1)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        u  = (1 - μ - ν) * uprev + μ * gprev + ν * gprev2 + dt * μ̃ * (u - bj₋₁ * fsalfirst)
        th = μ * th1 + ν * th2 + μ̃ * (1 - bj₋₁)
        if j < mdeg
            gprev2 = gprev
            gprev  = u
            th2    = th1
            th1    = th
            Tj₋₂   = Tj₋₁
            Tj₋₁   = Tj
            bj₋₂   = bj₋₁
            bj₋₁   = bj
        end
    end

    # final solution combination of the stages 
    γs = bj₋₁ / (2 * mdeg * w1)
    δs = -bj₋₁ / (2 * (mdeg - 2) * w1)
    u = (1 - γs / bs - δs / bj₋₂) * uprev + (γs / bs) * u + (δs / bj₋₂) * gprev2 + dt * bj₋₁ * fsalfirst
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # error estimate according to the paper is (1/10)*(y₀ - y₁ + dt*f(t+h,y₁))
    if integrator.opts.adaptive
        tmp = (uprev - u + dt * integrator.fsallast) / 10
        atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RKMC2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    (; k, tmp, gprev, gprev2, atmp) = cache
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    mdeg = max(3, ceil(Int, -0.8306782178712795 + 1.8547887825836553 * (abs(dt) * integrator.eigen_est)^0.533871357807877))
    mdeg = min(max(mdeg, ccache.min_stage), ccache.max_stage)

    if mdeg != ccache.mdeg
        ccache.mdeg = mdeg
        s = mdeg
        sign_s = iseven(s) ? 1 : -1
        αlo, αhi = 1e-10, 2.0
        for _ in 1:50
            α = (αlo + αhi) / 2
            Ts    = cosh(s * α)
            Tsm1  = cosh((s - 1) * α)
            Tsm2  = cosh((s - 2) * α)
            dTsm1 = (s - 1) * sinh((s - 1) * α) / sinh(α)
            fval  = 1 + sign_s / (s * (s - 2)) + cosh(α) + Ts / (2s) - Tsm2 / (2 * (s - 2)) - (1 + Tsm1)^2 / dTsm1
            fval > 0 ? (αlo = α) : (αhi = α)
        end
        α     = (αlo + αhi) / 2
        Tsm1  = cosh((s - 1) * α)
        dTsm1 = (s - 1) * sinh((s - 1) * α) / sinh(α)
        ccache.w0 = cosh(α)
        ccache.w1 = (1 + Tsm1) / dTsm1
    end
    w0, w1 = ccache.w0, ccache.w1

    Tj₋₂ = one(w0)
    Tj₋₁ = w0
    bj₋₂ = 1 / (1 + Tj₋₂)
    bj₋₁ = 1 / (1 + Tj₋₁)

    @.. broadcast = false gprev2 = uprev
    μs = bj₋₁ * w1
    @.. broadcast = false gprev = uprev + dt * μs * fsalfirst
    th2 = zero(eltype(u))
    th1 = μs
    bs = bj₋₁

    for j in 2:mdeg
        Tj = 2 * w0 * Tj₋₁ - Tj₋₂
        bj = 1 / (1 + Tj)
        bs = bj
        μ  = 2 * w0 * bj / bj₋₁
        ν  = -bj / bj₋₂
        μ̃  = 2 * w1 * bj / bj₋₁
        f(k, gprev, p, t + dt * th1)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false u = (1 - μ - ν) * uprev + μ * gprev + ν * gprev2 + dt * μ̃ * (k - bj₋₁ * fsalfirst)
        th = μ * th1 + ν * th2 + μ̃ * (1 - bj₋₁)
        if j < mdeg
            @.. broadcast = false gprev2 = gprev
            @.. broadcast = false gprev  = u
            th2  = th1
            th1  = th
            Tj₋₂ = Tj₋₁
            Tj₋₁ = Tj
            bj₋₂ = bj₋₁
            bj₋₁ = bj
        end
    end

    # final solution combination of the stages
    γs = bj₋₁ / (2 * mdeg * w1)
    δs = -bj₋₁ / (2 * (mdeg - 2) * w1)
    @.. broadcast = false u = (1 - γs / bs - δs / bj₋₂) * uprev + (γs / bs) * u + (δs / bj₋₂) * gprev2 + dt * bj₋₁ * fsalfirst

    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        @.. broadcast = false tmp = (uprev - u + dt * integrator.fsallast) / 10
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end
