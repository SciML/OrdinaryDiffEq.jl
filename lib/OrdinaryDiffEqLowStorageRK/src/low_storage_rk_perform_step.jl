# 2N low storage methods
function initialize!(integrator, cache::LowStorageRK2NConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(
        integrator, cache::LowStorageRK2NConstantCache,
        repeat_step = false
    )
    (; t, dt, u, f, p) = integrator
    (; A2end, B1, B2end, c2end) = cache

    # u1
    tmp = dt * integrator.fsalfirst
    u = u + B1 * tmp

    # other stages
    for i in eachindex(A2end)
        k = f(u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tmp = A2end[i] * tmp + dt * k
        u = u + B2end[i] * tmp
    end

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.fsalfirst = f(u, p, t + dt) # For interpolation, then FSAL'd
    integrator.u = u
end

get_fsalfirstlast(cache::LowStorageRK2NCache, u) = (nothing, nothing)

function initialize!(integrator, cache::LowStorageRK2NCache)
    (; k, tmp, williamson_condition) = cache
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = k
    integrator.f(k, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::LowStorageRK2NCache, repeat_step = false)
    (; t, dt, u, f, p) = integrator
    (; k, tmp, williamson_condition, stage_limiter!, step_limiter!, thread) = cache
    (; A2end, B1, B2end, c2end) = cache.tab

    # u1
    f(k, u, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false thread = thread tmp = dt * k
    @.. broadcast = false thread = thread u = u + B1 * tmp
    # other stages
    for i in eachindex(A2end)
        if williamson_condition
            f(ArrayFuse(tmp, u, (A2end[i], dt, B2end[i])), u, p, t + c2end[i] * dt)
        else
            @.. broadcast = false thread = thread tmp = A2end[i] * tmp
            stage_limiter!(u, integrator, p, t + c2end[i] * dt)
            f(k, u, p, t + c2end[i] * dt)
            @.. broadcast = false thread = thread tmp = tmp + dt * k
            @.. broadcast = false thread = thread u = u + B2end[i] * tmp
        end
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
end

# 2C low storage methods
function initialize!(integrator, cache::LowStorageRK2CConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(
        integrator, cache::LowStorageRK2CConstantCache,
        repeat_step = false
    )
    (; t, dt, u, f, p) = integrator
    (; A2end, B1, B2end, c2end) = cache

    # u1
    k = integrator.fsalfirst = f(u, p, t)
    integrator.k[1] = integrator.fsalfirst
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    u = u + B1 * dt * k

    # other stages
    for i in eachindex(A2end)
        tmp = u + A2end[i] * dt * k
        k = f(tmp, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        u = u + B2end[i] * dt * k
    end

    integrator.u = u
end

function initialize!(integrator, cache::LowStorageRK2CCache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::LowStorageRK2CCache, repeat_step = false)
    (; t, dt, u, f, p) = integrator
    (; k, fsalfirst, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; A2end, B1, B2end, c2end) = cache.tab

    # u1
    @.. broadcast = false thread = thread k = integrator.fsalfirst
    @.. broadcast = false thread = thread u = u + B1 * dt * k

    # other stages
    for i in eachindex(A2end)
        @.. broadcast = false thread = thread tmp = u + A2end[i] * dt * k
        f(k, tmp, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false thread = thread u = u + B2end[i] * dt * k
    end
    step_limiter!(u, integrator, p, t + dt)
    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

# 3S low storage methods
function initialize!(integrator, cache::LowStorageRK3SConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(
        integrator, cache::LowStorageRK3SConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end) = cache

    # u1
    tmp = u
    u = tmp + β1 * dt * integrator.fsalfirst

    # other stages
    for i in eachindex(γ12end)
        k = f(u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tmp = tmp + δ2end[i] * u
        u = γ12end[i] * u + γ22end[i] * tmp + γ32end[i] * uprev + β2end[i] * dt * k
    end

    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.u = u
end

function initialize!(integrator, cache::LowStorageRK3SCache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::LowStorageRK3SCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end) = cache.tab

    # u1
    @.. broadcast = false thread = thread tmp = u
    @.. broadcast = false thread = thread u = tmp + β1 * dt * integrator.fsalfirst

    # other stages
    for i in eachindex(γ12end)
        f(k, u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false thread = thread tmp = tmp + δ2end[i] * u
        @.. broadcast = false thread = thread u = γ12end[i] * u + γ22end[i] * tmp +
            γ32end[i] * uprev +
            β2end[i] * dt * k
    end

    step_limiter!(u, integrator, p, t + dt)
    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

# 3S+ low storage methods: 3S methods adding another memory location for the embedded method (non-FSAL version)
function initialize!(integrator, cache::LowStorageRK3SpConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(
        integrator, cache::LowStorageRK3SpConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end) = cache

    # u1
    integrator.fsalfirst = f(uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    tmp = uprev
    u = tmp + β1 * dt * integrator.fsalfirst
    # Initialize utilde for JET
    utilde = u
    if integrator.opts.adaptive
        utilde = bhat1 * dt * integrator.fsalfirst
    end

    # other stages
    for i in eachindex(γ12end)
        k = f(u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tmp = tmp + δ2end[i] * u
        u = γ12end[i] * u + γ22end[i] * tmp + γ32end[i] * uprev + β2end[i] * dt * k
        if integrator.opts.adaptive
            utilde = utilde + bhat2end[i] * dt * k
        end
    end

    if integrator.opts.adaptive
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.u = u
end

function initialize!(integrator, cache::LowStorageRK3SpCache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    return integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::LowStorageRK3SpCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, tmp, utilde, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end) = cache.tab

    # u1
    f(integrator.fsalfirst, uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false thread = thread tmp = uprev
    @.. broadcast = false thread = thread u = tmp + β1 * dt * integrator.fsalfirst
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = bhat1 * dt * integrator.fsalfirst
    end

    # other stages
    for i in eachindex(γ12end)
        stage_limiter!(u, integrator, p, t + c2end[i] * dt)
        f(k, u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false thread = thread tmp = tmp + δ2end[i] * u
        @.. broadcast = false thread = thread u = γ12end[i] * u + γ22end[i] * tmp +
            γ32end[i] * uprev + β2end[i] * dt * k
        if integrator.opts.adaptive
            @.. broadcast = false thread = thread utilde = utilde + bhat2end[i] * dt * k
        end
    end

    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

# 3S+ FSAL low storage methods: 3S methods adding another memory location for the embedded method (FSAL version)
function initialize!(integrator, cache::LowStorageRK3SpFSALConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(
        integrator, cache::LowStorageRK3SpFSALConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end, bhatfsal) = cache

    # u1
    tmp = uprev
    u = tmp + β1 * dt * integrator.fsalfirst
    # Initialize utilde for JET
    utilde = u
    if integrator.opts.adaptive
        utilde = bhat1 * dt * integrator.fsalfirst
    end

    # other stages
    for i in eachindex(γ12end)
        k = f(u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tmp = tmp + δ2end[i] * u
        u = γ12end[i] * u + γ22end[i] * tmp + γ32end[i] * uprev + β2end[i] * dt * k
        if integrator.opts.adaptive
            utilde = utilde + bhat2end[i] * dt * k
        end
    end

    # FSAL
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        utilde = utilde + bhatfsal * dt * integrator.fsallast
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::LowStorageRK3SpFSALCache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(
        integrator, cache::LowStorageRK3SpFSALCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; k, tmp, utilde, atmp, stage_limiter!, step_limiter!, thread) = cache
    (;
        γ12end, γ22end, γ32end, δ2end, β1, β2end,
        c2end, bhat1, bhat2end, bhatfsal,
    ) = cache.tab

    # u1
    @.. broadcast = false thread = thread tmp = uprev
    @.. broadcast = false thread = thread u = tmp + β1 * dt * integrator.fsalfirst
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = bhat1 * dt * integrator.fsalfirst
    end

    # other stages
    for i in eachindex(γ12end)
        stage_limiter!(u, integrator, p, t + c2end[i] * dt)
        f(k, u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false thread = thread tmp = tmp + δ2end[i] * u
        @.. broadcast = false thread = thread u = γ12end[i] * u + γ22end[i] * tmp +
            γ32end[i] * uprev + β2end[i] * dt * k
        if integrator.opts.adaptive
            @.. broadcast = false thread = thread utilde = utilde + bhat2end[i] * dt * k
        end
    end

    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    # FSAL
    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = utilde + bhatfsal * dt * k
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

# 2R+ low storage methods
function initialize!(integrator, cache::LowStorageRK2RPConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(
        integrator, cache::LowStorageRK2RPConstantCache,
        repeat_step = false
    )
    (; t, dt, u, uprev, f, fsalfirst, p) = integrator
    (; Aᵢ, Bₗ, B̂ₗ, Bᵢ, B̂ᵢ, Cᵢ) = cache

    k = fsalfirst
    # Initialize tmp for JET
    tmp = uprev
    integrator.opts.adaptive && (tmp = zero(uprev))

    #stages 1 to s-1
    for i in eachindex(Aᵢ)
        integrator.opts.adaptive && (tmp = tmp + (Bᵢ[i] - B̂ᵢ[i]) * dt * k)
        gprev = u + Aᵢ[i] * dt * k
        u = u + Bᵢ[i] * dt * k
        k = f(gprev, p, t + Cᵢ[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    #last stage
    integrator.opts.adaptive && (tmp = tmp + (Bₗ - B̂ₗ) * dt * k)
    u = u + Bₗ * dt * k

    #Error estimate
    if integrator.opts.adaptive
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

function initialize!(integrator, cache::LowStorageRK2RPCache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::LowStorageRK2RPCache, repeat_step = false)
    (; t, dt, u, uprev, f, fsalfirst, p) = integrator
    (; k, gprev, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; Aᵢ, Bₗ, B̂ₗ, Bᵢ, B̂ᵢ, Cᵢ) = cache.tab

    @.. broadcast = false thread = thread k = fsalfirst
    integrator.opts.adaptive && (@.. broadcast = false tmp = zero(uprev))

    #stages 1 to s-1
    for i in eachindex(Aᵢ)
        integrator.opts.adaptive &&
            (@.. broadcast = false thread = thread tmp = tmp + (Bᵢ[i] - B̂ᵢ[i]) * dt * k)
        @.. broadcast = false thread = thread gprev = u + Aᵢ[i] * dt * k
        @.. broadcast = false thread = thread u = u + Bᵢ[i] * dt * k
        f(k, gprev, p, t + Cᵢ[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    #last stage
    integrator.opts.adaptive &&
        (@.. broadcast = false thread = thread tmp = tmp + (Bₗ - B̂ₗ) * dt * k)
    @.. broadcast = false thread = thread u = u + Bₗ * dt * k

    #Error estimate
    if integrator.opts.adaptive
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    step_limiter!(u, integrator, p, t + dt)
    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

# 3R+ low storage methods
function initialize!(integrator, cache::LowStorageRK3RPConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(
        integrator, cache::LowStorageRK3RPConstantCache,
        repeat_step = false
    )
    (; t, dt, u, uprev, f, fsalfirst, p) = integrator
    (; Aᵢ₁, Aᵢ₂, Bₗ, B̂ₗ, Bᵢ, B̂ᵢ, Cᵢ) = cache

    fᵢ₋₂ = zero(fsalfirst)
    k = fsalfirst
    uᵢ₋₁ = uprev
    uᵢ₋₂ = uprev
    # Initialize tmp for JET
    tmp = uprev
    integrator.opts.adaptive && (tmp = zero(uprev))

    #stages 1 to s-1
    for i in eachindex(Aᵢ₁)
        integrator.opts.adaptive && (tmp = tmp + (Bᵢ[i] - B̂ᵢ[i]) * dt * k)
        gprev = uᵢ₋₂ + (Aᵢ₁[i] * k + Aᵢ₂[i] * fᵢ₋₂) * dt
        u = u + Bᵢ[i] * dt * k
        fᵢ₋₂ = k
        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = u
        k = f(gprev, p, t + Cᵢ[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    #last stage
    integrator.opts.adaptive && (tmp = tmp + (Bₗ - B̂ₗ) * dt * k)
    u = u + Bₗ * dt * k

    #Error estimate
    if integrator.opts.adaptive
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

function initialize!(integrator, cache::LowStorageRK3RPCache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::LowStorageRK3RPCache, repeat_step = false)
    (; t, dt, u, uprev, f, fsalfirst, p) = integrator
    (; k, uᵢ₋₁, uᵢ₋₂, gprev, fᵢ₋₂, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; Aᵢ₁, Aᵢ₂, Bₗ, B̂ₗ, Bᵢ, B̂ᵢ, Cᵢ) = cache.tab

    @.. broadcast = false thread = thread fᵢ₋₂ = zero(fsalfirst)
    @.. broadcast = false thread = thread k = fsalfirst
    integrator.opts.adaptive && (@.. broadcast = false thread = thread tmp = zero(uprev))
    @.. broadcast = false thread = thread uᵢ₋₁ = uprev
    @.. broadcast = false thread = thread uᵢ₋₂ = uprev

    #stages 1 to s-1
    for i in eachindex(Aᵢ₁)
        integrator.opts.adaptive &&
            (@.. broadcast = false thread = thread tmp = tmp + (Bᵢ[i] - B̂ᵢ[i]) * dt * k)
        @.. broadcast = false thread = thread gprev = uᵢ₋₂ + (Aᵢ₁[i] * k + Aᵢ₂[i] * fᵢ₋₂) * dt
        @.. broadcast = false thread = thread u = u + Bᵢ[i] * dt * k
        @.. broadcast = false thread = thread fᵢ₋₂ = k
        @.. broadcast = false thread = thread uᵢ₋₂ = uᵢ₋₁
        @.. broadcast = false thread = thread uᵢ₋₁ = u
        f(k, gprev, p, t + Cᵢ[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    #last stage
    integrator.opts.adaptive &&
        (@.. broadcast = false thread = thread tmp = tmp + (Bₗ - B̂ₗ) * dt * k)
    @.. broadcast = false thread = thread u = u + Bₗ * dt * k

    step_limiter!(u, integrator, p, t + dt)

    #Error estimate
    if integrator.opts.adaptive
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

# 4R+ low storage methods
function initialize!(integrator, cache::LowStorageRK4RPConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(
        integrator, cache::LowStorageRK4RPConstantCache,
        repeat_step = false
    )
    (; t, dt, u, uprev, f, fsalfirst, p) = integrator
    (; Aᵢ₁, Aᵢ₂, Aᵢ₃, Bₗ, B̂ₗ, Bᵢ, B̂ᵢ, Cᵢ) = cache

    fᵢ₋₂ = zero(fsalfirst)
    fᵢ₋₃ = zero(fsalfirst)
    k = fsalfirst
    uᵢ₋₁ = uprev
    uᵢ₋₂ = uprev
    uᵢ₋₃ = uprev
    # Initialize tmp for JET
    tmp = uprev
    integrator.opts.adaptive && (tmp = zero(uprev))

    #stages 1 to s-1
    for i in eachindex(Aᵢ₁)
        integrator.opts.adaptive && (tmp = tmp + (Bᵢ[i] - B̂ᵢ[i]) * dt * k)
        gprev = uᵢ₋₃ + (Aᵢ₁[i] * k + Aᵢ₂[i] * fᵢ₋₂ + Aᵢ₃[i] * fᵢ₋₃) * dt
        u = u + Bᵢ[i] * dt * k
        fᵢ₋₃ = fᵢ₋₂
        fᵢ₋₂ = k
        uᵢ₋₃ = uᵢ₋₂
        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = u
        k = f(gprev, p, t + Cᵢ[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    #last stage
    integrator.opts.adaptive && (tmp = tmp + (Bₗ - B̂ₗ) * dt * k)
    u = u + Bₗ * dt * k

    #Error estimate
    if integrator.opts.adaptive
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

function initialize!(integrator, cache::LowStorageRK4RPCache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::LowStorageRK4RPCache, repeat_step = false)
    (; t, dt, u, uprev, f, fsalfirst, p) = integrator
    (;
        k, uᵢ₋₁, uᵢ₋₂, uᵢ₋₃, gprev, fᵢ₋₂, fᵢ₋₃, tmp, atmp,
        stage_limiter!, step_limiter!, thread,
    ) = cache
    (; Aᵢ₁, Aᵢ₂, Aᵢ₃, Bₗ, B̂ₗ, Bᵢ, B̂ᵢ, Cᵢ) = cache.tab

    @.. broadcast = false thread = thread fᵢ₋₂ = zero(fsalfirst)
    @.. broadcast = false thread = thread fᵢ₋₃ = zero(fsalfirst)
    @.. broadcast = false thread = thread k = fsalfirst
    integrator.opts.adaptive && (@.. broadcast = false thread = thread tmp = zero(uprev))
    @.. broadcast = false thread = thread uᵢ₋₁ = uprev
    @.. broadcast = false thread = thread uᵢ₋₂ = uprev
    @.. broadcast = false thread = thread uᵢ₋₃ = uprev

    #stages 1 to s-1
    for i in eachindex(Aᵢ₁)
        integrator.opts.adaptive &&
            (@.. broadcast = false thread = thread tmp = tmp + (Bᵢ[i] - B̂ᵢ[i]) * dt * k)
        @.. broadcast = false thread = thread gprev = uᵢ₋₃ +
            (
            Aᵢ₁[i] * k + Aᵢ₂[i] * fᵢ₋₂ +
                Aᵢ₃[i] * fᵢ₋₃
        ) *
            dt
        @.. broadcast = false thread = thread u = u + Bᵢ[i] * dt * k
        @.. broadcast = false thread = thread fᵢ₋₃ = fᵢ₋₂
        @.. broadcast = false thread = thread fᵢ₋₂ = k
        @.. broadcast = false thread = thread uᵢ₋₃ = uᵢ₋₂
        @.. broadcast = false thread = thread uᵢ₋₂ = uᵢ₋₁
        @.. broadcast = false thread = thread uᵢ₋₁ = u
        f(k, gprev, p, t + Cᵢ[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    #last stage
    integrator.opts.adaptive &&
        (@.. broadcast = false thread = thread tmp = tmp + (Bₗ - B̂ₗ) * dt * k)
    @.. broadcast = false thread = thread u = u + Bₗ * dt * k

    step_limiter!(u, integrator, p, t + dt)

    #Error estimate
    if integrator.opts.adaptive
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

# 5R+ low storage methods
function initialize!(integrator, cache::LowStorageRK5RPConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(
        integrator, cache::LowStorageRK5RPConstantCache,
        repeat_step = false
    )
    (; t, dt, u, uprev, f, fsalfirst, p) = integrator
    (; Aᵢ₁, Aᵢ₂, Aᵢ₃, Aᵢ₄, Bₗ, B̂ₗ, Bᵢ, B̂ᵢ, Cᵢ) = cache

    fᵢ₋₂ = zero(fsalfirst)
    fᵢ₋₃ = zero(fsalfirst)
    fᵢ₋₄ = zero(fsalfirst)
    k = fsalfirst
    uᵢ₋₁ = uprev
    uᵢ₋₂ = uprev
    uᵢ₋₃ = uprev
    uᵢ₋₄ = uprev
    # Initialize tmp for JET
    tmp = uprev
    integrator.opts.adaptive && (tmp = zero(uprev))

    #stages 1 to s-1
    for i in eachindex(Aᵢ₁)
        integrator.opts.adaptive && (tmp = tmp + (Bᵢ[i] - B̂ᵢ[i]) * dt * k)
        gprev = uᵢ₋₄ + (Aᵢ₁[i] * k + Aᵢ₂[i] * fᵢ₋₂ + Aᵢ₃[i] * fᵢ₋₃ + Aᵢ₄[i] * fᵢ₋₄) * dt
        u = u + Bᵢ[i] * dt * k
        fᵢ₋₄ = fᵢ₋₃
        fᵢ₋₃ = fᵢ₋₂
        fᵢ₋₂ = k
        uᵢ₋₄ = uᵢ₋₃
        uᵢ₋₃ = uᵢ₋₂
        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = u
        k = f(gprev, p, t + Cᵢ[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    #last stage
    integrator.opts.adaptive && (tmp = tmp + (Bₗ - B̂ₗ) * dt * k)
    u = u + Bₗ * dt * k

    #Error estimate
    if integrator.opts.adaptive
        atmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

function initialize!(integrator, cache::LowStorageRK5RPCache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::LowStorageRK5RPCache, repeat_step = false)
    (; t, dt, u, uprev, f, fsalfirst, p) = integrator
    (;
        k, uᵢ₋₁, uᵢ₋₂, uᵢ₋₃, uᵢ₋₄, gprev, fᵢ₋₂, fᵢ₋₃, fᵢ₋₄, tmp,
        atmp, stage_limiter!, step_limiter!, thread,
    ) = cache
    (; Aᵢ₁, Aᵢ₂, Aᵢ₃, Aᵢ₄, Bₗ, B̂ₗ, Bᵢ, B̂ᵢ, Cᵢ) = cache.tab

    @.. broadcast = false thread = thread fᵢ₋₂ = zero(fsalfirst)
    @.. broadcast = false thread = thread fᵢ₋₃ = zero(fsalfirst)
    @.. broadcast = false thread = thread fᵢ₋₄ = zero(fsalfirst)
    @.. broadcast = false thread = thread k = fsalfirst
    integrator.opts.adaptive && (@.. broadcast = false thread = thread tmp = zero(uprev))
    @.. broadcast = false thread = thread uᵢ₋₁ = uprev
    @.. broadcast = false thread = thread uᵢ₋₂ = uprev
    @.. broadcast = false thread = thread uᵢ₋₃ = uprev
    @.. broadcast = false thread = thread uᵢ₋₄ = uprev

    #stages 1 to s-1
    for i in eachindex(Aᵢ₁)
        integrator.opts.adaptive &&
            (@.. broadcast = false thread = thread tmp = tmp + (Bᵢ[i] - B̂ᵢ[i]) * dt * k)
        @.. broadcast = false thread = thread gprev = uᵢ₋₄ +
            (
            Aᵢ₁[i] * k + Aᵢ₂[i] * fᵢ₋₂ +
                Aᵢ₃[i] * fᵢ₋₃ +
                Aᵢ₄[i] * fᵢ₋₄
        ) * dt
        @.. broadcast = false thread = thread u = u + Bᵢ[i] * dt * k
        @.. broadcast = false thread = thread fᵢ₋₄ = fᵢ₋₃
        @.. broadcast = false thread = thread fᵢ₋₃ = fᵢ₋₂
        @.. broadcast = false thread = thread fᵢ₋₂ = k
        @.. broadcast = false thread = thread uᵢ₋₄ = uᵢ₋₃
        @.. broadcast = false thread = thread uᵢ₋₃ = uᵢ₋₂
        @.. broadcast = false thread = thread uᵢ₋₂ = uᵢ₋₁
        @.. broadcast = false thread = thread uᵢ₋₁ = u
        f(k, gprev, p, t + Cᵢ[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    #last stage
    integrator.opts.adaptive &&
        (@.. broadcast = false thread = thread tmp = tmp + (Bₗ - B̂ₗ) * dt * k)
    @.. broadcast = false thread = thread u = u + Bₗ * dt * k

    step_limiter!(u, integrator, p, t + dt)

    #Error estimate
    if integrator.opts.adaptive
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::RK46NLCache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::RK46NLCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; α2, α3, α4, α5, α6, β1, β2, β3, β4, β5, β6, c2, c3, c4, c5, c6) = cache.tab

    # u1
    @.. broadcast = false thread = thread tmp = dt * fsalfirst
    @.. broadcast = false thread = thread u = uprev + β1 * tmp
    stage_limiter!(u, integrator, p, t + c2 * dt)
    # u2
    f(k, u, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = α2 * tmp + dt * k
    @.. broadcast = false thread = thread u = u + β2 * tmp
    stage_limiter!(u, integrator, p, t + c3 * dt)
    # u3
    f(k, u, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = α3 * tmp + dt * k
    @.. broadcast = false thread = thread u = u + β3 * tmp
    stage_limiter!(u, integrator, p, t + c4 * dt)
    # u4
    f(k, u, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = α4 * tmp + dt * k
    @.. broadcast = false thread = thread u = u + β4 * tmp
    stage_limiter!(u, integrator, p, t + c5 * dt)
    # u5 = u
    f(k, u, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = α5 * tmp + dt * k
    @.. broadcast = false thread = thread u = u + β5 * tmp
    stage_limiter!(u, integrator, p, t + c6 * dt)

    f(k, u, p, t + c6 * dt)
    @.. broadcast = false thread = thread tmp = α6 * tmp + dt * k
    @.. broadcast = false thread = thread u = u + β6 * tmp
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
end

function initialize!(integrator, cache::RK46NLConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator, cache::RK46NLConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α2, α3, α4, α5, α6, β1, β2, β3, β4, β5, β6, c2, c3, c4, c5, c6) = cache

    # u1
    tmp = dt * integrator.fsalfirst
    u = uprev + β1 * tmp
    # u2
    tmp = α2 * tmp + dt * f(u, p, t + c2 * dt)
    u = u + β2 * tmp
    # u3
    tmp = α3 * tmp + dt * f(u, p, t + c3 * dt)
    u = u + β3 * tmp
    # u4
    tmp = α4 * tmp + dt * f(u, p, t + c4 * dt)
    u = u + β4 * tmp
    # u5 = u
    tmp = α5 * tmp + dt * f(u, p, t + c5 * dt)
    u = u + β5 * tmp
    # u6
    tmp = α6 * tmp + dt * f(u, p, t + c6 * dt)
    u = u + β6 * tmp

    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    integrator.k[1] = integrator.fsalfirst
    integrator.u = u
end

function initialize!(integrator, cache::SHLDDRK52ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(
        integrator, cache::SHLDDRK52ConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; α2, α3, α4, α5, β1, β2, β3, β4, β5, c2, c3, c4, c5) = cache

    # u1
    tmp = dt * integrator.fsalfirst
    u = uprev + β1 * tmp
    # u2
    tmp = α2 * tmp + dt * f(u, p, t + c2 * dt)
    u = u + β2 * tmp
    # u3
    tmp = α3 * tmp + dt * f(u, p, t + c3 * dt)
    u = u + β3 * tmp
    # u4
    tmp = α4 * tmp + dt * f(u, p, t + c4 * dt)
    u = u + β4 * tmp
    # u5 = u
    tmp = α5 * tmp + dt * f(u, p, t + c5 * dt)
    u = u + β5 * tmp

    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
    integrator.u = u
end

function initialize!(integrator, cache::SHLDDRK52Cache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::SHLDDRK52Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; α2, α3, α4, α5, β1, β2, β3, β4, β5, c2, c3, c4, c5) = cache.tab

    # u1
    @.. thread = thread tmp = dt * fsalfirst
    @.. thread = thread u = uprev + β1 * tmp
    stage_limiter!(u, integrator, p, t + c2 * dt)
    # u2
    f(k, u, p, t + c2 * dt)
    @.. thread = thread tmp = α2 * tmp + dt * k
    @.. thread = thread u = u + β2 * tmp
    stage_limiter!(u, integrator, p, t + c3 * dt)
    # u3
    f(k, u, p, t + c3 * dt)
    @.. thread = thread tmp = α3 * tmp + dt * k
    @.. thread = thread u = u + β3 * tmp
    stage_limiter!(u, integrator, p, t + c4 * dt)
    # u4
    f(k, u, p, t + c4 * dt)
    @.. thread = thread tmp = α4 * tmp + dt * k
    @.. thread = thread u = u + β4 * tmp
    stage_limiter!(u, integrator, p, t + c5 * dt)
    # u5 = u
    f(k, u, p, t + c5 * dt)
    @.. thread = thread tmp = α5 * tmp + dt * k
    @.. thread = thread u = u + β5 * tmp
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
end

function initialize!(integrator, cache::SHLDDRK_2NConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(
        integrator, cache::SHLDDRK_2NConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (;
        α21, α31, α41, α51, β11, β21, β31, β41, β51, c21, c31, c41, c51, α22, α32,
        α42, α52, α62, β12, β22, β32, β42, β52, β62, c22, c32, c42, c52, c62,
    ) = cache

    if integrator.u_modified
        cache.step = 1
    end
    # cnt = cache.step

    if cache.step % 2 == 1
        cache.step += 1
        # u1
        tmp = dt * integrator.fsalfirst
        u = uprev + β11 * tmp
        # u2
        tmp = α21 * tmp + dt * f(u, p, t + c21 * dt)
        u = u + β21 * tmp
        # u3
        tmp = α31 * tmp + dt * f(u, p, t + c31 * dt)
        u = u + β31 * tmp
        # u4
        tmp = α41 * tmp + dt * f(u, p, t + c41 * dt)
        u = u + β41 * tmp
        # u5 = u
        tmp = α51 * tmp + dt * f(u, p, t + c51 * dt)
        u = u + β51 * tmp
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
    else
        cache.step += 1
        # u1
        tmp = dt * integrator.fsalfirst
        u = uprev + β12 * tmp
        # u2
        tmp = α22 * tmp + dt * f(u, p, t + c22 * dt)
        u = u + β22 * tmp
        # u3
        tmp = α32 * tmp + dt * f(u, p, t + c32 * dt)
        u = u + β32 * tmp
        # u4
        tmp = α42 * tmp + dt * f(u, p, t + c42 * dt)
        u = u + β42 * tmp
        # u5 = u
        tmp = α52 * tmp + dt * f(u, p, t + c52 * dt)
        u = u + β52 * tmp
        tmp = α62 * tmp + dt * f(u, p, t + c62 * dt)
        u = u + β62 * tmp
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
    end

    integrator.fsallast = f(u, p, t + dt) # For interpolation, then FSAL'd
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
end

function initialize!(integrator, cache::SHLDDRK_2NCache)
    (; k, fsalfirst) = cache

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # FSAL for interpolation
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::SHLDDRK_2NCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, tmp, stage_limiter!, step_limiter!, thread) = cache
    (;
        α21, α31, α41, α51, β11, β21, β31, β41, β51, c21, c31, c41, c51, α22, α32, α42,
        α52, α62, β12, β22, β32, β42, β52, β62, c22, c32, c42, c52, c62,
    ) = cache.tab

    if integrator.u_modified
        cache.step = 1
    end

    if cache.step % 2 == 1
        # u1
        @.. thread = thread tmp = dt * fsalfirst
        @.. thread = thread u = uprev + β11 * tmp
        stage_limiter!(u, integrator, p, t + c21 * dt)
        # u2
        f(k, u, p, t + c21 * dt)
        @.. thread = thread tmp = α21 * tmp + dt * k
        @.. thread = thread u = u + β21 * tmp
        stage_limiter!(u, integrator, p, t + c31 * dt)
        # u3
        f(k, u, p, t + c31 * dt)
        @.. thread = thread tmp = α31 * tmp + dt * k
        @.. thread = thread u = u + β31 * tmp
        stage_limiter!(u, integrator, p, t + c41 * dt)
        # u4
        f(k, u, p, t + c41 * dt)
        @.. thread = thread tmp = α41 * tmp + dt * k
        @.. thread = thread u = u + β41 * tmp
        stage_limiter!(u, integrator, p, t + c51 * dt)
        # u5 = u
        f(k, u, p, t + c51 * dt)
        @.. thread = thread tmp = α51 * tmp + dt * k
        @.. thread = thread u = u + β51 * tmp
        stage_limiter!(u, integrator, p, t + dt)
        step_limiter!(u, integrator, p, t + dt)

        f(k, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
    else
        # u1
        @.. thread = thread tmp = dt * fsalfirst
        @.. thread = thread u = uprev + β12 * tmp
        stage_limiter!(u, integrator, p, t + c22 * dt)
        # u2
        f(k, u, p, t + c22 * dt)
        @.. thread = thread tmp = α22 * tmp + dt * k
        @.. thread = thread u = u + β22 * tmp
        stage_limiter!(u, integrator, p, t + c32 * dt)
        # u3
        f(k, u, p, t + c32 * dt)
        @.. thread = thread tmp = α32 * tmp + dt * k
        @.. thread = thread u = u + β32 * tmp
        stage_limiter!(u, integrator, p, t + c42 * dt)
        # u4
        f(k, u, p, t + c42 * dt)
        @.. thread = thread tmp = α42 * tmp + dt * k
        @.. thread = thread u = u + β42 * tmp
        stage_limiter!(u, integrator, p, t + c52 * dt)
        # u5 = u
        f(k, u, p, t + c52 * dt)
        @.. thread = thread tmp = α52 * tmp + dt * k
        @.. thread = thread u = u + β52 * tmp
        stage_limiter!(u, integrator, p, t + c62 * dt)
        # u6 = u
        f(k, u, p, t + c62 * dt)
        @.. thread = thread tmp = α62 * tmp + dt * k
        @.. thread = thread u = u + β62 * tmp
        stage_limiter!(u, integrator, p, t + dt)
        step_limiter!(u, integrator, p, t + dt)

        f(k, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    end
end
