function initialize!(integrator, cache::BS3ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::BS3ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; a21, a32, a41, a42, a43, c1, c2, btilde1, btilde2, btilde3, btilde4) = cache
    k1 = integrator.fsalfirst
    a1 = dt * a21
    k2 = f(uprev + a1 * k1, p, t + c1 * dt)
    a2 = dt * a32
    k3 = f(uprev + a2 * k2, p, t + c2 * dt)
    u = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    k4 = f(u, p, t + dt)
    integrator.fsallast = k4
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    if integrator.opts.adaptive
        utilde = dt * (btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4)
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

get_fsalfirstlast(cache::BS3Cache, u) = (cache.fsalfirst, cache.k4)
function initialize!(integrator, cache::BS3Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::BS3Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k2, k3, k4, utilde, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a32, a41, a42, a43, c1, c2, btilde1, btilde2, btilde3, btilde4) = cache.tab
    # k1 = cache.fsalfirst
    k1 = integrator.fsalfirst
    a1 = dt * a21
    @.. broadcast = false thread = thread tmp = uprev + a1 * k1
    stage_limiter!(tmp, integrator, p, t + c1 * dt)
    f(k2, tmp, p, t + c1 * dt)
    a2 = dt * a32
    @.. broadcast = false thread = thread tmp = uprev + a2 * k2
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k3, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread u = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k4, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = dt * (
            btilde1 * k1 + btilde2 * k2 +
                btilde3 * k3 + btilde4 * k4
        )
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

function initialize!(integrator, cache::OwrenZen3ConstantCache)
    integrator.kshortsize = 4
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(
        integrator, cache::OwrenZen3ConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; a21, a31, a32, a41, a42, a43, c1, c2, btilde1, btilde2, btilde3) = cache
    k1 = integrator.fsalfirst
    a1 = dt * a21
    k2 = f(uprev + a1 * k1, p, t + c1 * dt)
    tmp = uprev + dt * (a31 * k1 + a32 * k2)
    k3 = f(tmp, p, t + c2 * dt)
    u = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    k4 = f(u, p, t + dt)
    integrator.fsallast = k4
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    if integrator.opts.adaptive
        utilde = dt * (btilde1 * k1 + btilde2 * k2 + btilde3 * k3)
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.u = u
end

get_fsalfirstlast(cache::OwrenZen3Cache, u) = (cache.k1, cache.k4)
function initialize!(integrator, cache::OwrenZen3Cache)
    integrator.kshortsize = 4
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::OwrenZen3Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k1, k2, k3, k4, utilde, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a31, a32, a41, a42, a43, c1, c2, btilde1, btilde2, btilde3) = cache.tab
    a1 = dt * a21
    @.. broadcast = false thread = thread tmp = uprev + a1 * k1
    stage_limiter!(tmp, integrator, p, t + c1 * dt)
    f(k2, tmp, p, t + c1 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k3, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread u = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k4, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = dt * (
            btilde1 * k1 + btilde2 * k2 +
                btilde3 * k3
        )
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

function initialize!(integrator, cache::OwrenZen4ConstantCache)
    integrator.kshortsize = 6
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(
        integrator, cache::OwrenZen4ConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, c1, c2, c3, c4, btilde1, btilde3, btilde4, btilde5) = cache
    k1 = integrator.fsalfirst
    a = dt * a21
    k2 = f(uprev + a * k1, p, t + c1 * dt)
    k3 = f(uprev + dt * (a31 * k1 + a32 * k2), p, t + c2 * dt)
    k4 = f(uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3), p, t + c3 * dt)
    k5 = f(uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4), p, t + c4 * dt)
    u = uprev + dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = f(u, p, t + dt)
    integrator.fsallast = k6
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
    if integrator.opts.adaptive
        utilde = dt * (btilde1 * k1 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5)
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.u = u
end

get_fsalfirstlast(cache::OwrenZen4Cache, u) = (cache.k1, cache.k6)
function initialize!(integrator, cache::OwrenZen4Cache)
    integrator.kshortsize = 6
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::OwrenZen4Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k1, k2, k3, k4, k5, k6, utilde, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, c1, c2, c3, c4, btilde1, btilde3, btilde4, btilde5) = cache.tab
    a = dt * a21
    @.. broadcast = false thread = thread tmp = uprev + a * k1
    stage_limiter!(tmp, integrator, p, t + c1 * dt)
    f(k2, tmp, p, t + c1 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k3, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k4, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k5, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k6, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = dt * (
            btilde1 * k1 + btilde3 * k3 +
                btilde4 * k4 +
                btilde5 * k5
        )
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    return nothing
end

function initialize!(integrator, cache::OwrenZen5ConstantCache)
    integrator.kshortsize = 8
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(
        integrator, cache::OwrenZen5ConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; a21, a31, a32, a41, a42, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, c1, c2, c3, c4, c5, c6, btilde1, btilde3, btilde4, btilde5, btilde6, btilde7) = cache
    k1 = integrator.fsalfirst
    a = dt * a21
    k2 = f(uprev + a * k1, p, t + c1 * dt)
    k3 = f(uprev + dt * (a31 * k1 + a32 * k2), p, t + c2 * dt)
    k4 = f(uprev + dt * (a41 * k1 + a42 * k2 + k3), p, t + c3 * dt)
    k5 = f(uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4), p, t + c4 * dt)
    k6 = f(
        uprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5), p,
        t + c5 * dt
    )
    k7 = f(
        uprev + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6),
        p, t + c6 * dt
    )
    u = uprev + dt * (a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 + a86 * k6 + a87 * k7)
    k8 = f(u, p, t + dt)
    integrator.fsallast = k8
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 7)
    if integrator.opts.adaptive
        utilde = dt *
            (
            btilde1 * k1 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5 + btilde6 * k6 +
                btilde7 * k7
        )
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.k[7] = k7
    integrator.k[8] = k8
    integrator.u = u
end

get_fsalfirstlast(cache::OwrenZen5Cache, u) = (cache.k1, cache.k8)
function initialize!(integrator, cache::OwrenZen5Cache)
    integrator.kshortsize = 8
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6
    integrator.k[7] = cache.k7
    integrator.k[8] = cache.k8
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::OwrenZen5Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k1, k2, k3, k4, k5, k6, k7, k8, utilde, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a31, a32, a41, a42, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, c1, c2, c3, c4, c5, c6, btilde1, btilde3, btilde4, btilde5, btilde6, btilde7) = cache.tab
    a = dt * a21
    @.. broadcast = false thread = thread tmp = uprev + a * k1
    stage_limiter!(tmp, integrator, p, t + c1 * dt)
    f(k2, tmp, p, t + c1 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k3, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a42 * k2 + k3)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k4, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k5, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 +
            a65 * k5
    )
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k6, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 +
            a75 * k5 +
            a76 * k6
    )
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k7, tmp, p, t + c6 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (
        a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 +
            a86 * k6 + a87 * k7
    )
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k8, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 7)
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = dt * (
            btilde1 * k1 + btilde3 * k3 +
                btilde4 * k4 +
                btilde5 * k5 + btilde6 * k6 +
                btilde7 * k7
        )
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    return nothing
end

function initialize!(integrator, cache::BS5ConstantCache)
    alg = unwrap_alg(integrator, false)
    alg.lazy ? (integrator.kshortsize = 8) : (integrator.kshortsize = 11)
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:7
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    integrator.k[integrator.kshortsize] = integrator.fsallast

    return if !alg.lazy
        @inbounds for i in 9:11
            integrator.k[i] = zero(integrator.fsalfirst)
        end
    end
end

@muladd function perform_step!(integrator, cache::BS5ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, bhat1, bhat3, bhat4, bhat5, bhat6, btilde1, btilde3, btilde4, btilde5, btilde6, btilde7, btilde8) = cache
    k1 = integrator.fsalfirst
    a = dt * a21
    k2 = f(uprev + a * k1, p, t + c1 * dt)
    k3 = f(uprev + dt * (a31 * k1 + a32 * k2), p, t + c2 * dt)
    k4 = f(uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3), p, t + c3 * dt)
    k5 = f(uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4), p, t + c4 * dt)
    k6 = f(
        uprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5), p,
        t + c5 * dt
    )
    k7 = f(
        uprev + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6),
        p, t + dt
    )
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    u = uprev + dt * (a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 + a86 * k6 + a87 * k7)
    integrator.fsallast = f(u, p, t + dt)
    k8 = integrator.fsallast
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    if integrator.opts.adaptive
        uhat = dt * (bhat1 * k1 + bhat3 * k3 + bhat4 * k4 + bhat5 * k5 + bhat6 * k6)
        utilde = dt *
            (
            btilde1 * k1 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5 + btilde6 * k6 +
                btilde7 * k7 + btilde8 * k8
        )
        atmp = calculate_residuals(
            uhat, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        EEst1 = integrator.opts.internalnorm(atmp, t)
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        EEst2 = integrator.opts.internalnorm(atmp, t)
        integrator.EEst = max(EEst1, EEst2)
    end
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.k[7] = k7
    integrator.k[8] = k8
    integrator.u = u

    alg = unwrap_alg(integrator, false)
    if !alg.lazy && (
            integrator.opts.adaptive == false ||
                accept_step_controller(integrator, integrator.opts.controller)
        )
        (; c6, c7, c8, a91, a92, a93, a94, a95, a96, a97, a98, a101, a102, a103, a104, a105, a106, a107, a108, a109, a111, a112, a113, a114, a115, a116, a117, a118, a119, a1110) = cache
        k = integrator.k
        k[9] = f(
            uprev +
                dt * (
                a91 * k[1] + a92 * k[2] + a93 * k[3] + a94 * k[4] + a95 * k[5] +
                    a96 * k[6] + a97 * k[7] + a98 * k[8]
            ),
            p,
            t + c6 * dt
        )
        k[10] = f(
            uprev +
                dt *
                (
                a101 * k[1] + a102 * k[2] + a103 * k[3] + a104 * k[4] + a105 * k[5] +
                    a106 * k[6] + a107 * k[7] + a108 * k[8] + a109 * k[9]
            ),
            p,
            t + c7 * dt
        )
        k[11] = f(
            uprev +
                dt *
                (
                a111 * k[1] + a112 * k[2] + a113 * k[3] + a114 * k[4] + a115 * k[5] +
                    a116 * k[6] + a117 * k[7] + a118 * k[8] + a119 * k[9] + a1110 * k[10]
            ),
            p, t + c8 * dt
        )
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    end
end

get_fsalfirstlast(cache::BS5Cache, u) = (cache.k1, cache.k8)
function initialize!(integrator, cache::BS5Cache)
    alg = unwrap_alg(integrator, false)
    alg.lazy ? (integrator.kshortsize = 8) : (integrator.kshortsize = 11)
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6
    integrator.k[7] = cache.k7
    integrator.k[8] = cache.k8

    if !alg.lazy
        integrator.k[9] = similar(cache.k1)
        integrator.k[10] = similar(cache.k1)
        integrator.k[11] = similar(cache.k1)
    end
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::BS5Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; k1, k2, k3, k4, k5, k6, k7, k8, utilde, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, bhat1, bhat3, bhat4, bhat5, bhat6, btilde1, btilde3, btilde4, btilde5, btilde6, btilde7, btilde8) = cache.tab
    a = dt * a21
    @.. broadcast = false thread = thread tmp = uprev + a * k1
    stage_limiter!(tmp, integrator, p, t + c1 * dt)
    f(k2, tmp, p, t + c1 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k3, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k4, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k5, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 +
            a65 * k5
    )
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k6, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 +
            a75 * k5 +
            a76 * k6
    )
    stage_limiter!(tmp, integrator, p, t + dt)
    f(k7, tmp, p, t + dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (
        a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 +
            a86 * k6 + a87 * k7
    )
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k8, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 7)
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = dt *
            (
            bhat1 * k1 + bhat3 * k3 + bhat4 * k4 +
                bhat5 * k5 +
                bhat6 * k6
        )
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        EEst1 = integrator.opts.internalnorm(atmp, t)
        @.. broadcast = false thread = thread utilde = dt * (
            btilde1 * k1 + btilde3 * k3 +
                btilde4 * k4 +
                btilde5 * k5 + btilde6 * k6 +
                btilde7 * k7 +
                btilde8 * k8
        )
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        EEst2 = integrator.opts.internalnorm(atmp, t)
        integrator.EEst = max(EEst1, EEst2)
    end
    alg = unwrap_alg(integrator, false)
    if !alg.lazy && (
            integrator.opts.adaptive == false ||
                accept_step_controller(integrator, integrator.opts.controller)
        )
        k = integrator.k
        (; c6, c7, c8, a91, a92, a93, a94, a95, a96, a97, a98, a101, a102, a103, a104, a105, a106, a107, a108, a109, a111, a112, a113, a114, a115, a116, a117, a118, a119, a1110) = cache.tab
        @.. broadcast = false thread = thread tmp = uprev +
            dt * (
            a91 * k[1] + a92 * k[2] + a93 * k[3] +
                a94 * k[4] +
                a95 * k[5] + a96 * k[6] + a97 * k[7] +
                a98 * k[8]
        )
        f(k[9], tmp, p, t + c6 * dt)
        @.. broadcast = false thread = thread tmp = uprev +
            dt *
            (
            a101 * k[1] + a102 * k[2] + a103 * k[3] +
                a104 * k[4] +
                a105 * k[5] + a106 * k[6] + a107 * k[7] +
                a108 * k[8] +
                a109 * k[9]
        )
        f(k[10], tmp, p, t + c7 * dt)
        @.. broadcast = false thread = thread tmp = uprev +
            dt *
            (
            a111 * k[1] + a112 * k[2] + a113 * k[3] +
                a114 * k[4] +
                a115 * k[5] + a116 * k[6] + a117 * k[7] +
                a118 * k[8] +
                a119 * k[9] + a1110 * k[10]
        )
        f(k[11], tmp, p, t + c8 * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    end
    return nothing
end

function initialize!(integrator, cache::DP5ConstantCache)
    integrator.kshortsize = 4
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    return @inbounds for i in eachindex(integrator.k)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
end

@muladd function perform_step!(integrator, cache::DP5ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract DP5ConstantCacheActual T T2
    k1 = integrator.fsalfirst
    a = dt * a21
    k2 = f(uprev + a * k1, p, t + c1 * dt)
    k3 = f(uprev + dt * (a31 * k1 + a32 * k2), p, t + c2 * dt)
    k4 = f(uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3), p, t + c3 * dt)
    k5 = f(uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4), p, t + c4 * dt)
    g6 = uprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = f(g6, p, t + dt)
    update = a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6
    u = uprev + dt * update
    integrator.fsallast = f(u, p, t + dt)
    k7 = integrator.fsallast
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    if integrator.alg isa CompositeAlgorithm
        g7 = u
        # Hairer II, page 22 modified to use the Inf norm
        integrator.eigen_est = integrator.opts.internalnorm(
            maximum(abs.((k7 .- k6) ./ (g7 .- g6))), t
        )
    end
    if integrator.opts.adaptive
        utilde = dt *
            (
            btilde1 * k1 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5 + btilde6 * k6 +
                btilde7 * k7
        )
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.k[1] = update
    bspl = k1 - update
    integrator.k[2] = bspl
    integrator.k[3] = update - k7 - bspl
    integrator.k[4] = d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6 + d7 * k7
    integrator.u = u
end

get_fsalfirstlast(cache::DP5Cache, u) = (cache.k1, cache.k7)
function initialize!(integrator, cache::DP5Cache)
    integrator.kshortsize = 4
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.update
    integrator.k[2] = cache.bspl
    integrator.k[3] = cache.dense_tmp3
    integrator.k[4] = cache.dense_tmp4
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::DP5Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract DP5ConstantCacheActual T T2
    (; k1, k2, k3, k4, k5, k6, k7, dense_tmp3, dense_tmp4, update, bspl, utilde, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    a = dt * a21
    @.. broadcast = false thread = thread tmp = uprev + a * k1
    stage_limiter!(tmp, integrator, p, t + c1 * dt)
    f(k2, tmp, p, t + c1 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k3, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k4, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k5, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 +
            a65 * k5
    )
    stage_limiter!(tmp, integrator, p, t + dt)
    f(k6, tmp, p, t + dt)
    @.. broadcast = false thread = thread update = a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 +
        a76 * k6
    @.. broadcast = false thread = thread u = uprev + dt * update
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k7, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    if integrator.alg isa CompositeAlgorithm
        g6 = tmp
        g7 = u
        # Hairer II, page 22 modified to use Inf norm
        @.. broadcast = false thread = thread utilde = abs((k7 - k6) / (g7 - g6))
        integrator.eigen_est = integrator.opts.internalnorm(
            norm(utilde, Inf) * oneunit(t), t
        )
    end
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = dt * (
            btilde1 * k1 + btilde3 * k3 +
                btilde4 * k4 +
                btilde5 * k5 + btilde6 * k6 +
                btilde7 * k7
        )
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    if integrator.opts.calck
        #integrator.k[4] == k5
        @.. broadcast = false thread = thread integrator.k[4] = d1 * k1 + d3 * k3 + d4 * k4 +
            d5 * k5 +
            d6 * k6 + d7 * k7
        #bspl == k3
        @.. broadcast = false thread = thread bspl = k1 - update
        # k6 === integrator.k[3] === k2
        @.. broadcast = false thread = thread integrator.k[3] = update - k7 - bspl
    end
    return nothing
end

function initialize!(integrator, cache::RKO65ConstantCache)
    integrator.kshortsize = 6
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::RKO65ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α21, α31, α41, α51, α32, α42, α52, α62, α43, α53, α63, α54, α64, α65, β2, β3, β4, β5, β6, c1, c2, c3, c4, c5, c6) = cache

    #k1=integrator.fsalfirst #f(uprev,p,t)
    k1 = f(uprev, p, t + c1 * dt)
    k2 = f(uprev + α21 * dt * k1, p, t + c2 * dt)
    k3 = f(uprev + α31 * dt * k1 + α32 * dt * k2, p, t + c3 * dt)
    k4 = f(uprev + α41 * dt * k1 + α42 * dt * k2 + α43 * dt * k3, p, t + c4 * dt)
    k5 = f(
        uprev + α51 * dt * k1 + α52 * dt * k2 + α53 * dt * k3 + α54 * dt * k4, p,
        t + c5 * dt
    )
    k6 = f(
        uprev + α62 * dt * k2 + α63 * dt * k3 + α64 * dt * k4 + α65 * dt * k5, p,
        t + c6 * dt
    )
    u = uprev + dt * (β2 * k2 + β3 * k3 + β4 * k4 + β5 * k5 + β6 * k6)

    integrator.fsallast = f(u, p, t + dt)  # For interpolation, then FSAL'd

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.u = u
end

get_fsalfirstlast(cache::RKO65Cache, u) = (cache.k1, cache.k6)
function initialize!(integrator, cache::RKO65Cache)
    (; k, fsalfirst) = cache
    integrator.kshortsize = 6
    resize!(integrator.k, integrator.kshortsize)

    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6

    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::RKO65Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, k, k1, k2, k3, k4, k5, k6, stage_limiter!, step_limiter!, thread) = cache
    (; α21, α31, α41, α51, α32, α42, α52, α62, α43, α53, α63, α54, α64, α65, β2, β3, β4, β5, β6, c1, c2, c3, c4, c5, c6) = cache.tab
    #println("L221: tmp", tmp)
    f(k1, uprev, p, t + c1 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α21 * dt * k1
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    #println("L224: tmp/k", tmp, k1)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α31 * dt * k1 + α32 * dt * k2
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α41 * dt * k1 + α42 * dt * k2 +
        α43 * dt * k3
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α51 * dt * k1 + α52 * dt * k2 +
        α53 * dt * k3 +
        α54 * dt * k4
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α62 * dt * k2 + α63 * dt * k3 +
        α64 * dt * k4 +
        α65 * dt * k5
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k6, tmp, p, t + c6 * dt)

    @.. broadcast = false thread = thread u = uprev +
        dt *
        (β2 * k2 + β3 * k3 + β4 * k4 + β5 * k5 + β6 * k6)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    #println("L238: tmp/u", tmp, u)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)

    #return nothing
end

function initialize!(integrator, cache::FRK65ConstantCache)
    integrator.kshortsize = 9
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::FRK65ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α21, α31, α41, α51, α61, α71, α81, α91, α32, α43, α53, α63, α73, α83, α54, α64, α74, α84, α94, α65, α75, α85, α95, α76, α86, α96, α87, α97, α98, β1, β7, β8, β1tilde, β4tilde, β5tilde, β6tilde, β7tilde, β8tilde, β9tilde, c2, c3, c4, c5, c6, c7, c8, c9, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11) = cache
    alg = unwrap_alg(integrator, false)
    ν = alg.omega * dt
    νsq = ν^2
    β4 = (d1 + νsq * (d2 + νsq * (d3 + νsq * (d4 + νsq * (d5 + νsq * (d6 + +νsq * d7)))))) /
        (
        1 +
            νsq * (d8 + νsq * (d9 + νsq * (d10 + νsq * (d11 + νsq * (d12 + +νsq * d13)))))
    )
    β5 = (e1 + νsq * (e2 + νsq * (e3 + νsq * (e4 + νsq * (e5 + νsq * e6))))) /
        (1 + νsq * (e8 + νsq * (e9 + νsq * (e10 + νsq * e11))))
    β6 = (f1 + νsq * (f2 + νsq * (f3 + νsq * (f4 + νsq * (f5 + νsq * f6))))) /
        (1 + νsq * (f8 + νsq * (f9 + νsq * (f10 + νsq * f11))))

    k1 = integrator.fsalfirst
    k2 = f(uprev + α21 * dt * k1, p, t + c2 * dt)
    k3 = f(uprev + α31 * dt * k1 + α32 * dt * k2, p, t + c3 * dt)
    k4 = f(uprev + α41 * dt * k1 + α43 * dt * k3, p, t + c4 * dt)
    k5 = f(uprev + α51 * dt * k1 + α53 * dt * k3 + α54 * dt * k4, p, t + c5 * dt)
    k6 = f(
        uprev + α61 * dt * k1 + α63 * dt * k3 + α64 * dt * k4 + α65 * dt * k5, p,
        t + c6 * dt
    )
    k7 = f(
        uprev + α71 * dt * k1 + α73 * dt * k3 + α74 * dt * k4 + α75 * dt * k5 +
            α76 * dt * k6,
        p,
        t + c7 * dt
    )
    k8 = f(
        uprev + α81 * dt * k1 + α83 * dt * k3 + α84 * dt * k4 + α85 * dt * k5 +
            α86 * dt * k6 + α87 * dt * k7,
        p,
        t + c8 * dt
    )
    u = uprev + dt * (β1 * k1 + β4 * k4 + β5 * k5 + β6 * k6 + β7 * k7 + β8 * k8)
    integrator.fsallast = f(u, p, t + dt)
    k9 = integrator.fsallast
    if integrator.opts.adaptive
        utilde = dt *
            (
            β1tilde * k1 + β4tilde * k4 + β5tilde * k5 + β6tilde * k6 + β7tilde * k7 +
                β8tilde * k8 + β9tilde * k9
        )
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 8)
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.k[7] = k7
    integrator.k[8] = k8
    integrator.k[9] = k9
    integrator.u = u
end

get_fsalfirstlast(cache::FRK65Cache, u) = (cache.k1, cache.k9)
function initialize!(integrator, cache::FRK65Cache)
    integrator.kshortsize = 9

    resize!(integrator.k, integrator.kshortsize)

    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6
    integrator.k[7] = cache.k7
    integrator.k[8] = cache.k8
    integrator.k[9] = cache.k9

    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::FRK65Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, k1, k2, k3, k4, k5, k6, k7, k8, k9, utilde, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; α21, α31, α41, α51, α61, α71, α81, α91, α32, α43, α53, α63, α73, α83, α54, α64, α74, α84, α94, α65, α75, α85, α95, α76, α86, α96, α87, α97, α98, β1, β7, β8, β1tilde, β4tilde, β5tilde, β6tilde, β7tilde, β8tilde, β9tilde, c2, c3, c4, c5, c6, c7, c8, c9, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11) = cache.tab
    alg = unwrap_alg(integrator, false)

    ν = alg.omega * dt
    νsq = ν^2
    β4 = (d1 + νsq * (d2 + νsq * (d3 + νsq * (d4 + νsq * (d5 + νsq * (d6 + +νsq * d7)))))) /
        (
        1 +
            νsq * (d8 + νsq * (d9 + νsq * (d10 + νsq * (d11 + νsq * (d12 + +νsq * d13)))))
    )
    β5 = (e1 + νsq * (e2 + νsq * (e3 + νsq * (e4 + νsq * (e5 + νsq * e6))))) /
        (1 + νsq * (e8 + νsq * (e9 + νsq * (e10 + νsq * e11))))
    β6 = (f1 + νsq * (f2 + νsq * (f3 + νsq * (f4 + νsq * (f5 + νsq * f6))))) /
        (1 + νsq * (f8 + νsq * (f9 + νsq * (f10 + νsq * f11))))

    @.. broadcast = false thread = thread tmp = uprev + α21 * dt * k1
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α31 * dt * k1 + α32 * dt * k2
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α41 * dt * k1 + α43 * dt * k3
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α51 * dt * k1 + α53 * dt * k3 +
        α54 * dt * k4
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α61 * dt * k1 + α63 * dt * k3 +
        α64 * dt * k4 +
        α65 * dt * k5
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k6, tmp, p, t + c6 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α71 * dt * k1 + α73 * dt * k3 +
        α74 * dt * k4 +
        α75 * dt * k5 + α76 * dt * k6
    stage_limiter!(tmp, integrator, p, t + c7 * dt)
    f(k7, tmp, p, t + c7 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α81 * dt * k1 + α83 * dt * k3 +
        α84 * dt * k4 +
        α85 * dt * k5 + α86 * dt * k6 + α87 * dt * k7
    stage_limiter!(tmp, integrator, p, t + c8 * dt)
    f(k8, tmp, p, t + c8 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (
        β1 * k1 + β4 * k4 + β5 * k5 + β6 * k6 + β7 * k7 +
            β8 * k8
    )
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k9, u, p, t + dt)

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 8)
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = dt * (
            β1tilde * k1 + β4tilde * k4 +
                β5tilde * k5 +
                β6tilde * k6 + β7tilde * k7 +
                β8tilde * k8 +
                β9tilde * k9
        )
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    return nothing
end

function initialize!(integrator, cache::RKMConstantCache)
    integrator.kshortsize = 6
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::RKMConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; α2, α3, α4, α5, α6, β1, β2, β3, β4, β6, c2, c3, c4, c5, c6) = cache

    #k1 = f(uprev, p, t)
    k1 = integrator.fsalfirst
    k2 = f(uprev + α2 * dt * k1, p, t + c2 * dt)
    k3 = f(uprev + α3 * dt * k2, p, t + c3 * dt)
    k4 = f(uprev + α4 * dt * k3, p, t + c4 * dt)
    k5 = f(uprev + α5 * dt * k4, p, t + c5 * dt)
    k6 = f(uprev + α6 * dt * k5, p, t + c6 * dt)
    u = uprev + dt * (β1 * k1 + β2 * k2 + β3 * k3 + β4 * k4 + β6 * k6)

    integrator.fsallast = f(u, p, t + dt) #interpolation then FSAL'd
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    # integrator.k[1] = integrator.fsalfirst
    # integrator.k[2] = integrator.fsallast
    integrator.u = u
end

get_fsalfirstlast(cache::RKMCache, u) = (cache.k1, zero(cache.k1))
function initialize!(integrator, cache::RKMCache)
    (; k, fsalfirst) = cache
    integrator.kshortsize = 6
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6

    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::RKMCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, fsalfirst, k, k1, k2, k3, k4, k5, k6, stage_limiter!, step_limiter!, thread) = cache
    (; α2, α3, α4, α5, α6, β1, β2, β3, β4, β6, c2, c3, c4, c5, c6) = cache.tab

    @.. broadcast = false thread = thread tmp = uprev + α2 * dt * k1
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α3 * dt * k2
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α4 * dt * k3
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α5 * dt * k4
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = uprev + α6 * dt * k5
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k6, tmp, p, t + c6 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (β1 * k1 + β2 * k2 + β3 * k3 + β4 * k4 + β6 * k6)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    return nothing
end

function initialize!(integrator, cache::PSRK4p7q6ConstantCache)
    integrator.kshortsize = 6
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = zero(integrator.fsalfirst)
    integrator.k[3] = zero(integrator.fsalfirst)
    integrator.k[4] = zero(integrator.fsalfirst)
    integrator.k[5] = zero(integrator.fsalfirst)
    return integrator.k[6] = integrator.fsallast
end

function perform_step!(integrator, cache::PSRK4p7q6ConstantCache, repeat_step = false)
    (; u, uprev, f, p, dt, t) = integrator
    (; a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, b1, b2, b3, b4, b5, b6, c2, c3, c4, c5, c6) = cache

    k1 = f(uprev, p, t)
    tmp = uprev + dt * (a21 * k1)
    k2 = f(tmp, p, t + c2 * dt)
    tmp = uprev + dt * (a31 * k1 + a32 * k2)
    k3 = f(tmp, p, t + c3 * dt)
    tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    k4 = f(tmp, p, t + dt * c4)
    tmp = uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    k5 = f(tmp, p, t + dt * c5)
    tmp = uprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = f(tmp, p, t + dt * c6)
    u = uprev + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6)

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    integrator.fsallast = k6

    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    return integrator.u = u
end

get_fsalfirstlast(cache::PSRK4p7q6Cache, u) = (cache.k1, cache.k6)
function initialize!(integrator, cache::PSRK4p7q6Cache)
    (; uprev, f, p, t) = integrator

    integrator.kshortsize = 6
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    return integrator.k[6] = cache.k6
end

function perform_step!(integrator, cache::PSRK4p7q6Cache, repeat_step = false)
    (; k1, k2, k3, k4, k5, k6, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, b1, b2, b3, b4, b5, b6, c2, c3, c4, c5, c6) = cache.tab
    (; u, uprev, t, dt, f, p) = integrator

    f(k1, uprev, p, t)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a21 * k1)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 +
            a65 * k5
    )
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k6, tmp, p, t + c6 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (
        b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 +
            b6 * k6
    )
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    integrator.fsallast = k6
    return nothing
end

function initialize!(integrator, cache::PSRK3p6q5ConstantCache)
    integrator.kshortsize = 5
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = zero(integrator.fsalfirst)
    integrator.k[3] = zero(integrator.fsalfirst)
    integrator.k[4] = zero(integrator.fsalfirst)
    return integrator.k[5] = integrator.fsallast
end

function perform_step!(integrator, cache::PSRK3p6q5ConstantCache, repeat_step = false)
    (; u, uprev, f, p, dt, t) = integrator
    (; a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, b1, b2, b3, b4, b5, c2, c3, c4, c5) = cache

    k1 = f(uprev, p, t)
    tmp = uprev + dt * (a21 * k1)
    k2 = f(tmp, p, t + c2 * dt)
    tmp = uprev + dt * (a31 * k1 + a32 * k2)
    k3 = f(tmp, p, t + c3 * dt)
    tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    k4 = f(tmp, p, t + dt * c4)
    tmp = uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    k5 = f(tmp, p, t + dt * c5)
    u = uprev + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5)

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
    integrator.fsallast = k5

    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    return integrator.u = u
end

get_fsalfirstlast(cache::PSRK3p6q5Cache, u) = (cache.k1, cache.k5)
function initialize!(integrator, cache::PSRK3p6q5Cache)
    (; uprev, f, p, t) = integrator

    integrator.kshortsize = 5
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    return integrator.k[5] = cache.k5
end

function perform_step!(integrator, cache::PSRK3p6q5Cache, repeat_step = false)
    (; k1, k2, k3, k4, k5, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, b1, b2, b3, b4, b5, c2, c3, c4, c5) = cache.tab
    (; u, uprev, t, dt, f, p) = integrator

    f(k1, uprev, p, t)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a21 * k1)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 5)
    integrator.fsallast = k5
    return nothing
end

function initialize!(integrator, cache::PSRK3p5q4ConstantCache)
    integrator.kshortsize = 4
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = zero(integrator.fsalfirst)
    integrator.k[3] = zero(integrator.fsalfirst)
    return integrator.k[4] = integrator.fsallast
end

function perform_step!(integrator, cache::PSRK3p5q4ConstantCache, repeat_step = false)
    (; u, uprev, f, p, dt, t) = integrator
    (; a21, a31, a32, a41, a42, a43, b1, b2, b3, b4, c2, c3, c4) = cache

    k1 = f(uprev, p, t)
    tmp = uprev + dt * (a21 * k1)
    k2 = f(tmp, p, t + c2 * dt)
    tmp = uprev + dt * (a31 * k1 + a32 * k2)
    k3 = f(tmp, p, t + c3 * dt)
    tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    k4 = f(tmp, p, t + dt * c4)
    u = uprev + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4)

    integrator.fsallast = k4
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)

    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    return integrator.u = u
end

get_fsalfirstlast(cache::PSRK3p5q4Cache, u) = (cache.k1, cache.k4)
function initialize!(integrator, cache::PSRK3p5q4Cache)
    (; uprev, f, p, t) = integrator

    integrator.kshortsize = 4
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.fsalfirst = cache.k1
    return integrator.fsallast = cache.k4
end

function perform_step!(integrator, cache::PSRK3p5q4Cache, repeat_step = false)
    (; k1, k2, k3, k4, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a31, a32, a41, a42, a43, b1, b2, b3, b4, c2, c3, c4) = cache.tab
    (; u, uprev, t, dt, f, p) = integrator

    f(k1, uprev, p, t)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a21 * k1)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
    integrator.fsallast = k4
    return nothing
end

function initialize!(integrator, cache::MSRK5ConstantCache)
    integrator.kshortsize = 9
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.fsallast = zero(integrator.fsalfirst)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

function perform_step!(integrator, cache::MSRK5ConstantCache, repeat_step = false)
    (; u, uprev, f, p, dt, t) = integrator
    (; a21, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, b1, b4, b5, b6, b7, b8, c2, c3, c4, c5, c6, c7, c8) = cache

    k1 = integrator.fsalfirst
    tmp = uprev + dt * (a21 * k1)
    k2 = f(tmp, p, t + c2 * dt)
    tmp = uprev + dt * (a31 * k1 + a32 * k2)
    k3 = f(tmp, p, t + c3 * dt)
    tmp = uprev + dt * (a41 * k1 + a43 * k3)
    k4 = f(tmp, p, t + dt * c4)
    tmp = uprev + dt * (a51 * k1 + a53 * k3 + a54 * k4)
    k5 = f(tmp, p, t + dt * c5)
    tmp = uprev + dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = f(tmp, p, t + dt * c6)
    tmp = uprev + dt * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
    k7 = f(tmp, p, t + dt * c7)
    tmp = uprev + dt * (a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 + a86 * k6 + a87 * k7)
    k8 = f(tmp, p, t + dt * c8)
    u = uprev + dt * (b1 * k1 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7 + b8 * k8)
    k9 = f(u, p, t + dt)
    integrator.fsallast = k9
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 8)

    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.k[7] = k7
    integrator.k[8] = k8
    integrator.k[9] = k9
    return integrator.u = u
end

get_fsalfirstlast(cache::MSRK5Cache, u) = (cache.k1, cache.k9)
function initialize!(integrator, cache::MSRK5Cache)
    (; uprev, f, p, t) = integrator

    integrator.kshortsize = 9
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6
    integrator.k[7] = cache.k7
    integrator.k[8] = cache.k8
    integrator.k[9] = cache.k9
    integrator.fsalfirst = cache.k1
    integrator.fsallast = cache.k9

    f(integrator.fsalfirst, uprev, p, t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::MSRK5Cache, repeat_step = false)
    (; k1, k2, k3, k4, k5, k6, k7, k8, k9, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, b1, b4, b5, b6, b7, b8, c2, c3, c4, c5, c6, c7, c8) = cache.tab
    (; u, uprev, t, dt, f, p) = integrator

    @.. broadcast = false thread = thread tmp = uprev + dt * (a21 * k1)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a51 * k1 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k6, tmp, p, t + c6 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 +
            a76 * k6
    )
    stage_limiter!(tmp, integrator, p, t + c7 * dt)
    f(k7, tmp, p, t + c7 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 +
            a86 * k6 +
            a87 * k7
    )

    stage_limiter!(tmp, integrator, p, t + c8 * dt)
    f(k8, tmp, p, t + c8 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (
        b1 * k1 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7 +
            b8 * k8
    )
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k9, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 8)
    integrator.fsallast = k9

    return nothing
end

function initialize!(integrator, cache::MSRK6ConstantCache)
    integrator.kshortsize = 9
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.fsallast = zero(integrator.fsalfirst)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

function perform_step!(integrator, cache::MSRK6ConstantCache, repeat_step = false)
    (; u, uprev, f, p, dt, t) = integrator
    (; a21, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, b1, b4, b5, b6, b7, b8, c2, c3, c4, c5, c6, c7, c8) = cache

    k1 = integrator.fsalfirst
    tmp = uprev + dt * (a21 * k1)
    k2 = f(tmp, p, t + c2 * dt)
    tmp = uprev + dt * (a32 * k2)
    k3 = f(tmp, p, t + c3 * dt)
    tmp = uprev + dt * (a41 * k1 + a43 * k3)
    k4 = f(tmp, p, t + dt * c4)
    tmp = uprev + dt * (a51 * k1 + a53 * k3 + a54 * k4)
    k5 = f(tmp, p, t + dt * c5)
    tmp = uprev + dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = f(tmp, p, t + dt * c6)
    tmp = uprev + dt * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
    k7 = f(tmp, p, t + dt * c7)
    tmp = uprev + dt * (a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 + a86 * k6 + a87 * k7)
    k8 = f(tmp, p, t + dt * c8)
    u = uprev + dt * (b1 * k1 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7 + b8 * k8)
    k9 = f(u, p, t + dt)
    integrator.fsallast = k9
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 8)

    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.k[7] = k7
    integrator.k[8] = k8
    integrator.k[9] = k9
    return integrator.u = u
end

get_fsalfirstlast(cache::MSRK6Cache, u) = (cache.k1, cache.k9)
function initialize!(integrator, cache::MSRK6Cache)
    (; uprev, f, p, t) = integrator

    integrator.kshortsize = 9
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6
    integrator.k[7] = cache.k7
    integrator.k[8] = cache.k8
    integrator.k[9] = cache.k9
    integrator.fsalfirst = cache.k1
    integrator.fsallast = cache.k9

    f(integrator.fsalfirst, uprev, p, t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::MSRK6Cache, repeat_step = false)
    (; k1, k2, k3, k4, k5, k6, k7, k8, k9, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, b1, b4, b5, b6, b7, b8, c2, c3, c4, c5, c6, c7, c8) = cache.tab
    (; u, uprev, t, dt, f, p) = integrator

    @.. broadcast = false thread = thread tmp = uprev + dt * (a21 * k1)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a51 * k1 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k6, tmp, p, t + c6 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 +
            a76 * k6
    )
    stage_limiter!(tmp, integrator, p, t + c7 * dt)
    f(k7, tmp, p, t + c7 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 +
            a86 * k6 +
            a87 * k7
    )
    stage_limiter!(tmp, integrator, p, t + c8 * dt)
    f(k8, tmp, p, t + c8 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (
        b1 * k1 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7 +
            b8 * k8
    )
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k9, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 8)
    integrator.fsallast = k9

    return nothing
end

function initialize!(integrator, cache::Stepanov5ConstantCache)
    integrator.kshortsize = 7
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.fsallast = zero(integrator.fsalfirst)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

function perform_step!(integrator, cache::Stepanov5ConstantCache, repeat_step = false)
    (; u, uprev, f, p, dt, t) = integrator
    (; a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, b1, b3, b4, b5, b6, btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7, c2, c3, c4, c5, c6) = cache

    k1 = integrator.fsalfirst
    tmp = uprev + dt * (a21 * k1)
    k2 = f(tmp, p, t + c2 * dt)
    tmp = uprev + dt * (a31 * k1 + a32 * k2)
    k3 = f(tmp, p, t + c3 * dt)
    tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    k4 = f(tmp, p, t + dt * c4)
    tmp = uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    k5 = f(tmp, p, t + dt * c5)
    tmp = uprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = f(tmp, p, t + dt * c6)
    u = uprev + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6)
    k7 = f(u, p, t + dt)
    integrator.fsallast = k7
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)

    if integrator.opts.adaptive
        utilde = dt * (
            btilde1 * k1 + btilde2 * k2 +
                btilde3 * k3 + btilde4 * k4 +
                btilde5 * k5 + btilde6 * k6 +
                btilde7 * k7
        )
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.k[7] = k7
    return integrator.u = u
end

get_fsalfirstlast(cache::Stepanov5Cache, u) = (cache.k1, cache.k7)
function initialize!(integrator, cache::Stepanov5Cache)
    (; uprev, f, p, t) = integrator

    integrator.kshortsize = 7
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6
    integrator.k[7] = cache.k7
    integrator.fsalfirst = cache.k1
    integrator.fsallast = cache.k7

    f(integrator.fsalfirst, uprev, p, t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::Stepanov5Cache, repeat_step = false)
    (; k1, k2, k3, k4, k5, k6, k7, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, b1, b3, b4, b5, b6, btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7, c2, c3, c4, c5, c6) = cache.tab
    (; u, uprev, t, dt, f, p) = integrator

    @.. broadcast = false thread = thread tmp = uprev + dt * (a21 * k1)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 +
            a65 * k5
    )
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k6, tmp, p, t + c6 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k7, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    integrator.fsallast = k7

    if integrator.opts.adaptive
        utilde = dt *
            (
            btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5 +
                btilde6 * k6 +
                btilde7 * k7
        )
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    return nothing
end

function initialize!(integrator, cache::SIR54ConstantCache)
    integrator.kshortsize = 8
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.fsallast = zero(integrator.fsalfirst)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

function perform_step!(integrator, cache::SIR54ConstantCache, repeat_step = false)
    (; u, uprev, f, p, dt, t) = integrator
    (; a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, b1, b2, b3, b4, b5, b6, btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7, c2, c3, c4, c5, c6, c7) = cache

    k1 = integrator.fsalfirst
    tmp = uprev + dt * (a21 * k1)
    k2 = f(tmp, p, t + c2 * dt)
    tmp = uprev + dt * (a31 * k1 + a32 * k2)
    k3 = f(tmp, p, t + c3 * dt)
    tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    k4 = f(tmp, p, t + dt * c4)
    tmp = uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    k5 = f(tmp, p, t + dt * c5)
    tmp = uprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = f(tmp, p, t + dt * c6)
    tmp = uprev + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
    k7 = f(tmp, p, t + dt * c7)
    u = uprev + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6)
    k8 = f(u, p, t + dt)
    integrator.fsallast = k8
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 7)

    if integrator.opts.adaptive
        utilde = dt * (
            btilde1 * k1 + btilde2 * k2 + btilde3 * k3 +
                btilde4 * k4 +
                btilde5 * k5 + btilde6 * k6 +
                btilde7 * k7
        )
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.k[7] = k7
    integrator.k[8] = k8
    return integrator.u = u
end

get_fsalfirstlast(cache::SIR54Cache, u) = (cache.k1, cache.k8)
function initialize!(integrator, cache::SIR54Cache)
    (; uprev, f, p, t) = integrator

    integrator.kshortsize = 8
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6
    integrator.k[7] = cache.k7
    integrator.k[8] = cache.k8
    integrator.fsalfirst = cache.k1
    integrator.fsallast = cache.k8

    f(integrator.fsalfirst, uprev, p, t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::SIR54Cache, repeat_step = false)
    (; k1, k2, k3, k4, k5, k6, k7, k8, utilde, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, b1, b2, b3, b4, b5, b6, btilde1, btilde2, btilde3, btilde4, btilde5, btilde6, btilde7, c2, c3, c4, c5, c6, c7) = cache.tab
    (; u, uprev, t, dt, f, p) = integrator

    @.. broadcast = false thread = thread tmp = uprev + dt * (a21 * k1)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 +
            a65 * k5
    )
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k6, tmp, p, t + c6 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 +
            a75 * k5 + a76 * k6
    )
    stage_limiter!(tmp, integrator, p, t + c7 * dt)
    f(k7, tmp, p, t + c7 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (
        b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 +
            b6 * k6
    )
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k8, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 7)
    integrator.fsallast = k8

    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = dt *
            (
            btilde1 * k1 + btilde2 * k2 +
                btilde3 * k3 + btilde4 * k4 +
                btilde5 * k5 + btilde6 * k6 +
                btilde7 * k7
        )
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    return nothing
end

function initialize!(integrator, cache::Alshina2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::Alshina2ConstantCache, repeat_step = false)
    (; u, uprev, f, p, dt, t) = integrator
    (; a21, b1, b2, b1tilde, c2) = cache

    k1 = f(uprev, p, t)
    tmp = uprev + dt * (a21 * k1)
    k2 = f(tmp, p, t + c2 * dt)
    u = uprev + dt * (b1 * k1 + b2 * k2)

    if integrator.opts.adaptive
        utilde = dt * (b1tilde * k1)
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    integrator.fsallast = k2

    integrator.k[1] = k1
    integrator.k[2] = k2
    return integrator.u = u
end

get_fsalfirstlast(cache::Alshina2Cache, u) = (cache.k1, cache.k2)
function initialize!(integrator, cache::Alshina2Cache)
    (; uprev, f, p, t) = integrator

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.fsalfirst = cache.k1
    integrator.fsallast = cache.k2

    f(integrator.fsalfirst, uprev, p, t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::Alshina2Cache, repeat_step = false)
    (; k1, k2, utilde, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, b1, b2, b1tilde, c2) = cache.tab
    (; u, uprev, t, dt, f, p) = integrator

    f(k1, uprev, p, t)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a21 * integrator.fsalfirst)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)

    @.. broadcast = false thread = thread u = uprev +
        dt * (b1 * k1 + b2 * k2)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = dt * (b1tilde * k1)
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    integrator.fsallast = k2

    return nothing
end

function initialize!(integrator, cache::Alshina3ConstantCache)
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = zero(integrator.fsalfirst)
    return integrator.k[3] = integrator.fsallast
end

function perform_step!(integrator, cache::Alshina3ConstantCache, repeat_step = false)
    (; u, uprev, f, p, dt, t) = integrator
    (; a21, a32, b1, b2, b3, b2tilde, c2, c3) = cache

    k1 = f(uprev, p, t)
    tmp = uprev + dt * (a21 * k1)
    k2 = f(tmp, p, t + c2 * dt)
    tmp = uprev + dt * (a32 * k2)
    k3 = f(tmp, p, t + c3 * dt)
    u = uprev + dt * (b1 * k1 + b2 * k2 + b3 * k3)

    if integrator.opts.adaptive
        utilde = dt * (b2tilde * k2)
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    integrator.fsallast = k3

    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    return integrator.u = u
end

get_fsalfirstlast(cache::Alshina3Cache, u) = (cache.k1, cache.k3)
function initialize!(integrator, cache::Alshina3Cache)
    (; uprev, f, p, t) = integrator

    integrator.kshortsize = 3
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.fsalfirst = cache.k1
    integrator.fsallast = cache.k3

    f(integrator.fsalfirst, uprev, p, t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::Alshina3Cache, repeat_step = false)
    (; k1, k2, k3, utilde, tmp, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; a21, a32, b1, b2, b3, b2tilde, c2, c3) = cache.tab
    (; u, uprev, t, dt, f, p) = integrator

    f(k1, uprev, p, t)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a21 * k1)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt * (b1 * k1 + b2 * k2 + b3 * k3)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = dt * (b2tilde * k2)
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    integrator.fsallast = k3

    return nothing
end

function initialize!(integrator, cache::Alshina6ConstantCache)
    integrator.kshortsize = 7
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = zero(integrator.fsalfirst)
    integrator.k[3] = zero(integrator.fsalfirst)
    integrator.k[4] = zero(integrator.fsalfirst)
    integrator.k[5] = zero(integrator.fsalfirst)
    integrator.k[6] = zero(integrator.fsalfirst)
    return integrator.k[7] = integrator.fsallast
end

function perform_step!(integrator, cache::Alshina6ConstantCache, repeat_step = false)
    (; u, uprev, f, p, dt, t) = integrator
    (;
        a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76,
        b1, b5, b6, b7, c2, c3, c4, c5, c6, c7,
    ) = cache

    k1 = f(uprev, p, t)
    tmp = uprev + dt * (a21 * k1)
    k2 = f(tmp, p, t + c2 * dt)
    tmp = uprev + dt * (a31 * k1 + a32 * k2)
    k3 = f(tmp, p, t + c3 * dt)
    tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    k4 = f(tmp, p, t + dt * c4)
    tmp = uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    k5 = f(tmp, p, t + dt * c5)
    tmp = uprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = f(tmp, p, t + dt * c6)
    tmp = uprev + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
    k7 = f(tmp, p, t + dt * c7)

    integrator.fsallast = k7

    u = uprev + dt * (b1 * k1 + b5 * k5 + b6 * k6 + b7 * k7)

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 7)

    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.k[7] = k7
    return integrator.u = u
end

get_fsalfirstlast(cache::Alshina6Cache, u) = (cache.k1, cache.k7)
function initialize!(integrator, cache::Alshina6Cache)
    (; uprev, f, p, t) = integrator

    integrator.kshortsize = 7
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.k1
    integrator.k[2] = cache.k2
    integrator.k[3] = cache.k3
    integrator.k[4] = cache.k4
    integrator.k[5] = cache.k5
    integrator.k[6] = cache.k6
    integrator.k[7] = cache.k7
    integrator.fsalfirst = cache.k1
    return integrator.fsallast = cache.k7
end

function perform_step!(integrator, cache::Alshina6Cache, repeat_step = false)
    (; k1, k2, k3, k4, k5, k6, k7, tmp, stage_limiter!, step_limiter!, thread) = cache
    (;
        a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76,
        b1, b5, b6, b7, c2, c3, c4, c5, c6, c7,
    ) = cache.tab
    (; u, uprev, t, dt, f, p) = integrator

    f(k1, uprev, p, t)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a21 * k1)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast = false thread = thread tmp = uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 +
            a65 * k5
    )
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k6, tmp, p, t + c6 * dt)
    @.. broadcast = false thread = thread tmp = uprev +
        dt * (
        a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 +
            a75 * k5 + a76 * k6
    )
    stage_limiter!(tmp, integrator, p, t + c7 * dt)
    f(k7, tmp, p, t + c7 * dt)
    @.. broadcast = false thread = thread u = uprev +
        dt *
        (b1 * k1 + b5 * k5 + b6 * k6 + b7 * k7)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 7)
    integrator.fsallast = k7
    return nothing
end
