function initialize!(integrator, ::QPRK98ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, ::QPRK98ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract QPRK98Tableau T T2

    k1 = integrator.fsalfirst
    k2 = f(uprev + b21 * k1 * dt, p, t + d2 * dt)
    k3 = f(uprev + dt * (b31 * k1 + b32 * k2), p, t + d3 * dt)
    k4 = f(uprev + dt * (b41 * k1 + b43 * k3), p, t + d4 * dt)
    k5 = f(uprev + dt * (b51 * k1 + b53 * k3 + b54 * k4), p, t + d5 * dt)
    k6 = f(uprev + dt * (b61 * k1 + b64 * k4 + b65 * k5), p, t + d6 * dt)
    k7 = f(uprev + dt * (b71 * k1 + b74 * k4 + b75 * k5 + b76 * k6), p, t + d7 * dt)
    k8 = f(uprev + dt * (b81 * k1 + b86 * k6 + b87 * k7), p, t + d8 * dt)
    k9 = f(uprev + dt * (b91 * k1 + b96 * k6 + b97 * k7 + b98 * k8), p, t + d9 * dt)
    k10 = f(uprev + dt * (b10_1 * k1 + b10_6 * k6 + b10_7 * k7 + b10_8 * k8 + b10_9 * k9),
        p, t + d10 * dt)
    k11 = f(
        uprev +
        dt * (b11_1 * k1 + b11_6 * k6 + b11_7 * k7 + b11_8 * k8 + b11_9 * k9
         + b11_10 * k10),
        p,
        t + d11 * dt)
    k12 = f(
        uprev +
        dt * (b12_1 * k1 + b12_6 * k6 + b12_7 * k7 + b12_8 * k8 + b12_9 * k9
         + b12_10 * k10 + b12_11 * k11),
        p,
        t + d12 * dt)
    k13 = f(
        uprev +
        dt * (b13_1 * k1 + b13_6 * k6 + b13_7 * k7 + b13_8 * k8 + b13_9 * k9
         + b13_10 * k10 + b13_11 * k11 + b13_12 * k12),
        p,
        t + d13 * dt)
    k14 = f(
        uprev +
        dt * (b14_1 * k1 + b14_6 * k6 + b14_7 * k7 + b14_8 * k8 + b14_9 * k9
         + b14_10 * k10 + b14_11 * k11 + b14_12 * k12 + b14_13 * k13),
        p,
        t + d14 * dt)
    k15 = f(
        uprev +
        dt * (b15_1 * k1 + b15_6 * k6 + b15_7 * k7 + b15_8 * k8 + b15_9 * k9
         + b15_10 * k10 + b15_11 * k11 + b15_12 * k12 + b15_13 * k13 + b15_14 * k14),
        p, t + dt)
    k16 = f(
        uprev +
        dt * (b16_1 * k1 + b16_6 * k6 + b16_7 * k7 + b16_8 * k8 + b16_9 * k9
         + b16_10 * k10 + b16_11 * k11 + b16_12 * k12 + b16_13 * k13 + b16_14 * k14),
        p, t + dt)
    integrator.stats.nf += 15
    u = uprev +
        dt * (w1 * k1 + w8 * k8 + w9 * k9 + w10 * k10 + w11 * k11 + w12 * k12 +
         w13 * k13 + w14 * k14 + w15 * k15 + w16 * k16)

    if integrator.opts.adaptive
        utilde = dt * (ϵ1 * k1 + ϵ8 * k8 + ϵ9 * k9 + ϵ10 * k10 + ϵ11 * k11 + ϵ12 * k12 +
                  ϵ13 * k13 + ϵ14 * k14 + ϵ15 * k15)
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    integrator.fsallast = f(u, p, t + dt)
    integrator.stats.nf += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::QPRK98Cache)
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.k
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.stats.nf += 1
end

@muladd function perform_step!(integrator, cache::QPRK98Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract QPRK98Tableau T T2
    @unpack fsalfirst, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16,
    utilde, tmp, atmp, k, stage_limiter!, step_limiter!, thread = cache
    k1 = fsalfirst
    f(k1, uprev, p, t)
    @.. broadcast=false thread=thread tmp=uprev + dt * b21 * k1
    stage_limiter!(tmp, integrator, p, t + d2 * dt)
    f(k2, tmp, p, t + d2 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (b31 * k1 + b32 * k2)
    stage_limiter!(tmp, integrator, p, t + d3 * dt)
    f(k3, tmp, p, t + d3 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (b41 * k1 + b43 * k3)
    stage_limiter!(tmp, integrator, p, t + d4 * dt)
    f(k4, tmp, p, t + d4 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (b51 * k1 + b53 * k3 + b54 * k4)
    stage_limiter!(uprev, integrator, p, t + d5 * dt)
    f(k5, tmp, p, t + d5 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (b61 * k1 + b64 * k4 + b65 * k5)
    stage_limiter!(tmp, integrator, p, t + d6 * dt)
    f(k6, tmp, p, t + d6 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (b71 * k1 + b74 * k4 + b75 * k5
                                                + b76 * k6)
    stage_limiter!(tmp, integrator, p, t + d7 * dt)
    f(k7, tmp, p, t + d7 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (b81 * k1 + b86 * k6 + b87 * k7)
    stage_limiter!(tmp, integrator, p, t + d8 * dt)
    f(k8, tmp, p, t + d8 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (b91 * k1 + b96 * k6 + b97 * k7
                                                + b98 * k8)
    stage_limiter!(tmp, integrator, p, t + d9 * dt)
    f(k9, tmp, p, t + d9 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (b10_1 * k1 + b10_6 * k6
                                           + b10_7 * k7 + b10_8 * k8 + b10_9 * k9)
    stage_limiter!(tmp, integrator, p, t + d10 * dt)
    f(k10, tmp, p, t + d10 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (b11_1 * k1 + b11_6 * k6
                                           + b11_7 * k7 + b11_8 * k8 + b11_9 * k9
                                           + b11_10 * k10)
    stage_limiter!(tmp, integrator, p, t + d11 * dt)
    f(k11, tmp, p, t + d11 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (b12_1 * k1 + b12_6 * k6
                                           + b12_7 * k7 + b12_8 * k8 + b12_9 * k9
                                           + b12_10 * k10 + b12_11 * k11)
    stage_limiter!(tmp, integrator, p, t + d12 * dt)
    f(k12, tmp, p, t + d12 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (b13_1 * k1 + b13_6 * k6
                                           + b13_7 * k7 + b13_8 * k8 + b13_9 * k9
                                           + b13_10 * k10 + b13_11 * k11 + b13_12 * k12)
    stage_limiter!(tmp, integrator, p, t + d13 * dt)
    f(k13, tmp, p, t + d13 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (b14_1 * k1 + b14_6 * k6
                                           + b14_7 * k7 + b14_8 * k8 + b14_9 * k9
                                           + b14_10 * k10 + b14_11 * k11 + b14_12 * k12
                                           + b14_13 * k13)
    stage_limiter!(tmp, integrator, p, t + d14 * dt)
    f(k14, tmp, p, t + d14 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (b15_1 * k1 + b15_6 * k6
                                           + b15_7 * k7 + b15_8 * k8 + b15_9 * k9
                                           + b15_10 * k10 + b15_11 * k11 + b15_12 * k12
                                           + b15_13 * k13 + b15_14 * k14)
    stage_limiter!(tmp, integrator, p, t + dt)
    f(k15, tmp, p, t + dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (b16_1 * k1 + b16_6 * k6
                                           + b16_7 * k7 + b16_8 * k8 + b16_9 * k9
                                           + b16_10 * k10 + b16_11 * k11 + b16_12 * k12
                                           + b16_13 * k13 + b16_14 * k14)
    stage_limiter!(u, integrator, p, t + dt)
    f(k16, tmp, p, t + dt)

    integrator.stats.nf += 16

    @.. broadcast=false thread=thread u=uprev +
                                        dt * (w1 * k1 + w8 * k8 + w9 * k9
                                         + w10 * k10 + w11 * k11 + w12 * k12 + w13 * k13
                                         + w14 * k14 + w15 * k15 + w16 * k16)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        @.. broadcast=false thread=thread utilde=dt * (ϵ1 * k1 + ϵ8 * k8 +
                                                  ϵ9 * k9 +
                                                  ϵ10 * k10 + ϵ11 * k11 +
                                                  ϵ12 * k12 +
                                                  ϵ13 * k13 + ϵ14 * k14 +
                                                  ϵ15 * k15)
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    f(k, u, p, t + dt)
    integrator.stats.nf += 1
    return nothing
end