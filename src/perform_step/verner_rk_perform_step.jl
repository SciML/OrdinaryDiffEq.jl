function initialize!(integrator, cache::Vern6ConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.stats.nf += 1
    alg = unwrap_alg(integrator, false)
    alg.lazy ? (integrator.kshortsize = 9) : (integrator.kshortsize = 12)
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:8
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    integrator.k[integrator.kshortsize] = integrator.fsallast

    if !alg.lazy
        @inbounds for i in 10:12
            integrator.k[i] = zero(integrator.fsalfirst)
        end
    end
end

@muladd function perform_step!(integrator, cache::Vern6ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack c1, c2, c3, c4, c5, c6, a21, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, a91, a94, a95, a96, a97, a98, btilde1, btilde4, btilde5, btilde6, btilde7, btilde8, btilde9 = cache.tab
    k1 = integrator.fsalfirst
    a = dt * a21
    k2 = f(uprev + a * k1, p, t + c1 * dt)
    k3 = f(uprev + dt * (a31 * k1 + a32 * k2), p, t + c2 * dt)
    k4 = f(uprev + dt * (a41 * k1 + a43 * k3), p, t + c3 * dt)
    k5 = f(uprev + dt * (a51 * k1 + a53 * k3 + a54 * k4), p, t + c4 * dt)
    k6 = f(uprev + dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5), p, t + c5 * dt)
    k7 = f(uprev + dt * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6), p,
        t + c6 * dt)
    g8 = uprev + dt * (a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 + a86 * k6 + a87 * k7)
    k8 = f(g8, p, t + dt)
    u = uprev + dt * (a91 * k1 + a94 * k4 + a95 * k5 + a96 * k6 + a97 * k7 + a98 * k8)
    integrator.fsallast = f(u, p, t + dt)
    k9 = integrator.fsallast
    integrator.stats.nf += 8
    if typeof(integrator.alg) <: CompositeAlgorithm
        g9 = u
        ϱu = integrator.opts.internalnorm(k9 - k8, t)
        ϱd = integrator.opts.internalnorm(g9 - g8, t)
        integrator.eigen_est = ϱu / ϱd
    end
    if integrator.opts.adaptive
        utilde = dt *
                 (btilde1 * k1 + btilde4 * k4 + btilde5 * k5 + btilde6 * k6 + btilde7 * k7 +
                  btilde8 * k8 + btilde9 * k9)
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
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
    integrator.k[9] = k9

    alg = unwrap_alg(integrator, false)
    if !alg.lazy && (integrator.opts.adaptive == false ||
        accept_step_controller(integrator, integrator.opts.controller))
        k = integrator.k
        @unpack c10, a1001, a1004, a1005, a1006, a1007, a1008, a1009, c11, a1101, a1104, a1105, a1106, a1107, a1108, a1109, a1110, c12, a1201, a1204, a1205, a1206, a1207, a1208, a1209, a1210, a1211 = cache.tab.extra
        k[10] = f(uprev +
                  dt * (a1001 * k[1] + a1004 * k[4] + a1005 * k[5] + a1006 * k[6] +
                   a1007 * k[7] + a1008 * k[8] + a1009 * k[9]), p, t + c10 * dt)
        k[11] = f(uprev +
                  dt * (a1101 * k[1] + a1104 * k[4] + a1105 * k[5] + a1106 * k[6] +
                   a1107 * k[7] + a1108 * k[8] + a1109 * k[9] + a1110 * k[10]), p,
            t + c11 * dt)
        k[12] = f(uprev +
                  dt * (a1201 * k[1] + a1204 * k[4] + a1205 * k[5] + a1206 * k[6] +
                   a1207 * k[7] + a1208 * k[8] + a1209 * k[9] + a1210 * k[10] +
                   a1211 * k[11]), p, t + c12 * dt)
        integrator.stats.nf += 3
    end

    integrator.u = u
end

function initialize!(integrator, cache::Vern6Cache)
    alg = unwrap_alg(integrator, false)
    alg.lazy ? (integrator.kshortsize = 9) : (integrator.kshortsize = 12)
    integrator.fsalfirst = cache.k1
    integrator.fsallast = cache.k9
    @unpack k = integrator
    resize!(k, integrator.kshortsize)
    k[1] = cache.k1
    k[2] = cache.k2
    k[3] = cache.k3
    k[4] = cache.k4
    k[5] = cache.k5
    k[6] = cache.k6
    k[7] = cache.k7
    k[8] = cache.k8
    k[9] = cache.k9 # Set the pointers

    if !alg.lazy
        k[10] = similar(cache.k1)
        k[11] = similar(cache.k1)
        k[12] = similar(cache.k1)
    end

    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.stats.nf += 1
end

@muladd function perform_step!(integrator, cache::Vern6Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    uidx = eachindex(integrator.uprev)
    @unpack c1, c2, c3, c4, c5, c6, a21, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a83, a84, a85, a86, a87, a91, a94, a95, a96, a97, a98, btilde1, btilde4, btilde5, btilde6, btilde7, btilde8, btilde9 = cache.tab
    @unpack k1, k2, k3, k4, k5, k6, k7, k8, k9, utilde, tmp, rtmp, atmp, stage_limiter!, step_limiter!, thread = cache
    a = dt * a21
    @.. broadcast=false thread=thread tmp=uprev + a * k1
    stage_limiter!(tmp, integrator, p, t + c1 * dt)
    f(k2, tmp, p, t + c1 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (a31 * k1 + a32 * k2)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k3, tmp, p, t + c2 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (a41 * k1 + a43 * k3)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k4, tmp, p, t + c3 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (a51 * k1 + a53 * k3 + a54 * k4)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k5, tmp, p, t + c4 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k6, tmp, p, t + c5 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 +
                                           a76 * k6)
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k7, tmp, p, t + c6 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a81 * k1 + a83 * k3 + a84 * k4 + a85 * k5 +
                                           a86 * k6 +
                                           a87 * k7)
    stage_limiter!(tmp, integrator, p, t + dt)
    f(k8, tmp, p, t + dt)
    @.. broadcast=false thread=thread u=uprev +
                                        dt *
                                        (a91 * k1 + a94 * k4 + a95 * k5 + a96 * k6 +
                                         a97 * k7 + a98 * k8)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k9, u, p, t + dt)
    integrator.stats.nf += 8
    if integrator.alg isa CompositeAlgorithm
        g9 = u
        g8 = tmp
        @.. broadcast=false thread=thread rtmp=k9 - k8
        ϱu = integrator.opts.internalnorm(rtmp, t)
        @.. broadcast=false thread=thread utilde=g9 - g8
        ϱd = integrator.opts.internalnorm(utilde, t)
        integrator.eigen_est = ϱu / ϱd
    end
    if integrator.opts.adaptive
        @.. broadcast=false thread=thread utilde=dt * (btilde1 * k1 + btilde4 * k4 +
                                                  btilde5 * k5 +
                                                  btilde6 * k6 + btilde7 * k7 +
                                                  btilde8 * k8 +
                                                  btilde9 * k9)
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    alg = unwrap_alg(integrator, false)
    if !alg.lazy && (integrator.opts.adaptive == false ||
        accept_step_controller(integrator, integrator.opts.controller))
        k = integrator.k
        @unpack c10, a1001, a1004, a1005, a1006, a1007, a1008, a1009, c11, a1101, a1104, a1105, a1106, a1107, a1108, a1109, a1110, c12, a1201, a1204, a1205, a1206, a1207, a1208, a1209, a1210, a1211 = cache.tab.extra
        @unpack tmp = cache
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1001 * k[1] + a1004 * k[4] + a1005 * k[5] +
                                               a1006 * k[6] +
                                               a1007 * k[7] + a1008 * k[8] + a1009 * k[9])
        f(k[10], tmp, p, t + c10 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1101 * k[1] + a1104 * k[4] + a1105 * k[5] +
                                               a1106 * k[6] +
                                               a1107 * k[7] + a1108 * k[8] + a1109 * k[9] +
                                               a1110 * k[10])
        f(k[11], tmp, p, t + c11 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1201 * k[1] + a1204 * k[4] + a1205 * k[5] +
                                               a1206 * k[6] +
                                               a1207 * k[7] + a1208 * k[8] + a1209 * k[9] +
                                               a1210 * k[10] + a1211 * k[11])
        integrator.stats.nf += 3
        f(k[12], tmp, p, t + c12 * dt)
    end
    return nothing
end

function initialize!(integrator, cache::Vern7ConstantCache)
    alg = unwrap_alg(integrator, false)
    alg.lazy ? (integrator.kshortsize = 10) : (integrator.kshortsize = 16)
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    @inbounds for i in eachindex(integrator.k)
        integrator.k[i] = zero(integrator.uprev) ./ oneunit(integrator.t)
    end
end

@muladd function perform_step!(integrator, cache::Vern7ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, k, f, p = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract Vern7Tableau T T2
    k1 = f(uprev, p, t)
    a = dt * a021
    k2 = f(uprev + a * k1, p, t + c2 * dt)
    k3 = f(uprev + dt * (a031 * k1 + a032 * k2), p, t + c3 * dt)
    k4 = f(uprev + dt * (a041 * k1 + a043 * k3), p, t + c4 * dt)
    k5 = f(uprev + dt * (a051 * k1 + a053 * k3 + a054 * k4), p, t + c5 * dt)
    k6 = f(uprev + dt * (a061 * k1 + a063 * k3 + a064 * k4 + a065 * k5), p, t + c6 * dt)
    k7 = f(uprev + dt * (a071 * k1 + a073 * k3 + a074 * k4 + a075 * k5 + a076 * k6), p,
        t + c7 * dt)
    k8 = f(uprev +
           dt * (a081 * k1 + a083 * k3 + a084 * k4 + a085 * k5 + a086 * k6 + a087 * k7), p,
        t + c8 * dt)
    g9 = uprev +
         dt *
         (a091 * k1 + a093 * k3 + a094 * k4 + a095 * k5 + a096 * k6 + a097 * k7 + a098 * k8)
    g10 = uprev +
          dt * (a101 * k1 + a103 * k3 + a104 * k4 + a105 * k5 + a106 * k6 + a107 * k7)
    k9 = f(g9, p, t + dt)
    k10 = f(g10, p, t + dt)
    integrator.stats.nf += 10
    u = uprev + dt * (b1 * k1 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7 + b8 * k8 + b9 * k9)
    if typeof(integrator.alg) <: CompositeAlgorithm
        ϱu = integrator.opts.internalnorm(k10 - k9, t)
        ϱd = integrator.opts.internalnorm(g10 - g9, t)
        integrator.eigen_est = ϱu / ϱd
    end
    if integrator.opts.adaptive
        utilde = dt *
                 (btilde1 * k1 + btilde4 * k4 + btilde5 * k5 + btilde6 * k6 + btilde7 * k7 +
                  btilde8 * k8 + btilde9 * k9 + btilde10 * k10)
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
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
    integrator.k[9] = k9
    integrator.k[10] = k10
    integrator.u = u

    alg = unwrap_alg(integrator, false)
    if !alg.lazy && (integrator.opts.adaptive == false ||
        accept_step_controller(integrator, integrator.opts.controller))
        k = integrator.k
        @OnDemandTableauExtract Vern7ExtraStages T T2
        k[11] = f(uprev +
                  dt * (a1101 * k[1] + a1104 * k[4] + a1105 * k[5] + a1106 * k[6] +
                   a1107 * k[7] + a1108 * k[8] + a1109 * k[9]), p, t + c11 * dt)
        k[12] = f(uprev +
                  dt * (a1201 * k[1] + a1204 * k[4] + a1205 * k[5] + a1206 * k[6] +
                   a1207 * k[7] + a1208 * k[8] + a1209 * k[9] + a1211 * k[11]), p,
            t + c12 * dt)
        k[13] = f(uprev +
                  dt * (a1301 * k[1] + a1304 * k[4] + a1305 * k[5] + a1306 * k[6] +
                   a1307 * k[7] + a1308 * k[8] + a1309 * k[9] + a1311 * k[11] +
                   a1312 * k[12]), p, t + c13 * dt)
        k[14] = f(uprev +
                  dt * (a1401 * k[1] + a1404 * k[4] + a1405 * k[5] + a1406 * k[6] +
                   a1407 * k[7] + a1408 * k[8] + a1409 * k[9] + a1411 * k[11] +
                   a1412 * k[12] + a1413 * k[13]), p, t + c14 * dt)
        k[15] = f(uprev +
                  dt * (a1501 * k[1] + a1504 * k[4] + a1505 * k[5] + a1506 * k[6] +
                   a1507 * k[7] + a1508 * k[8] + a1509 * k[9] + a1511 * k[11] +
                   a1512 * k[12] + a1513 * k[13]), p, t + c15 * dt)
        k[16] = f(uprev +
                  dt * (a1601 * k[1] + a1604 * k[4] + a1605 * k[5] + a1606 * k[6] +
                   a1607 * k[7] + a1608 * k[8] + a1609 * k[9] + a1611 * k[11] +
                   a1612 * k[12] + a1613 * k[13]), p, t + c16 * dt)
        integrator.stats.nf += 6
    end
end

function initialize!(integrator, cache::Vern7Cache)
    @unpack k1, k2, k3, k4, k5, k6, k7, k8, k9, k10 = cache
    @unpack k = integrator
    alg = unwrap_alg(integrator, false)
    alg.lazy ? (integrator.kshortsize = 10) : (integrator.kshortsize = 16)
    resize!(k, integrator.kshortsize)
    k[1] = k1
    k[2] = k2
    k[3] = k3
    k[4] = k4
    k[5] = k5
    k[6] = k6
    k[7] = k7
    k[8] = k8
    k[9] = k9
    k[10] = k10 # Setup pointers

    if !alg.lazy
        k[11] = similar(cache.k1)
        k[12] = similar(cache.k1)
        k[13] = similar(cache.k1)
        k[14] = similar(cache.k1)
        k[15] = similar(cache.k1)
        k[16] = similar(cache.k1)
    end
end

@muladd function perform_step!(integrator, cache::Vern7Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract Vern7Tableau T T2
    @unpack k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, utilde, tmp, rtmp, atmp, stage_limiter!, step_limiter!, thread = cache
    f(k1, uprev, p, t)
    a = dt * a021
    @.. broadcast=false thread=thread tmp=uprev + a * k1
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (a031 * k1 + a032 * k2)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (a041 * k1 + a043 * k3)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (a051 * k1 + a053 * k3 + a054 * k4)
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt *
                                          (a061 * k1 + a063 * k3 + a064 * k4 + a065 * k5)
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k6, tmp, p, t + c6 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt *
                                          (a071 * k1 + a073 * k3 + a074 * k4 + a075 * k5 +
                                           a076 * k6)
    stage_limiter!(tmp, integrator, p, t + c7 * dt)
    f(k7, tmp, p, t + c7 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt *
                                          (a081 * k1 + a083 * k3 + a084 * k4 + a085 * k5 +
                                           a086 * k6 +
                                           a087 * k7)
    stage_limiter!(tmp, integrator, p, t + c8 * dt)
    f(k8, tmp, p, t + c8 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt *
                                          (a091 * k1 + a093 * k3 + a094 * k4 + a095 * k5 +
                                           a096 * k6 +
                                           a097 * k7 + a098 * k8)
    stage_limiter!(tmp, integrator, p, t + dt)
    f(k9, tmp, p, t + dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt *
                                          (a101 * k1 + a103 * k3 + a104 * k4 + a105 * k5 +
                                           a106 * k6 +
                                           a107 * k7)
    stage_limiter!(tmp, integrator, p, t + dt)
    f(k10, tmp, p, t + dt)
    @.. broadcast=false thread=thread u=uprev +
                                        dt *
                                        (b1 * k1 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7 +
                                         b8 * k8 +
                                         b9 * k9)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    integrator.stats.nf += 10
    if integrator.alg isa CompositeAlgorithm
        g10 = u
        g9 = tmp
        @.. broadcast=false thread=thread rtmp=k10 - k9
        ϱu = integrator.opts.internalnorm(rtmp, t)
        @.. broadcast=false thread=thread utilde=g10 - g9
        ϱd = integrator.opts.internalnorm(utilde, t)
        integrator.eigen_est = ϱu / ϱd
    end
    if integrator.opts.adaptive
        @.. broadcast=false thread=thread utilde=dt * (btilde1 * k1 + btilde4 * k4 +
                                                  btilde5 * k5 +
                                                  btilde6 * k6 + btilde7 * k7 +
                                                  btilde8 * k8 +
                                                  btilde9 * k9 + btilde10 * k10)
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    alg = unwrap_alg(integrator, false)
    if !alg.lazy && (integrator.opts.adaptive == false ||
        accept_step_controller(integrator, integrator.opts.controller))
        k = integrator.k
        @unpack tmp = cache
        @OnDemandTableauExtract Vern7ExtraStages T T2
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1101 * k[1] + a1104 * k[4] + a1105 * k[5] +
                                               a1106 * k[6] +
                                               a1107 * k[7] + a1108 * k[8] + a1109 * k[9])
        f(k[11], tmp, p, t + c11 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1201 * k[1] + a1204 * k[4] + a1205 * k[5] +
                                               a1206 * k[6] +
                                               a1207 * k[7] + a1208 * k[8] + a1209 * k[9] +
                                               a1211 * k[11])
        f(k[12], tmp, p, t + c12 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1301 * k[1] + a1304 * k[4] + a1305 * k[5] +
                                               a1306 * k[6] +
                                               a1307 * k[7] + a1308 * k[8] + a1309 * k[9] +
                                               a1311 * k[11] + a1312 * k[12])
        f(k[13], tmp, p, t + c13 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1401 * k[1] + a1404 * k[4] + a1405 * k[5] +
                                               a1406 * k[6] +
                                               a1407 * k[7] + a1408 * k[8] + a1409 * k[9] +
                                               a1411 * k[11] + a1412 * k[12] +
                                               a1413 * k[13])
        f(k[14], tmp, p, t + c14 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1501 * k[1] + a1504 * k[4] + a1505 * k[5] +
                                               a1506 * k[6] +
                                               a1507 * k[7] + a1508 * k[8] + a1509 * k[9] +
                                               a1511 * k[11] + a1512 * k[12] +
                                               a1513 * k[13])
        f(k[15], tmp, p, t + c15 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1601 * k[1] + a1604 * k[4] + a1605 * k[5] +
                                               a1606 * k[6] +
                                               a1607 * k[7] + a1608 * k[8] + a1609 * k[9] +
                                               a1611 * k[11] + a1612 * k[12] +
                                               a1613 * k[13])
        f(k[16], tmp, p, t + c16 * dt)
        integrator.stats.nf += 6
    end
    return nothing
end

function initialize!(integrator, cache::Vern8ConstantCache)
    alg = unwrap_alg(integrator, false)
    alg.lazy ? (integrator.kshortsize = 13) : (integrator.kshortsize = 21)
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    @inbounds for i in eachindex(integrator.k)
        integrator.k[i] = zero(integrator.uprev) ./ oneunit(integrator.t)
    end
end

@muladd function perform_step!(integrator, cache::Vern8ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, a0201, a0301, a0302, a0401, a0403, a0501, a0503, a0504, a0601, a0604, a0605, a0701, a0704, a0705, a0706, a0801, a0804, a0805, a0806, a0807, a0901, a0904, a0905, a0906, a0907, a0908, a1001, a1004, a1005, a1006, a1007, a1008, a1009, a1101, a1104, a1105, a1106, a1107, a1108, a1109, a1110, a1201, a1204, a1205, a1206, a1207, a1208, a1209, a1210, a1211, a1301, a1304, a1305, a1306, a1307, a1308, a1309, a1310, b1, b6, b7, b8, b9, b10, b11, b12, btilde1, btilde6, btilde7, btilde8, btilde9, btilde10, btilde11, btilde12, btilde13 = cache.tab
    k1 = f(uprev, p, t)
    a = dt * a0201
    k2 = f(uprev + a * k1, p, t + c2 * dt)
    k3 = f(uprev + dt * (a0301 * k1 + a0302 * k2), p, t + c3 * dt)
    k4 = f(uprev + dt * (a0401 * k1 + a0403 * k3), p, t + c4 * dt)
    k5 = f(uprev + dt * (a0501 * k1 + a0503 * k3 + a0504 * k4), p, t + c5 * dt)
    k6 = f(uprev + dt * (a0601 * k1 + a0604 * k4 + a0605 * k5), p, t + c6 * dt)
    k7 = f(uprev + dt * (a0701 * k1 + a0704 * k4 + a0705 * k5 + a0706 * k6), p, t + c7 * dt)
    k8 = f(uprev + dt * (a0801 * k1 + a0804 * k4 + a0805 * k5 + a0806 * k6 + a0807 * k7), p,
        t + c8 * dt)
    k9 = f(uprev +
           dt *
           (a0901 * k1 + a0904 * k4 + a0905 * k5 + a0906 * k6 + a0907 * k7 + a0908 * k8), p,
        t + c9 * dt)
    k10 = f(uprev +
            dt *
            (a1001 * k1 + a1004 * k4 + a1005 * k5 + a1006 * k6 + a1007 * k7 + a1008 * k8 +
             a1009 * k9), p, t + c10 * dt)
    k11 = f(uprev +
            dt *
            (a1101 * k1 + a1104 * k4 + a1105 * k5 + a1106 * k6 + a1107 * k7 + a1108 * k8 +
             a1109 * k9 + a1110 * k10), p, t + c11 * dt)
    g12 = uprev +
          dt *
          (a1201 * k1 + a1204 * k4 + a1205 * k5 + a1206 * k6 + a1207 * k7 + a1208 * k8 +
           a1209 * k9 + a1210 * k10 + a1211 * k11)
    g13 = uprev +
          dt *
          (a1301 * k1 + a1304 * k4 + a1305 * k5 + a1306 * k6 + a1307 * k7 + a1308 * k8 +
           a1309 * k9 + a1310 * k10)
    k12 = f(g12, p, t + dt)
    k13 = f(g13, p, t + dt)
    integrator.stats.nf += 13
    u = uprev +
        dt * (b1 * k1 + b6 * k6 + b7 * k7 + b8 * k8 + b9 * k9 + b10 * k10 + b11 * k11 +
         b12 * k12)
    if typeof(integrator.alg) <: CompositeAlgorithm
        ϱu = integrator.opts.internalnorm(k13 - k12, t)
        ϱd = integrator.opts.internalnorm(g13 - g12, t)
        integrator.eigen_est = ϱu / ϱd
    end
    if integrator.opts.adaptive
        utilde = dt *
                 (btilde1 * k1 + btilde6 * k6 + btilde7 * k7 + btilde8 * k8 + btilde9 * k9 +
                  btilde10 * k10 + btilde11 * k11 + btilde12 * k12 + btilde13 * k13)
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
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
    integrator.k[9] = k9
    integrator.k[10] = k10
    integrator.k[11] = k11
    integrator.k[12] = k12
    integrator.k[13] = k13
    integrator.u = u

    alg = unwrap_alg(integrator, false)
    if !alg.lazy && (integrator.opts.adaptive == false ||
        accept_step_controller(integrator, integrator.opts.controller))
        k = integrator.k
        @unpack c14, a1401, a1406, a1407, a1408, a1409, a1410, a1411, a1412, c15, a1501, a1506, a1507, a1508, a1509, a1510, a1511, a1512, a1514, c16, a1601, a1606, a1607, a1608, a1609, a1610, a1611, a1612, a1614, a1615, c17, a1701, a1706, a1707, a1708, a1709, a1710, a1711, a1712, a1714, a1715, a1716, c18, a1801, a1806, a1807, a1808, a1809, a1810, a1811, a1812, a1814, a1815, a1816, a1817, c19, a1901, a1906, a1907, a1908, a1909, a1910, a1911, a1912, a1914, a1915, a1916, a1917, c20, a2001, a2006, a2007, a2008, a2009, a2010, a2011, a2012, a2014, a2015, a2016, a2017, c21, a2101, a2106, a2107, a2108, a2109, a2110, a2111, a2112, a2114, a2115, a2116, a2117 = cache.tab.extra
        k[14] = f(uprev +
                  dt * (a1401 * k[1] + a1406 * k[6] + a1407 * k[7] + a1408 * k[8] +
                   a1409 * k[9] + a1410 * k[10] + a1411 * k[11] + a1412 * k[12]), p,
            t + c14 * dt)
        k[15] = f(uprev +
                  dt * (a1501 * k[1] + a1506 * k[6] + a1507 * k[7] + a1508 * k[8] +
                   a1509 * k[9] + a1510 * k[10] + a1511 * k[11] + a1512 * k[12] +
                   a1514 * k[14]), p, t + c15 * dt)
        k[16] = f(uprev +
                  dt * (a1601 * k[1] + a1606 * k[6] + a1607 * k[7] + a1608 * k[8] +
                   a1609 * k[9] + a1610 * k[10] + a1611 * k[11] + a1612 * k[12] +
                   a1614 * k[14] + a1615 * k[15]), p, t + c16 * dt)
        k[17] = f(uprev +
                  dt * (a1701 * k[1] + a1706 * k[6] + a1707 * k[7] + a1708 * k[8] +
                   a1709 * k[9] + a1710 * k[10] + a1711 * k[11] + a1712 * k[12] +
                   a1714 * k[14] + a1715 * k[15] + a1716 * k[16]), p, t + c17 * dt)
        k[18] = f(uprev +
                  dt * (a1801 * k[1] + a1806 * k[6] + a1807 * k[7] + a1808 * k[8] +
                   a1809 * k[9] + a1810 * k[10] + a1811 * k[11] + a1812 * k[12] +
                   a1814 * k[14] + a1815 * k[15] + a1816 * k[16] + a1817 * k[17]), p,
            t + c18 * dt)
        k[19] = f(uprev +
                  dt * (a1901 * k[1] + a1906 * k[6] + a1907 * k[7] + a1908 * k[8] +
                   a1909 * k[9] + a1910 * k[10] + a1911 * k[11] + a1912 * k[12] +
                   a1914 * k[14] + a1915 * k[15] + a1916 * k[16] + a1917 * k[17]), p,
            t + c19 * dt)
        k[20] = f(uprev +
                  dt * (a2001 * k[1] + a2006 * k[6] + a2007 * k[7] + a2008 * k[8] +
                   a2009 * k[9] + a2010 * k[10] + a2011 * k[11] + a2012 * k[12] +
                   a2014 * k[14] + a2015 * k[15] + a2016 * k[16] + a2017 * k[17]), p,
            t + c20 * dt)
        k[21] = f(uprev +
                  dt * (a2101 * k[1] + a2106 * k[6] + a2107 * k[7] + a2108 * k[8] +
                   a2109 * k[9] + a2110 * k[10] + a2111 * k[11] + a2112 * k[12] +
                   a2114 * k[14] + a2115 * k[15] + a2116 * k[16] + a2117 * k[17]), p,
            t + c21 * dt)
        integrator.stats.nf += 8
    end
end

function initialize!(integrator, cache::Vern8Cache)
    @unpack k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13 = cache
    @unpack k = integrator
    alg = unwrap_alg(integrator, false)
    alg.lazy ? (integrator.kshortsize = 13) : (integrator.kshortsize = 21)
    resize!(k, integrator.kshortsize)
    k[1] = k1
    k[2] = k2
    k[3] = k3
    k[4] = k4
    k[5] = k5
    k[6] = k6
    k[7] = k7
    k[8] = k8
    k[9] = k9
    k[10] = k10
    k[11] = k11
    k[12] = k12
    k[13] = k13 # Setup pointers

    if !alg.lazy
        for i in 14:21
            k[i] = similar(cache.k1)
        end
    end
end

@muladd function perform_step!(integrator, cache::Vern8Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    uidx = eachindex(integrator.uprev)
    @unpack c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, a0201, a0301, a0302, a0401, a0403, a0501, a0503, a0504, a0601, a0604, a0605, a0701, a0704, a0705, a0706, a0801, a0804, a0805, a0806, a0807, a0901, a0904, a0905, a0906, a0907, a0908, a1001, a1004, a1005, a1006, a1007, a1008, a1009, a1101, a1104, a1105, a1106, a1107, a1108, a1109, a1110, a1201, a1204, a1205, a1206, a1207, a1208, a1209, a1210, a1211, a1301, a1304, a1305, a1306, a1307, a1308, a1309, a1310, b1, b6, b7, b8, b9, b10, b11, b12, btilde1, btilde6, btilde7, btilde8, btilde9, btilde10, btilde11, btilde12, btilde13 = cache.tab
    @unpack k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, utilde, tmp, rtmp, atmp, stage_limiter!, step_limiter!, thread = cache
    f(k1, uprev, p, t)
    a = dt * a0201
    @.. broadcast=false thread=thread tmp=uprev + a * k1
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k2, tmp, p, t + c2 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (a0301 * k1 + a0302 * k2)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k3, tmp, p, t + c3 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (a0401 * k1 + a0403 * k3)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k4, tmp, p, t + c4 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a0501 * k1 + a0503 * k3 + a0504 * k4)
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k5, tmp, p, t + c5 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a0601 * k1 + a0604 * k4 + a0605 * k5)
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k6, tmp, p, t + c6 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a0701 * k1 + a0704 * k4 + a0705 * k5 +
                                           a0706 * k6)
    stage_limiter!(tmp, integrator, p, t + c7 * dt)
    f(k7, tmp, p, t + c7 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt *
                                          (a0801 * k1 + a0804 * k4 + a0805 * k5 +
                                           a0806 * k6 + a0807 * k7)
    stage_limiter!(tmp, integrator, p, t + c8 * dt)
    f(k8, tmp, p, t + c8 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a0901 * k1 + a0904 * k4 + a0905 * k5 +
                                           a0906 * k6 +
                                           a0907 * k7 + a0908 * k8)
    stage_limiter!(tmp, integrator, p, t + c9 * dt)
    f(k9, tmp, p, t + c9 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a1001 * k1 + a1004 * k4 + a1005 * k5 +
                                           a1006 * k6 +
                                           a1007 * k7 + a1008 * k8 + a1009 * k9)
    stage_limiter!(tmp, integrator, p, t + c10 * dt)
    f(k10, tmp, p, t + c10 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a1101 * k1 + a1104 * k4 + a1105 * k5 +
                                           a1106 * k6 +
                                           a1107 * k7 + a1108 * k8 + a1109 * k9 +
                                           a1110 * k10)
    stage_limiter!(tmp, integrator, p, t + c11 * dt)
    f(k11, tmp, p, t + c11 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a1201 * k1 + a1204 * k4 + a1205 * k5 +
                                           a1206 * k6 +
                                           a1207 * k7 + a1208 * k8 + a1209 * k9 +
                                           a1210 * k10 +
                                           a1211 * k11)
    stage_limiter!(tmp, integrator, p, t + dt)
    f(k12, tmp, p, t + dt)
    @.. broadcast=false thread=thread u=uprev +
                                        dt *
                                        (a1301 * k1 + a1304 * k4 + a1305 * k5 + a1306 * k6 +
                                         a1307 * k7 +
                                         a1308 * k8 + a1309 * k9 + a1310 * k10)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k13, u, p, t + dt)
    integrator.stats.nf += 13
    if integrator.alg isa CompositeAlgorithm
        g13 = u
        g12 = tmp
        @.. broadcast=false thread=thread rtmp=k13 - k12
        ϱu = integrator.opts.internalnorm(rtmp, t)
        @.. broadcast=false thread=thread utilde=g13 - g12
        ϱd = integrator.opts.internalnorm(utilde, t)
        integrator.eigen_est = ϱu / ϱd
    end
    @.. broadcast=false thread=thread u=uprev +
                                        dt *
                                        (b1 * k1 + b6 * k6 + b7 * k7 + b8 * k8 + b9 * k9 +
                                         b10 * k10 +
                                         b11 * k11 + b12 * k12)
    if integrator.opts.adaptive
        @.. broadcast=false thread=thread utilde=dt * (btilde1 * k1 + btilde6 * k6 +
                                                  btilde7 * k7 +
                                                  btilde8 * k8 + btilde9 * k9 +
                                                  btilde10 * k10 +
                                                  btilde11 * k11 + btilde12 * k12 +
                                                  btilde13 * k13)
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    alg = unwrap_alg(integrator, false)
    if !alg.lazy && (integrator.opts.adaptive == false ||
        accept_step_controller(integrator, integrator.opts.controller))
        k = integrator.k
        @unpack c14, a1401, a1406, a1407, a1408, a1409, a1410, a1411, a1412, c15, a1501, a1506, a1507, a1508, a1509, a1510, a1511, a1512, a1514, c16, a1601, a1606, a1607, a1608, a1609, a1610, a1611, a1612, a1614, a1615, c17, a1701, a1706, a1707, a1708, a1709, a1710, a1711, a1712, a1714, a1715, a1716, c18, a1801, a1806, a1807, a1808, a1809, a1810, a1811, a1812, a1814, a1815, a1816, a1817, c19, a1901, a1906, a1907, a1908, a1909, a1910, a1911, a1912, a1914, a1915, a1916, a1917, c20, a2001, a2006, a2007, a2008, a2009, a2010, a2011, a2012, a2014, a2015, a2016, a2017, c21, a2101, a2106, a2107, a2108, a2109, a2110, a2111, a2112, a2114, a2115, a2116, a2117 = cache.tab.extra
        @unpack tmp = cache
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1401 * k[1] + a1406 * k[6] + a1407 * k[7] +
                                               a1408 * k[8] +
                                               a1409 * k[9] + a1410 * k[10] +
                                               a1411 * k[11] +
                                               a1412 * k[12])
        f(k[14], tmp, p, t + c14 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1501 * k[1] + a1506 * k[6] + a1507 * k[7] +
                                               a1508 * k[8] +
                                               a1509 * k[9] + a1510 * k[10] +
                                               a1511 * k[11] +
                                               a1512 * k[12] + a1514 * k[14])
        f(k[15], tmp, p, t + c15 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1601 * k[1] + a1606 * k[6] + a1607 * k[7] +
                                               a1608 * k[8] +
                                               a1609 * k[9] + a1610 * k[10] +
                                               a1611 * k[11] +
                                               a1612 * k[12] + a1614 * k[14] +
                                               a1615 * k[15])
        f(k[16], tmp, p, t + c16 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1701 * k[1] + a1706 * k[6] + a1707 * k[7] +
                                               a1708 * k[8] +
                                               a1709 * k[9] + a1710 * k[10] +
                                               a1711 * k[11] +
                                               a1712 * k[12] + a1714 * k[14] +
                                               a1715 * k[15] +
                                               a1716 * k[16])
        f(k[17], tmp, p, t + c17 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1801 * k[1] + a1806 * k[6] + a1807 * k[7] +
                                               a1808 * k[8] +
                                               a1809 * k[9] + a1810 * k[10] +
                                               a1811 * k[11] +
                                               a1812 * k[12] + a1814 * k[14] +
                                               a1815 * k[15] +
                                               a1816 * k[16] + a1817 * k[17])
        f(k[18], tmp, p, t + c18 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1901 * k[1] + a1906 * k[6] + a1907 * k[7] +
                                               a1908 * k[8] +
                                               a1909 * k[9] + a1910 * k[10] +
                                               a1911 * k[11] +
                                               a1912 * k[12] + a1914 * k[14] +
                                               a1915 * k[15] +
                                               a1916 * k[16] + a1917 * k[17])
        f(k[19], tmp, p, t + c19 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a2001 * k[1] + a2006 * k[6] + a2007 * k[7] +
                                               a2008 * k[8] +
                                               a2009 * k[9] + a2010 * k[10] +
                                               a2011 * k[11] +
                                               a2012 * k[12] + a2014 * k[14] +
                                               a2015 * k[15] +
                                               a2016 * k[16] + a2017 * k[17])
        f(k[20], tmp, p, t + c20 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a2101 * k[1] + a2106 * k[6] + a2107 * k[7] +
                                               a2108 * k[8] +
                                               a2109 * k[9] + a2110 * k[10] +
                                               a2111 * k[11] +
                                               a2112 * k[12] + a2114 * k[14] +
                                               a2115 * k[15] +
                                               a2116 * k[16] + a2117 * k[17])
        integrator.stats.nf += 8
        f(k[21], tmp, p, t + c21 * dt)
    end
    return nothing
end

function initialize!(integrator, cache::Vern9ConstantCache)
    alg = unwrap_alg(integrator, false)
    alg.lazy ? (integrator.kshortsize = 10) : (integrator.kshortsize = 20)
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Avoid undefined entries if k is an array of arrays
    @inbounds for i in eachindex(integrator.k)
        integrator.k[i] = zero(integrator.uprev) ./ oneunit(integrator.t)
    end
end

@muladd function perform_step!(integrator, cache::Vern9ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract Vern9Tableau T T2
    k1 = f(uprev, p, t)
    a = dt * a0201
    k2 = f(uprev + a * k1, p, t + c1 * dt)
    k3 = f(uprev + dt * (a0301 * k1 + a0302 * k2), p, t + c2 * dt)
    k4 = f(uprev + dt * (a0401 * k1 + a0403 * k3), p, t + c3 * dt)
    k5 = f(uprev + dt * (a0501 * k1 + a0503 * k3 + a0504 * k4), p, t + c4 * dt)
    k6 = f(uprev + dt * (a0601 * k1 + a0604 * k4 + a0605 * k5), p, t + c5 * dt)
    k7 = f(uprev + dt * (a0701 * k1 + a0704 * k4 + a0705 * k5 + a0706 * k6), p, t + c6 * dt)
    k8 = f(uprev + dt * (a0801 * k1 + a0806 * k6 + a0807 * k7), p, t + c7 * dt)
    k9 = f(uprev + dt * (a0901 * k1 + a0906 * k6 + a0907 * k7 + a0908 * k8), p, t + c8 * dt)
    k10 = f(uprev + dt * (a1001 * k1 + a1006 * k6 + a1007 * k7 + a1008 * k8 + a1009 * k9),
        p, t + c9 * dt)
    k11 = f(uprev +
            dt *
            (a1101 * k1 + a1106 * k6 + a1107 * k7 + a1108 * k8 + a1109 * k9 + a1110 * k10),
        p, t + c10 * dt)
    k12 = f(uprev +
            dt *
            (a1201 * k1 + a1206 * k6 + a1207 * k7 + a1208 * k8 + a1209 * k9 + a1210 * k10 +
             a1211 * k11), p, t + c11 * dt)
    k13 = f(uprev +
            dt *
            (a1301 * k1 + a1306 * k6 + a1307 * k7 + a1308 * k8 + a1309 * k9 + a1310 * k10 +
             a1311 * k11 + a1312 * k12), p, t + c12 * dt)
    k14 = f(uprev +
            dt *
            (a1401 * k1 + a1406 * k6 + a1407 * k7 + a1408 * k8 + a1409 * k9 + a1410 * k10 +
             a1411 * k11 + a1412 * k12 + a1413 * k13), p, t + c13 * dt)
    g15 = uprev +
          dt *
          (a1501 * k1 + a1506 * k6 + a1507 * k7 + a1508 * k8 + a1509 * k9 + a1510 * k10 +
           a1511 * k11 + a1512 * k12 + a1513 * k13 + a1514 * k14)
    g16 = uprev +
          dt *
          (a1601 * k1 + a1606 * k6 + a1607 * k7 + a1608 * k8 + a1609 * k9 + a1610 * k10 +
           a1611 * k11 + a1612 * k12 + a1613 * k13)
    k15 = f(g15, p, t + dt)
    k16 = f(g16, p, t + dt)
    integrator.stats.nf += 16
    u = uprev +
        dt * (b1 * k1 + b8 * k8 + b9 * k9 + b10 * k10 + b11 * k11 + b12 * k12 + b13 * k13 +
         b14 * k14 + b15 * k15)
    if typeof(integrator.alg) <: CompositeAlgorithm
        ϱu = integrator.opts.internalnorm(k16 - k15, t)
        ϱd = integrator.opts.internalnorm(g16 - g15, t)
        integrator.eigen_est = ϱu / ϱd
    end
    if integrator.opts.adaptive
        utilde = dt * (btilde1 * k1 + btilde8 * k8 + btilde9 * k9 + btilde10 * k10 +
                  btilde11 * k11 + btilde12 * k12 + btilde13 * k13 + btilde14 * k14 +
                  btilde15 * k15 + btilde16 * k16)
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    # k2, k3,k4,k5,k6,k7 are not used in the code (not even in interpolations), we dont need their pointers.
    # So we mapped k[2] (from integrator) with k8 (from cache), k[3] with k9 and so on.
    integrator.k[1] = k1
    integrator.k[2] = k8
    integrator.k[3] = k9
    integrator.k[4] = k10
    integrator.k[5] = k11
    integrator.k[6] = k12
    integrator.k[7] = k13
    integrator.k[8] = k14
    integrator.k[9] = k15
    integrator.k[10] = k16
    integrator.u = u

    alg = unwrap_alg(integrator, false)
    if !alg.lazy && (integrator.opts.adaptive == false ||
        accept_step_controller(integrator, integrator.opts.controller))
        k = integrator.k
        @OnDemandTableauExtract Vern9ExtraStages T T2
        k[11] = f(uprev +
                  dt * (a1701 * k[1] + a1708 * k[2] + a1709 * k[3] + a1710 * k[4] +
                   a1711 * k[5] + a1712 * k[6] + a1713 * k[7] + a1714 * k[8] +
                   a1715 * k[9]),
            p, t + c17 * dt)
        k[12] = f(uprev +
                  dt * (a1801 * k[1] + a1808 * k[2] + a1809 * k[3] + a1810 * k[4] +
                   a1811 * k[5] + a1812 * k[6] + a1813 * k[7] + a1814 * k[8] +
                   a1815 * k[9] + a1817 * k[11]), p, t + c18 * dt)
        k[13] = f(uprev +
                  dt * (a1901 * k[1] + a1908 * k[2] + a1909 * k[3] + a1910 * k[4] +
                   a1911 * k[5] + a1912 * k[6] + a1913 * k[7] + a1914 * k[8] +
                   a1915 * k[9] + a1917 * k[11] + a1918 * k[12]), p, t + c19 * dt)
        k[14] = f(uprev +
                  dt * (a2001 * k[1] + a2008 * k[2] + a2009 * k[3] + a2010 * k[4] +
                   a2011 * k[5] + a2012 * k[6] + a2013 * k[7] + a2014 * k[8] +
                   a2015 * k[9] + a2017 * k[11] + a2018 * k[12] + a2019 * k[13]), p,
            t + c20 * dt)
        k[15] = f(uprev +
                  dt * (a2101 * k[1] + a2108 * k[2] + a2109 * k[3] + a2110 * k[4] +
                   a2111 * k[5] + a2112 * k[6] + a2113 * k[7] + a2114 * k[8] +
                   a2115 * k[9] + a2117 * k[11] + a2118 * k[12] + a2119 * k[13] +
                   a2120 * k[14]), p, t + c21 * dt)
        k[16] = f(uprev +
                  dt * (a2201 * k[1] + a2208 * k[2] + a2209 * k[3] + a2210 * k[4] +
                   a2211 * k[5] + a2212 * k[6] + a2213 * k[7] + a2214 * k[8] +
                   a2215 * k[9] + a2217 * k[11] + a2218 * k[12] + a2219 * k[13] +
                   a2220 * k[14] + a2221 * k[15]), p, t + c22 * dt)
        k[17] = f(uprev +
                  dt * (a2301 * k[1] + a2308 * k[2] + a2309 * k[3] + a2310 * k[4] +
                   a2311 * k[5] + a2312 * k[6] + a2313 * k[7] + a2314 * k[8] +
                   a2315 * k[9] + a2317 * k[11] + a2318 * k[12] + a2319 * k[13] +
                   a2320 * k[14] + a2321 * k[15]), p, t + c23 * dt)
        k[18] = f(uprev +
                  dt * (a2401 * k[1] + a2408 * k[2] + a2409 * k[3] + a2410 * k[4] +
                   a2411 * k[5] + a2412 * k[6] + a2413 * k[7] + a2414 * k[8] +
                   a2415 * k[9] + a2417 * k[11] + a2418 * k[12] + a2419 * k[13] +
                   a2420 * k[14] + a2421 * k[15]), p, t + c24 * dt)
        k[19] = f(uprev +
                  dt * (a2501 * k[1] + a2508 * k[2] + a2509 * k[3] + a2510 * k[4] +
                   a2511 * k[5] + a2512 * k[6] + a2513 * k[7] + a2514 * k[8] +
                   a2515 * k[9] + a2517 * k[11] + a2518 * k[12] + a2519 * k[13] +
                   a2520 * k[14] + a2521 * k[15]), p, t + c25 * dt)
        k[20] = f(uprev +
                  dt * (a2601 * k[1] + a2608 * k[2] + a2609 * k[3] + a2610 * k[4] +
                   a2611 * k[5] + a2612 * k[6] + a2613 * k[7] + a2614 * k[8] +
                   a2615 * k[9] + a2617 * k[11] + a2618 * k[12] + a2619 * k[13] +
                   a2620 * k[14] + a2621 * k[15]), p, t + c26 * dt)
        integrator.stats.nf += 10
    end
end

function initialize!(integrator, cache::Vern9Cache)
    @unpack k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16 = cache
    @unpack k = integrator
    alg = unwrap_alg(integrator, false)
    alg.lazy ? (integrator.kshortsize = 10) : (integrator.kshortsize = 20)
    resize!(k, integrator.kshortsize)
    # k2, k3,k4,k5,k6,k7 are not used in the code (not even in interpolations), we dont need their pointers.
    # So we mapped k[2] (from integrator) with k8 (from cache), k[3] with k9 and so on.
    k[1] = k1
    k[2] = k8
    k[3] = k9
    k[4] = k10
    k[5] = k11
    k[6] = k12
    k[7] = k13
    k[8] = k14
    k[9] = k15
    k[10] = k16 # Setup pointers

    if !alg.lazy
        for i in 11:20
            k[i] = similar(cache.k1)
        end
    end
end

@muladd function perform_step!(integrator, cache::Vern9Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    uidx = eachindex(integrator.uprev)
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract Vern9Tableau T T2
    @unpack k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, utilde, tmp, rtmp, atmp, stage_limiter!, step_limiter!, thread = cache
    f(k1, uprev, p, t)
    a = dt * a0201
    @.. broadcast=false thread=thread tmp=uprev + a * k1
    stage_limiter!(tmp, integrator, p, t + c1 * dt)
    f(k2, tmp, p, t + c1 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (a0301 * k1 + a0302 * k2)
    stage_limiter!(tmp, integrator, p, t + c2 * dt)
    f(k3, tmp, p, t + c2 * dt)
    @.. broadcast=false thread=thread tmp=uprev + dt * (a0401 * k1 + a0403 * k3)
    stage_limiter!(tmp, integrator, p, t + c3 * dt)
    f(k4, tmp, p, t + c3 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a0501 * k1 + a0503 * k3 + a0504 * k4)
    stage_limiter!(tmp, integrator, p, t + c4 * dt)
    f(k5, tmp, p, t + c4 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a0601 * k1 + a0604 * k4 + a0605 * k5)
    stage_limiter!(tmp, integrator, p, t + c5 * dt)
    f(k6, tmp, p, t + c5 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a0701 * k1 + a0704 * k4 + a0705 * k5 +
                                           a0706 * k6)
    stage_limiter!(tmp, integrator, p, t + c6 * dt)
    f(k7, tmp, p, t + c6 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a0801 * k1 + a0806 * k6 + a0807 * k7)
    stage_limiter!(tmp, integrator, p, t + c7 * dt)
    f(k8, tmp, p, t + c7 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a0901 * k1 + a0906 * k6 + a0907 * k7 +
                                           a0908 * k8)
    stage_limiter!(tmp, integrator, p, t + c8 * dt)
    f(k9, tmp, p, t + c8 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt *
                                          (a1001 * k1 + a1006 * k6 + a1007 * k7 +
                                           a1008 * k8 + a1009 * k9)
    stage_limiter!(tmp, integrator, p, t + c9 * dt)
    f(k10, tmp, p, t + c9 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a1101 * k1 + a1106 * k6 + a1107 * k7 +
                                           a1108 * k8 +
                                           a1109 * k9 + a1110 * k10)
    stage_limiter!(tmp, integrator, p, t + c10 * dt)
    f(k11, tmp, p, t + c10 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a1201 * k1 + a1206 * k6 + a1207 * k7 +
                                           a1208 * k8 +
                                           a1209 * k9 + a1210 * k10 + a1211 * k11)
    stage_limiter!(tmp, integrator, p, t + c11 * dt)
    f(k12, tmp, p, t + c11 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a1301 * k1 + a1306 * k6 + a1307 * k7 +
                                           a1308 * k8 +
                                           a1309 * k9 + a1310 * k10 + a1311 * k11 +
                                           a1312 * k12)
    stage_limiter!(tmp, integrator, p, t + c12 * dt)
    f(k13, tmp, p, t + c12 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a1401 * k1 + a1406 * k6 + a1407 * k7 +
                                           a1408 * k8 +
                                           a1409 * k9 + a1410 * k10 + a1411 * k11 +
                                           a1412 * k12 +
                                           a1413 * k13)
    stage_limiter!(tmp, integrator, p, t + c13 * dt)
    f(k14, tmp, p, t + c13 * dt)
    @.. broadcast=false thread=thread tmp=uprev +
                                          dt * (a1501 * k1 + a1506 * k6 + a1507 * k7 +
                                           a1508 * k8 +
                                           a1509 * k9 + a1510 * k10 + a1511 * k11 +
                                           a1512 * k12 +
                                           a1513 * k13 + a1514 * k14)
    stage_limiter!(tmp, integrator, p, t + dt)
    f(k15, tmp, p, t + dt)
    @.. broadcast=false thread=thread u=uprev +
                                        dt *
                                        (a1601 * k1 + a1606 * k6 + a1607 * k7 + a1608 * k8 +
                                         a1609 * k9 +
                                         a1610 * k10 + a1611 * k11 + a1612 * k12 +
                                         a1613 * k13)
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    f(k16, u, p, t + dt)
    integrator.stats.nf += 16
    if integrator.alg isa CompositeAlgorithm
        g16 = u
        g15 = tmp
        @.. broadcast=false thread=thread rtmp=k16 - k15
        ϱu = integrator.opts.internalnorm(rtmp, t)
        @.. broadcast=false thread=thread utilde=g16 - g15
        ϱd = integrator.opts.internalnorm(utilde, t)
        integrator.eigen_est = ϱu / ϱd
    end
    @.. broadcast=false thread=thread u=uprev +
                                        dt *
                                        (b1 * k1 + b8 * k8 + b9 * k9 + b10 * k10 +
                                         b11 * k11 + b12 * k12 +
                                         b13 * k13 + b14 * k14 + b15 * k15)
    if integrator.opts.adaptive
        @.. broadcast=false thread=thread utilde=dt * (btilde1 * k1 + btilde8 * k8 +
                                                  btilde9 * k9 +
                                                  btilde10 * k10 + btilde11 * k11 +
                                                  btilde12 * k12 +
                                                  btilde13 * k13 + btilde14 * k14 +
                                                  btilde15 * k15 +
                                                  btilde16 * k16)
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    alg = unwrap_alg(integrator, false)
    if !alg.lazy && (integrator.opts.adaptive == false ||
        accept_step_controller(integrator, integrator.opts.controller))
        k = integrator.k
        @unpack tmp = cache
        @OnDemandTableauExtract Vern9ExtraStages T T2
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1701 * k[1] + a1708 * k[2] + a1709 * k[3] +
                                               a1710 * k[4] +
                                               a1711 * k[5] + a1712 * k[6] + a1713 * k[7] +
                                               a1714 * k[8] +
                                               a1715 * k[9])
        f(k[11], tmp, p, t + c17 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1801 * k[1] + a1808 * k[2] + a1809 * k[3] +
                                               a1810 * k[4] +
                                               a1811 * k[5] + a1812 * k[6] + a1813 * k[7] +
                                               a1814 * k[8] +
                                               a1815 * k[9] + a1817 * k[11])
        f(k[12], tmp, p, t + c18 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a1901 * k[1] + a1908 * k[2] + a1909 * k[3] +
                                               a1910 * k[4] +
                                               a1911 * k[5] + a1912 * k[6] + a1913 * k[7] +
                                               a1914 * k[8] +
                                               a1915 * k[9] + a1917 * k[11] + a1918 * k[12])
        f(k[13], tmp, p, t + c19 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a2001 * k[1] + a2008 * k[2] + a2009 * k[3] +
                                               a2010 * k[4] +
                                               a2011 * k[5] + a2012 * k[6] + a2013 * k[7] +
                                               a2014 * k[8] +
                                               a2015 * k[9] + a2017 * k[11] +
                                               a2018 * k[12] +
                                               a2019 * k[13])
        f(k[14], tmp, p, t + c20 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a2101 * k[1] + a2108 * k[2] + a2109 * k[3] +
                                               a2110 * k[4] +
                                               a2111 * k[5] + a2112 * k[6] + a2113 * k[7] +
                                               a2114 * k[8] +
                                               a2115 * k[9] + a2117 * k[11] +
                                               a2118 * k[12] +
                                               a2119 * k[13] + a2120 * k[14])
        f(k[15], tmp, p, t + c21 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a2201 * k[1] + a2208 * k[2] + a2209 * k[3] +
                                               a2210 * k[4] +
                                               a2211 * k[5] + a2212 * k[6] + a2213 * k[7] +
                                               a2214 * k[8] +
                                               a2215 * k[9] + a2217 * k[11] +
                                               a2218 * k[12] +
                                               a2219 * k[13] + a2220 * k[14] +
                                               a2221 * k[15])
        f(k[16], tmp, p, t + c22 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a2301 * k[1] + a2308 * k[2] + a2309 * k[3] +
                                               a2310 * k[4] +
                                               a2311 * k[5] + a2312 * k[6] + a2313 * k[7] +
                                               a2314 * k[8] +
                                               a2315 * k[9] + a2317 * k[11] +
                                               a2318 * k[12] +
                                               a2319 * k[13] + a2320 * k[14] +
                                               a2321 * k[15])
        f(k[17], tmp, p, t + c23 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a2401 * k[1] + a2408 * k[2] + a2409 * k[3] +
                                               a2410 * k[4] +
                                               a2411 * k[5] + a2412 * k[6] + a2413 * k[7] +
                                               a2414 * k[8] +
                                               a2415 * k[9] + a2417 * k[11] +
                                               a2418 * k[12] +
                                               a2419 * k[13] + a2420 * k[14] +
                                               a2421 * k[15])
        f(k[18], tmp, p, t + c24 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a2501 * k[1] + a2508 * k[2] + a2509 * k[3] +
                                               a2510 * k[4] +
                                               a2511 * k[5] + a2512 * k[6] + a2513 * k[7] +
                                               a2514 * k[8] +
                                               a2515 * k[9] + a2517 * k[11] +
                                               a2518 * k[12] +
                                               a2519 * k[13] + a2520 * k[14] +
                                               a2521 * k[15])
        f(k[19], tmp, p, t + c25 * dt)
        @.. broadcast=false thread=thread tmp=uprev +
                                              dt *
                                              (a2601 * k[1] + a2608 * k[2] + a2609 * k[3] +
                                               a2610 * k[4] +
                                               a2611 * k[5] + a2612 * k[6] + a2613 * k[7] +
                                               a2614 * k[8] +
                                               a2615 * k[9] + a2617 * k[11] +
                                               a2618 * k[12] +
                                               a2619 * k[13] + a2620 * k[14] +
                                               a2621 * k[15])
        integrator.stats.nf += 10
        f(k[20], tmp, p, t + c26 * dt)
    end
    return nothing
end
