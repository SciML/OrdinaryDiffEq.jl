## y'' = f(t, y, y')
## y(t₀) = y₀; y'(t₀) = y'₀
## kᵢ' = f(t₀+cᵢh, y₀+cᵢhy'₀+h²∑āᵢⱼk'ⱼ, y'₀+h∑aᵢⱼk'ⱼ)
## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
## y'₁ = y'₀ + h∑bᵢk'ᵢ

const NystromCCDefaultInitialization = Union{Nystrom4ConstantCache, FineRKN4ConstantCache,
    FineRKN5ConstantCache,
    Nystrom4VelocityIndependentConstantCache,
    Nystrom5VelocityIndependentConstantCache,
    IRKN3ConstantCache, IRKN4ConstantCache,
    DPRKN4ConstantCache, DPRKN5ConstantCache,
    DPRKN6FMConstantCache, DPRKN8ConstantCache,
    DPRKN12ConstantCache, ERKN4ConstantCache,
    ERKN5ConstantCache, ERKN7ConstantCache}

function initialize!(integrator, cache::NystromCCDefaultInitialization)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    duprev, uprev = integrator.uprev.x
    kdu = integrator.f.f1(duprev, uprev, integrator.p, integrator.t)
    ku = integrator.f.f2(duprev, uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
    integrator.fsalfirst = ArrayPartition((kdu, ku))
end

const NystromDefaultInitialization = Union{Nystrom4Cache, FineRKN4Cache, FineRKN5Cache,
    Nystrom4VelocityIndependentCache,
    Nystrom5VelocityIndependentCache,
    IRKN3Cache, IRKN4Cache,
    DPRKN4Cache, DPRKN5Cache,
    DPRKN6FMCache, DPRKN8Cache,
    DPRKN12Cache, ERKN4Cache,
    ERKN5Cache, ERKN7Cache}

function initialize!(integrator, cache::NystromDefaultInitialization)
    @unpack fsalfirst, k = cache
    duprev, uprev = integrator.uprev.x

    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.k[1].x[1], duprev, uprev, integrator.p, integrator.t)
    integrator.f.f2(integrator.k[1].x[2], duprev, uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
end

@muladd function perform_step!(integrator, cache::Nystrom4ConstantCache,
        repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    k₁ = integrator.fsalfirst.x[1]
    halfdt = dt / 2
    dtsq = dt^2
    eighth_dtsq = dtsq / 8
    half_dtsq = dtsq / 2
    ttmp = t + halfdt

    ## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
    ku = uprev + halfdt * duprev + eighth_dtsq * k₁
    ## y'₁ = y'₀ + h∑bᵢk'ᵢ
    kdu = duprev + halfdt * k₁

    k₂ = f.f1(kdu, ku, p, ttmp)
    ku = uprev + halfdt * duprev + eighth_dtsq * k₁
    kdu = duprev + halfdt * k₂

    k₃ = f.f1(kdu, ku, p, ttmp)
    ku = uprev + dt * duprev + half_dtsq * k₃
    kdu = duprev + dt * k₃

    k₄ = f.f1(kdu, ku, p, t + dt)
    u = uprev + (dtsq / 6) * (k₁ + k₂ + k₃) + dt * duprev
    du = duprev + (dt / 6) * (k₁ + k₄ + 2 * (k₂ + k₃))

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Nystrom4Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, fsalfirst, k₂, k₃, k₄, k = cache
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    k₁ = integrator.fsalfirst.x[1]
    halfdt = dt / 2
    dtsq = dt^2
    eighth_dtsq = dtsq / 8
    half_dtsq = dtsq / 2
    ttmp = t + halfdt

    ## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
    @.. broadcast=false ku=uprev + halfdt * duprev + eighth_dtsq * k₁
    ## y'₁ = y'₀ + h∑bᵢk'ᵢ
    @.. broadcast=false kdu=duprev + halfdt * k₁

    f.f1(k₂, kdu, ku, p, ttmp)
    @.. broadcast=false ku=uprev + halfdt * duprev + eighth_dtsq * k₁
    @.. broadcast=false kdu=duprev + halfdt * k₂

    f.f1(k₃, kdu, ku, p, ttmp)
    @.. broadcast=false ku=uprev + dt * duprev + half_dtsq * k₃
    @.. broadcast=false kdu=duprev + dt * k₃

    f.f1(k₄, kdu, ku, p, t + dt)
    @.. broadcast=false u=uprev + (dtsq / 6) * (k₁ + k₂ + k₃) + dt * duprev
    @.. broadcast=false du=duprev + (dt / 6) * (k₁ + k₄ + 2 * (k₂ + k₃))

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
end

@muladd function perform_step!(integrator, cache::FineRKN4ConstantCache,
        repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c2, c3, c4, c5, a21, a31, a32, a41, a43, a51,
    a52, a53, a54, abar21, abar31, abar32, abar41, abar42, abar43, abar51,
    abar52, abar53, abar54, b1, b3, b4, b5, bbar1, bbar3, bbar4, bbar5, btilde1, btilde3, btilde4, btilde5, bptilde1,
    bptilde3, bptilde4, bptilde5 = cache
    k1 = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c2 * duprev + dt * (a21 * k1))
    kdu = duprev + dt * (abar21 * k1)

    k2 = f.f1(kdu, ku, p, t + dt * c2)
    ku = uprev + dt * (c3 * duprev + dt * (a31 * k1 + a32 * k2))
    kdu = duprev + dt * (abar31 * k1 + abar32 * k2)

    k3 = f.f1(kdu, ku, p, t + dt * c3)
    ku = uprev + dt * (c4 * duprev + dt * (a41 * k1 + a43 * k3)) # a42 = 0
    kdu = duprev + dt * (abar41 * k1 + abar42 * k2 + abar43 * k3)

    k4 = f.f1(kdu, ku, p, t + dt * c4)
    ku = uprev + dt * (c5 * duprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))
    kdu = duprev + dt * (abar51 * k1 + abar52 * k2 + abar53 * k3 + abar54 * k4)

    k5 = f.f1(kdu, ku, p, t + dt * c5)

    u = uprev + dt * (duprev + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5)) # b2 = 0
    du = duprev + dt * (bbar1 * k1 + bbar3 * k3 + bbar4 * k4 + bbar5 * k5) # bbar2 = 0

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 5
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq * (btilde1 * k1 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5) # btilde2 = 0
        duhat = dt * (bptilde1 * k1 + bptilde3 * k3 + bptilde4 * k4 + bptilde5 * k5) # bptilde2 = 0
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::FineRKN4Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, atmp, fsalfirst, k2, k3, k4, k5, k, utilde = cache
    @unpack c2, c3, c4, c5, a21, a31, a32, a41, a43, a51,
    a52, a53, a54, abar21, abar31, abar32, abar41, abar42, abar43, abar51,
    abar52, abar53, abar54, b1, b3, b4, b5, bbar1, bbar3, bbar4, bbar5, btilde1, btilde3, btilde4, btilde5, bptilde1,
    bptilde3, bptilde4, bptilde5 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a21 * k1))
    @.. broadcast=false kdu=duprev + dt * (abar21 * k1)

    f.f1(k2, kdu, ku, p, t + dt * c2)
    @.. broadcast=false ku=uprev + dt * (c3 * duprev + dt * (a31 * k1 + a32 * k2))
    @.. broadcast=false kdu=duprev + dt * (abar31 * k1 + abar32 * k2)

    f.f1(k3, kdu, ku, p, t + dt * c3)
    @.. broadcast=false ku=uprev +
                           dt * (c4 * duprev + dt * (a41 * k1 + a43 * k3)) # a42 = 0
    @.. broadcast=false kdu=duprev + dt * (abar41 * k1 + abar42 * k2 + abar43 * k3)

    f.f1(k4, kdu, ku, p, t + dt * c4)
    @.. broadcast=false ku=uprev +
                           dt *
                           (c5 * duprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))
    @.. broadcast=false kdu=duprev +
                            dt * (abar51 * k1 + abar52 * k2 + abar53 * k3 + abar54 * k4)

    f.f1(k5, kdu, ku, p, t + dt * c5)
    @.. broadcast=false u=uprev +
                          dt * (duprev + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5)) # b2 = 0
    @.. broadcast=false du=duprev +
                           dt *
                           (bbar1 * k1 + bbar3 * k3 + bbar4 * k4 + bbar5 * k5) # bbar2 = 0

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 5
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @.. broadcast=false uhat=dtsq *
                                 (btilde1 * k1 + btilde3 * k3 + btilde4 * k4 +
                                  btilde5 * k5) # btilde2 = 0
        @.. broadcast=false duhat=dt *
                                  (bptilde1 * k1 + bptilde3 * k3 + bptilde4 * k4 +
                                   bptilde5 * k5) # bptilde2 = 0

        calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::FineRKN5ConstantCache,
        repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c2, c3, c4, c5, c6, c7, a21, a31, a32, a41, a43, a51, a52, a53, a54, a61, a62, a63, a64, a71, a73, a74, a75, abar21, abar31, abar32, abar41, abar42, abar43, abar51, abar52, abar53, abar54, abar61, abar62, abar63, abar64, abar65, abar71, abar73, abar74, abar75, abar76, b1, b3, b4, b5, bbar1, bbar3, bbar4, bbar5, bbar6, btilde1, btilde3, btilde4, btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6, bptilde7 = cache
    k1 = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c2 * duprev + dt * (a21 * k1))
    kdu = duprev + dt * (abar21 * k1)

    k2 = f.f1(kdu, ku, p, t + dt * c2)
    ku = uprev + dt * (c3 * duprev + dt * (a31 * k1 + a32 * k2))
    kdu = duprev + dt * (abar31 * k1 + abar32 * k2)

    k3 = f.f1(kdu, ku, p, t + dt * c3)
    ku = uprev + dt * (c4 * duprev + dt * (a41 * k1 + a43 * k3)) # a42 = 0
    kdu = duprev + dt * (abar41 * k1 + abar42 * k2 + abar43 * k3)

    k4 = f.f1(kdu, ku, p, t + dt * c4)
    ku = uprev + dt * (c5 * duprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))
    kdu = duprev + dt * (abar51 * k1 + abar52 * k2 + abar53 * k3 + abar54 * k4)

    k5 = f.f1(kdu, ku, p, t + dt * c5)
    ku = uprev +
         dt * (c6 * duprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4)) # a65 = 0
    kdu = duprev +
          dt * (abar61 * k1 + abar62 * k2 + abar63 * k3 + abar64 * k4 + abar65 * k5)

    k6 = f.f1(kdu, ku, p, t + dt * c6)
    ku = uprev +
         dt * (c7 * duprev +
               dt * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5)) # a72 = a76 = 0
    kdu = duprev +
          dt * (abar71 * k1 + abar73 * k3 + abar74 * k4 + abar75 * k5 +
                abar76 * k6) # abar72 = 0

    k7 = f.f1(kdu, ku, p, t + dt * c7)
    u = uprev + dt * (duprev + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5)) # no b6, b7
    du = duprev + dt * (bbar1 * k1 + bbar3 * k3 + bbar4 * k4 + bbar5 * k5 + bbar6 * k6) # no b2, b7

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 7
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq * (btilde1 * k1 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5)
        duhat = dt * (bptilde1 * k1 + bptilde3 * k3 + bptilde4 * k4 + bptilde5 * k5 +
                 bptilde6 * k6 + bptilde7 * k7)
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::FineRKN5Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, atmp, fsalfirst, k2, k3, k4, k5, k6, k7, k, utilde = cache
    @unpack c1, c2, c3, c4, c5, c6, c7, a21, a31, a32, a41, a43, a51, a52, a53, a54, a61, a62, a63, a64, a71, a73, a74, a75, abar21, abar31, abar32, abar41, abar42, abar43, abar51, abar52, abar53, abar54, abar61, abar62, abar63, abar64, abar65, abar71, abar73, abar74, abar75, abar76, b1, b3, b4, b5, bbar1, bbar3, bbar4, bbar5, bbar6, btilde1, btilde3, btilde4, btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6, bptilde7 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a21 * k1))
    @.. broadcast=false kdu=duprev + dt * (abar21 * k1)

    f.f1(k2, kdu, ku, p, t + dt * c2)
    @.. broadcast=false ku=uprev + dt * (c3 * duprev + dt * (a31 * k1 + a32 * k2))
    @.. broadcast=false kdu=duprev + dt * (abar31 * k1 + abar32 * k2)

    f.f1(k3, kdu, ku, p, t + dt * c3)
    @.. broadcast=false ku=uprev +
                           dt * (c4 * duprev + dt * (a41 * k1 + a43 * k3)) # a42 = 0
    @.. broadcast=false kdu=duprev + dt * (abar41 * k1 + abar42 * k2 + abar43 * k3)

    f.f1(k4, kdu, ku, p, t + dt * c4)
    @.. broadcast=false ku=uprev +
                           dt *
                           (c5 * duprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))
    @.. broadcast=false kdu=duprev +
                            dt * (abar51 * k1 + abar52 * k2 + abar53 * k3 + abar54 * k4)

    f.f1(k5, kdu, ku, p, t + dt * c5)
    @.. broadcast=false ku=uprev +
                           dt * (c6 * duprev +
                            dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4)) # a65 = 0
    @.. broadcast=false kdu=duprev +
                            dt * (abar61 * k1 + abar62 * k2 + abar63 * k3 + abar64 * k4 +
                             abar65 * k5)

    f.f1(k6, kdu, ku, p, t + dt * c6)
    @.. broadcast=false ku=uprev +
                           dt * (c7 * duprev +
                            dt * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5)) # a72 = a76 = 0
    @.. broadcast=false kdu=duprev +
                            dt * (abar71 * k1 + abar73 * k3 + abar74 * k4 +
                             abar75 * k5 + abar76 * k6) # abar72 = 0

    f.f1(k7, kdu, ku, p, t + dt * c7)
    @.. broadcast=false u=uprev +
                          dt * (duprev + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5))
    @.. broadcast=false du=duprev +
                           dt *
                           (bbar1 * k1 + bbar3 * k3 + bbar4 * k4 + bbar5 * k5 + bbar6 * k6)

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 7
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @.. broadcast=false uhat=dtsq *
                                 (btilde1 * k1 + btilde3 * k3 + btilde4 * k4 +
                                  btilde5 * k5)
        @.. broadcast=false duhat=dt *
                                  (bptilde1 * k1 + bptilde3 * k3 + bptilde4 * k4 +
                                   bptilde5 * k5 + bptilde6 * k6 + bptilde7 * k7)

        calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::Nystrom4VelocityIndependentConstantCache,
        repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    k₁ = integrator.fsalfirst.x[1]
    halfdt = dt / 2
    dtsq = dt^2
    eighth_dtsq = dtsq / 8
    half_dtsq = dtsq / 2
    ttmp = t + halfdt

    ## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
    ku = uprev + halfdt * duprev + eighth_dtsq * k₁

    k₂ = f.f1(duprev, ku, p, ttmp)
    ku = uprev + dt * duprev + half_dtsq * k₂

    k₃ = f.f1(duprev, ku, p, t + dt)
    u = uprev + (dtsq / 6) * (k₁ + 2 * k₂) + dt * duprev
    du = duprev + (dt / 6) * (k₁ + k₃ + 4 * k₂)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 3
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Nystrom4VelocityIndependentCache,
        repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, fsalfirst, k₂, k₃, k = cache
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    k₁ = integrator.fsalfirst.x[1]
    halfdt = dt / 2
    dtsq = dt^2
    eighth_dtsq = dtsq / 8
    half_dtsq = dtsq / 2
    ttmp = t + halfdt

    ## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
    @.. broadcast=false ku=uprev + halfdt * duprev + eighth_dtsq * k₁

    f.f1(k₂, duprev, ku, p, ttmp)
    @.. broadcast=false ku=uprev + dt * duprev + half_dtsq * k₂

    f.f1(k₃, duprev, ku, p, t + dt)
    @.. broadcast=false u=uprev + (dtsq / 6) * (k₁ + 2 * k₂) + dt * duprev
    @.. broadcast=false du=duprev + (dt / 6) * (k₁ + k₃ + 4 * k₂)

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 3
    integrator.stats.nf2 += 1
end

@muladd function perform_step!(integrator, cache::IRKN3ConstantCache, repeat_step = false)
    @unpack t, dt, k, tprev, f, p = integrator
    duprev, uprev = integrator.uprev.x
    duprev2, uprev2 = integrator.uprev2.x
    @unpack bconst1, bconst2, c1, a21, b1, b2, bbar1, bbar2 = cache
    k₁ = integrator.fsalfirst
    # if there's a discontinuity or the solver is in the first step
    if integrator.iter < 2 && !integrator.u_modified
        perform_step!(integrator, Nystrom4VelocityIndependentConstantCache())
        k = integrator.fsallast
        k1cache = ArrayPartition((k.x[1], f.f1(duprev, uprev, p, t + c1 * dt)))
        kdu = uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[1])
        k₂.x[1] = f.f1(duprev, kdu, p, t + c1 * dt)
        integrator.stats.nf += 2
    else
        kdu = uprev2 + dt * (c1 * duprev2 + dt * a21 * k1cache.x[1])
        ku = uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[2])

        k₂x1 = f.f1(duprev, ku, p, t + c1 * dt)
        du = duprev +
             dt * (b1 * k1cache.x[1] + bbar1 * k1cache.x[1] + b2 * (k₂x1 - k₂.x[1]))
        u = uprev + bconst1 * dt * duprev +
            dt * (bconst2 * duprev2 + dt * bbar2 * (k₂x1 - k₂.x[1]))

        integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt),
            f.f2(du, u, p, t + dt)))
        integrator.stats.nf += 3
        integrator.stats.nf2 += 1
        copyto!(k₂.x[1], k₂.x[2])
        k1cache = ArrayPartition((k1cache.x[1], k.x[2]))
    end # end if
end

@muladd function perform_step!(integrator, cache::IRKN3Cache, repeat_step = false)
    @unpack t, dt, k, tprev, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    duprev2, uprev2 = integrator.uprev2.x
    uidx = eachindex(integrator.uprev.x[1])
    @unpack tmp, fsalfirst, k₂, k = cache
    @unpack bconst1, bconst2, c1, a21, b1, b2, bbar1, bbar2 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    k1cache = cache.tmp2
    k₁ = fsalfirst
    # if there's a discontinuity or the solver is in the first step
    if integrator.iter < 2 && !integrator.u_modified
        perform_step!(integrator, integrator.cache.onestep_cache)
        copyto!(k1cache.x[1], k.x[1])
        f.f1(k1cache.x[2], duprev, uprev, p, t + c1 * dt)
        @.. broadcast=false kdu=uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[2])
        f.f1(k₂.x[1], duprev, kdu, p, t + c1 * dt)
        integrator.stats.nf += 2
    else
        @.. broadcast=false kdu=uprev2 + dt * (c1 * duprev2 + dt * a21 * k1cache.x[1])
        @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[2])

        f.f1(k₂.x[2], duprev, ku, p, t + c1 * dt)
        @tight_loop_macros for i in uidx
            @inbounds u[i] = uprev[i] + bconst1 * dt * duprev[i] +
                             dt *
                             (bconst2 * duprev2[i] + dt * bbar2 * (k₂.x[2][i] - k₂.x[1][i]))
            @inbounds du[i] = duprev[i] +
                              dt * (b1 * k1cache.x[1][i] + bbar1 * k1cache.x[2][i] +
                               b2 * (k₂.x[2][i] - k₂.x[1][i]))
        end
        f.f1(k.x[1], du, u, p, t + dt)
        f.f2(k.x[2], du, u, p, t + dt)
        integrator.stats.nf += 3
        integrator.stats.nf2 += 1
        copyto!(k₂.x[1], k₂.x[2])
        copyto!(k1cache.x[2], k1cache.x[1])
        copyto!(k1cache.x[1], k.x[1])
    end # end if
end

@muladd function perform_step!(integrator, cache::IRKN4Cache, repeat_step = false)
    @unpack t, dt, k, tprev, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    duprev2, uprev2 = integrator.uprev2.x
    uidx = eachindex(integrator.uprev.x[1])
    @unpack tmp, tmp2, fsalfirst, k₂, k₃, k = cache
    @unpack bconst1, bconst2, c1, c2, a21, a32, b1, b2, b3, bbar1, bbar2, bbar3 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    k1cache = integrator.cache.tmp2
    k₁ = fsalfirst
    # if there's a discontinuity or the solver is in the first step
    if integrator.iter < 2 && !integrator.u_modified
        perform_step!(integrator, integrator.cache.onestep_cache)
        copyto!(k1cache.x[1], k.x[1])
        f.f1(k1cache.x[2], duprev, uprev, p, t + c1 * dt)
        @.. broadcast=false kdu=uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[1])
        f.f1(k₂.x[1], duprev, kdu, p, t + c1 * dt)
        @.. broadcast=false kdu=uprev + dt * (c2 * duprev + dt * a32 * k1cache.x[2])
        f.f1(k₃.x[1], duprev, kdu, p, t + c1 * dt)
        integrator.stats.nf += 3
    else
        @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[1])
        @.. broadcast=false kdu=uprev2 + dt * (c1 * duprev2 + dt * a21 * k1cache.x[2])

        f.f1(k₂.x[2], duprev, ku, p, t + c1 * dt)
        @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * a32 * k₂.x[2])
        @.. broadcast=false kdu=uprev2 + dt * (c2 * duprev2 + dt * a32 * k₂.x[1])

        f.f1(k₃.x[2], duprev, ku, p, t + c2 * dt)
        @tight_loop_macros for i in uidx
            @inbounds u[i] = uprev[i] + dt * bconst1 * duprev[i] +
                             dt * (bconst2 * duprev2[i] +
                              dt * (bbar2 * (k₂.x[2][i] - k₂.x[1][i]) +
                               bbar3 * (k₃.x[2][i] - k₃.x[1][i])))
            @inbounds du[i] = duprev[i] +
                              dt * (b1 * k1cache.x[1][i] + bbar1 * k1cache.x[2][i] +
                               b2 * (k₂.x[2][i] - k₂.x[1][i]) +
                               b3 * (k₃.x[2][i] - k₃.x[1][i]))
        end
        f.f1(k.x[1], du, u, p, t + dt)
        f.f2(k.x[2], du, u, p, t + dt)
        integrator.stats.nf += 4
        integrator.stats.nf2 += 1
        copyto!(k₂.x[1], k₂.x[2])
        copyto!(k₃.x[1], k₃.x[2])
        copyto!(k1cache.x[2], k1cache.x[1])
        copyto!(k1cache.x[1], k.x[1])
    end # end if
end

@muladd function perform_step!(integrator, cache::Nystrom5VelocityIndependentConstantCache,
        repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c1, c2, a21, a31, a32, a41, a42, a43, bbar1, bbar2, bbar3, b1, b2, b3, b4 = cache
    k₁ = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c1 * duprev + dt * a21 * k₁)

    k₂ = f.f1(duprev, ku, p, t + c1 * dt)
    ku = uprev + dt * (c2 * duprev + dt * (a31 * k₁ + a32 * k₂))

    k₃ = f.f1(duprev, ku, p, t + c2 * dt)
    ku = uprev + dt * (duprev + dt * (a41 * k₁ + a42 * k₂ + a43 * k₃))

    k₄ = f.f1(duprev, ku, p, t + dt)
    u = uprev + dt * (duprev + dt * (bbar1 * k₁ + bbar2 * k₂ + bbar3 * k₃))
    du = duprev + dt * (b1 * k₁ + b2 * k₂ + b3 * k₃ + b4 * k₄)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Nystrom5VelocityIndependentCache,
        repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    uidx = eachindex(integrator.uprev.x[1])
    @unpack tmp, fsalfirst, k₂, k₃, k₄, k = cache
    @unpack c1, c2, a21, a31, a32, a41, a42, a43, bbar1, bbar2, bbar3, b1, b2, b3, b4 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    k₁ = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k₁)

    f.f1(k₂, du, ku, p, t + c1 * dt)
    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a31 * k₁ + a32 * k₂))

    f.f1(k₃, du, ku, p, t + c2 * dt)
    #@tight_loop_macros for i in uidx
    #  @inbounds ku[i] = uprev[i] + dt*(duprev[i] + dt*(a41*k₁[i] + a42*k₂[i] + a43*k₃[i]))
    #end
    @.. broadcast=false ku=uprev + dt * (duprev + dt * (a41 * k₁ + a42 * k₂ + a43 * k₃))

    f.f1(k₄, duprev, ku, p, t + dt)
    #@tight_loop_macros for i in uidx
    #  @inbounds u[i]  = uprev[i] + dt*(duprev[i] + dt*(bbar1*k₁[i] + bbar2*k₂[i] + bbar3*k₃[i]))
    #  @inbounds du[i] = duprev[i] + dt*(b1*k₁[i] + b2*k₂[i] + b3*k₃[i] + b4*k₄[i])
    #end
    @.. broadcast=false u=uprev +
                          dt * (duprev + dt * (bbar1 * k₁ + bbar2 * k₂ + bbar3 * k₃))
    @.. broadcast=false du=duprev + dt * (b1 * k₁ + b2 * k₂ + b3 * k₃ + b4 * k₄)
    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
    return nothing
end

@muladd function perform_step!(integrator, cache::DPRKN4ConstantCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c1, c2, c3, a21, a31, a32, a41, a42, a43, b1, b2, b3, bp1, bp2, bp3, bp4, btilde1, btilde2, btilde3, btilde4, bptilde1, bptilde2, bptilde3, bptilde4 = cache
    k1 = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c1 * duprev + dt * a21 * k1)

    k2 = f.f1(duprev, ku, p, t + dt * c1)
    ku = uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    k3 = f.f1(duprev, ku, p, t + dt * c2)
    ku = uprev + dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    k4 = f.f1(duprev, ku, p, t + dt * c3)

    u = uprev + dt * (duprev + dt * (b1 * k1 + b2 * k2 + b3 * k3))
    du = duprev + dt * (bp1 * k1 + bp2 * k2 + bp3 * k3 + bp4 * k4)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq * (btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4)
        duhat = dt * (bptilde1 * k1 + bptilde2 * k2 + bptilde3 * k3 + bptilde4 * k4)
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::DPRKN4Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, atmp, fsalfirst, k2, k3, k4, k, utilde = cache
    @unpack c1, c2, c3, a21, a31, a32, a41, a42, a43, b1, b2, b3, bp1, bp2, bp3, bp4, btilde1, btilde2, btilde3, btilde4, bptilde1, bptilde2, bptilde3, bptilde4 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k1)

    f.f1(k2, duprev, ku, p, t + dt * c1)
    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    f.f1(k3, duprev, ku, p, t + dt * c2)
    @.. broadcast=false ku=uprev +
                           dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    f.f1(k4, duprev, ku, p, t + dt * c3)
    @tight_loop_macros for i in uidx
        @inbounds u[i] = uprev[i] +
                         dt * (duprev[i] +
                               dt * (b1 * k1[i] + b2 * k2[i] + b3 * k3[i]))
        @inbounds du[i] = duprev[i] +
                          dt * (bp1 * k1[i] + bp2 * k2[i] + bp3 * k3[i] + bp4 * k4[i])
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @tight_loop_macros for i in uidx
            @inbounds uhat[i] = dtsq *
                                (btilde1 * k1[i] + btilde2 * k2[i] + btilde3 * k3[i] +
                                 btilde4 * k4[i])
            @inbounds duhat[i] = dt *
                                 (bptilde1 * k1[i] + bptilde2 * k2[i] + bptilde3 * k3[i] +
                                  bptilde4 * k4[i])
        end
        calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::DPRKN5ConstantCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c1, c2, c3, c4, c5, a21, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1, bp3, bp4, bp5, bp6, btilde1, btilde3, btilde4, btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6 = cache
    k1 = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c1 * duprev + dt * a21 * k1)

    k2 = f.f1(duprev, ku, p, t + dt * c1)
    ku = uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    k3 = f.f1(duprev, ku, p, t + dt * c2)
    ku = uprev + dt * (c3 * duprev + dt * (a41 * k1 + a43 * k3))

    k4 = f.f1(duprev, ku, p, t + dt * c3)
    ku = uprev + dt * (c4 * duprev + dt * (a51 * k1 + a53 * k3 + a54 * k4))

    k5 = f.f1(duprev, ku, p, t + dt * c4)
    ku = uprev +
         dt * (c5 * duprev + dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5))

    k6 = f.f1(duprev, ku, p, t + dt * c5)
    u = uprev +
        dt * (duprev + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5))
    du = duprev +
         dt * (bp1 * k1 + bp3 * k3 + bp4 * k4 + bp5 * k5 + bp6 * k6)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 6
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq * (btilde1 * k1 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5)
        duhat = dt * (bptilde1 * k1 + bptilde3 * k3 + bptilde4 * k4 + bptilde5 * k5 +
                 bptilde6 * k6)
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::DPRKN5Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, atmp, fsalfirst, k2, k3, k4, k5, k6, k, utilde = cache
    @unpack c1, c2, c3, c4, c5, a21, a31, a32, a41, a43, a51, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1, bp3, bp4, bp5, bp6, btilde1, btilde3, btilde4, btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k1)

    f.f1(k2, duprev, ku, p, t + dt * c1)
    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    f.f1(k3, duprev, ku, p, t + dt * c2)
    @.. broadcast=false ku=uprev +
                           dt * (c3 * duprev + dt * (a41 * k1 + a43 * k3))

    f.f1(k4, duprev, ku, p, t + dt * c3)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c4 * duprev[i] +
                           dt * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]))
    end

    f.f1(k5, duprev, ku, p, t + dt * c4)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c5 * duprev[i] +
                           dt * (a61 * k1[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]))
    end

    f.f1(k6, duprev, ku, p, t + dt * c5)
    @tight_loop_macros for i in uidx
        @inbounds u[i] = uprev[i] +
                         dt * (duprev[i] +
                          dt * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i]))
        @inbounds du[i] = duprev[i] +
                          dt * (bp1 * k1[i] + bp3 * k3[i] + bp4 * k4[i] + bp5 * k5[i] +
                           bp6 * k6[i])
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 6
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @tight_loop_macros for i in uidx
            @inbounds uhat[i] = dtsq *
                                (btilde1 * k1[i] + btilde3 * k3[i] + btilde4 * k4[i] +
                                 btilde5 * k5[i])
            @inbounds duhat[i] = dt *
                                 (bptilde1 * k1[i] + bptilde3 * k3[i] + bptilde4 * k4[i] +
                                  bptilde5 * k5[i] + bptilde6 * k6[i])
        end
        calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

function initialize!(integrator, cache::DPRKN6ConstantCache)
    duprev, uprev = integrator.uprev.x
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    kdu = integrator.f.f1(duprev, uprev, integrator.p, integrator.t)
    ku = integrator.f.f2(duprev, uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
    integrator.fsalfirst = ArrayPartition((kdu, ku))
    integrator.fsallast = zero(integrator.fsalfirst)

    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::DPRKN6ConstantCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1, bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4, btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6 = cache
    k1 = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c1 * duprev + dt * a21 * k1)

    k2 = f.f1(duprev, ku, p, t + dt * c1)
    ku = uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    k3 = f.f1(duprev, ku, p, t + dt * c2)
    ku = uprev + dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    k4 = f.f1(duprev, ku, p, t + dt * c3)
    ku = uprev + dt * (c4 * duprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))

    k5 = f.f1(duprev, ku, p, t + dt * c4)
    ku = uprev + dt * (c5 * duprev + dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)) # no a62

    k6 = f.f1(duprev, ku, p, t + dt * c5)
    u = uprev + dt * (duprev + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5)) # b1 -- b5, no b2
    du = duprev + dt * (bp1 * k1 + bp3 * k3 + bp4 * k4 + bp5 * k5 + bp6 * k6) # bp1 -- bp6, no bp2

    #=
    @tight_loop_macros for i in uidx
      @inbounds u[i]  = uprev[i] + dt*(duprev[i] + dt*(bhat1*k1.x[2][i] + bhat2*k2.x[2][i] + bhat3*k3.x[2][i]))
      @inbounds du[i] = duprev[i]+ dt*(bphat1*k1.x[2][i] + bphat3*k3.x[2][i] + bphat4*k4.x[2][i] + bphat5*k5.x[2][i] + bphat6*k6.x[2][i])
    end
    =#

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.k[1] = ArrayPartition(integrator.fsalfirst.x[1], k2)
    integrator.k[2] = ArrayPartition(k3, k4)
    integrator.k[3] = ArrayPartition(k5, k6)
    integrator.stats.nf += 6
    integrator.stats.nf2 += 1

    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq *
               (btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5)
        duhat = dt * (bptilde1 * k1 + bptilde3 * k3 + bptilde4 * k4 + bptilde5 * k5 +
                 bptilde6 * k6)
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

function initialize!(integrator, cache::DPRKN6Cache)
    @unpack fsalfirst, k = cache
    duprev, uprev = integrator.uprev.x

    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    integrator.kshortsize = 3
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = ArrayPartition(cache.fsalfirst.x[1], cache.k2)
    integrator.k[2] = ArrayPartition(cache.k3, cache.k4)
    integrator.k[3] = ArrayPartition(cache.k5, cache.k6)
    integrator.f.f1(integrator.fsallast.x[1], duprev, uprev, integrator.p, integrator.t)
    integrator.f.f2(integrator.fsallast.x[2], duprev, uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
end

@muladd function perform_step!(integrator, cache::DPRKN6Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, atmp, fsalfirst, k2, k3, k4, k5, k6, k, utilde = cache
    @unpack c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1, bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4, btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k1)

    f.f1(k2, du, ku, p, t + dt * c1)
    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    f.f1(k3, du, ku, p, t + dt * c2)
    @.. broadcast=false ku=uprev +
                           dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    f.f1(k4, du, ku, p, t + dt * c3)
    @.. broadcast=false ku=uprev +
                           dt *
                           (c4 * duprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))

    f.f1(k5, du, ku, p, t + dt * c4)
    @.. broadcast=false ku=uprev +
                           dt *
                           (c5 * duprev + dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)) # no a62

    f.f1(k6, du, ku, p, t + dt * c5)

    @.. broadcast=false u=uprev +
                          dt * (duprev + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5)) # b1 -- b5, no b2
    @.. broadcast=false du=duprev +
                           dt * (bp1 * k1 + bp3 * k3 + bp4 * k4 + bp5 * k5 + bp6 * k6) # bp1 -- bp6, no bp2

    #=
    @tight_loop_macros for i in uidx
      @inbounds u[i]  = uprev[i] + dt*(duprev[i] + dt*(bhat1*k1.x[2][i] + bhat2*k2.x[2][i] + bhat3*k3.x[2][i]))
      @inbounds du[i] = duprev[i]+ dt*(bphat1*k1.x[2][i] + bphat3*k3.x[2][i] + bphat4*k4.x[2][i] + bphat5*k5.x[2][i] + bphat6*k6.x[2][i])
    end
    =#

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 6
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @.. broadcast=false uhat=dtsq * (btilde1 * k1 + btilde2 * k2 + btilde3 * k3 +
                                  btilde4 * k4 + btilde5 * k5)
        @.. broadcast=false duhat=dt * (bptilde1 * k1 + bptilde3 * k3 + bptilde4 * k4 +
                                   bptilde5 * k5 + bptilde6 * k6)
        calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::DPRKN6FMConstantCache,
        repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, b1, b2, b3, b4, b5, bp1, bp2, bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4, btilde5, bptilde1, bptilde2, bptilde3, bptilde4, bptilde5 = cache
    k1 = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c1 * duprev + dt * a21 * k1)

    k2 = f.f1(duprev, ku, p, t + dt * c1)
    ku = uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    k3 = f.f1(duprev, ku, p, t + dt * c2)
    ku = uprev + dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    k4 = f.f1(duprev, ku, p, t + dt * c3)
    ku = uprev + dt * (c4 * duprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))

    k5 = f.f1(duprev, ku, p, t + dt * c4)
    ku = uprev +
         dt * (c5 * duprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5))

    k6 = f.f1(duprev, ku, p, t + dt * c5)
    u = uprev +
        dt * (duprev + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5))
    du = duprev +
         dt * (bp1 * k1 + bp2 * k2 + bp3 * k3 + bp4 * k4 + bp5 * k5 + bp6 * k6)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 6
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq *
               (btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5)
        duhat = dt * (bptilde1 * k1 + bptilde2 * k2 + bptilde3 * k3 + bptilde4 * k4 +
                 bptilde5 * k5)
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::DPRKN6FMCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, atmp, fsalfirst, k2, k3, k4, k5, k6, k, utilde = cache
    @unpack c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, b1, b2, b3, b4, b5, bp1, bp2, bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4, btilde5, bptilde1, bptilde2, bptilde3, bptilde4, bptilde5 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k1)

    f.f1(k2, duprev, ku, p, t + dt * c1)
    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    f.f1(k3, duprev, ku, p, t + dt * c2)
    @.. broadcast=false ku=uprev +
                           dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    f.f1(k4, duprev, ku, p, t + dt * c3)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c4 * duprev[i] +
                           dt * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]))
    end

    f.f1(k5, duprev, ku, p, t + dt * c4)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c5 * duprev[i] +
                           dt * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] +
                            a65 * k5[i]))
    end

    f.f1(k6, duprev, ku, p, t + dt * c5)
    @tight_loop_macros for i in uidx
        @inbounds u[i] = uprev[i] +
                         dt * (duprev[i] +
                          dt *
                          (b1 * k1[i] + b2 * k2[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i]))
        @inbounds du[i] = duprev[i] +
                          dt * (bp1 * k1[i] + bp2 * k2[i] + bp3 * k3[i] + bp4 * k4[i] +
                           bp5 * k5[i] + bp6 * k6[i])
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 6
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @tight_loop_macros for i in uidx
            @inbounds uhat[i] = dtsq *
                                (btilde1 * k1[i] + btilde2 * k2[i] + btilde3 * k3[i] +
                                 btilde4 * k4[i] + btilde5 * k5[i])
            @inbounds duhat[i] = dt *
                                 (bptilde1 * k1[i] + bptilde2 * k2[i] + bptilde3 * k3[i] +
                                  bptilde4 * k4[i] + bptilde5 * k5[i])
        end
        calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::DPRKN8ConstantCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c1, c2, c3, c4, c5, c6, c7, c8, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a82, a83, a84, a85, a86, a87, a91, a93, a94, a95, a96, a97, b1, b3, b4, b5, b6, b7, bp1, bp3, bp4, bp5, bp6, bp7, bp8, btilde1, btilde3, btilde4, btilde5, btilde6, btilde7, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6, bptilde7, bptilde8, bptilde9 = cache
    k1 = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c1 * duprev + dt * a21 * k1)

    k2 = f.f1(duprev, ku, p, t + dt * c1)
    ku = uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    k3 = f.f1(duprev, ku, p, t + dt * c2)
    ku = uprev + dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    k4 = f.f1(duprev, ku, p, t + dt * c3)
    ku = uprev + dt * (c4 * duprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))

    k5 = f.f1(duprev, ku, p, t + dt * c4)
    ku = uprev +
         dt * (c5 * duprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5))

    k6 = f.f1(duprev, ku, p, t + dt * c5)
    ku = uprev +
         dt * (c6 * duprev +
          dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6))

    k7 = f.f1(duprev, ku, p, t + dt * c6)
    ku = uprev +
         dt * (c7 * duprev +
          dt * (a81 * k1 + a82 * k2 + a83 * k3 + a84 * k4 + a85 * k5 + a86 * k6 + a87 * k7))

    k8 = f.f1(duprev, ku, p, t + dt * c7)
    ku = uprev +
         dt * (c8 * duprev +
          dt * (a91 * k1 + a93 * k3 + a94 * k4 + a95 * k5 + a96 * k6 + a97 * k7)) # no a92 & a98

    k9 = f.f1(duprev, ku, p, t + dt * c8)
    u = uprev +
        dt * (duprev + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7)) # b1 -- b7, no b2
    du = duprev +
         dt * (bp1 * k1 + bp3 * k3 + bp4 * k4 + bp5 * k5 + bp6 * k6 + bp7 * k7 + bp8 * k8) # bp1 -- bp8, no bp2

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 9
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq *
               (btilde1 * k1 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5 + btilde6 * k6 +
                btilde7 * k7)
        duhat = dt * (bptilde1 * k1 + bptilde3 * k3 + bptilde4 * k4 + bptilde5 * k5 +
                 bptilde6 * k6 + bptilde7 * k7 + bptilde8 * k8 + bptilde9 * k9)
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::DPRKN8Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, atmp, fsalfirst, k2, k3, k4, k5, k6, k7, k8, k9, k, utilde = cache
    @unpack c1, c2, c3, c4, c5, c6, c7, c8, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a82, a83, a84, a85, a86, a87, a91, a93, a94, a95, a96, a97, b1, b3, b4, b5, b6, b7, bp1, bp3, bp4, bp5, bp6, bp7, bp8, btilde1, btilde3, btilde4, btilde5, btilde6, btilde7, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6, bptilde7, bptilde8, bptilde9 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k1)

    f.f1(k2, duprev, ku, p, t + dt * c1)
    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    f.f1(k3, duprev, ku, p, t + dt * c2)
    @.. broadcast=false ku=uprev +
                           dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    f.f1(k4, duprev, ku, p, t + dt * c3)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c4 * duprev[i] +
                           dt * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]))
    end

    f.f1(k5, duprev, ku, p, t + dt * c4)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c5 * duprev[i] +
                           dt * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] +
                            a65 * k5[i]))
    end

    f.f1(k6, duprev, ku, p, t + dt * c5)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c6 * duprev[i] +
                           dt * (a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] +
                            a75 * k5[i] + a76 * k6[i]))
    end

    f.f1(k7, duprev, ku, p, t + dt * c6)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c7 * duprev[i] +
                           dt * (a81 * k1[i] + a82 * k2[i] + a83 * k3[i] + a84 * k4[i] +
                            a85 * k5[i] + a86 * k6[i] + a87 * k7[i]))
    end

    f.f1(k8, duprev, ku, p, t + dt * c7)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c8 * duprev[i] +
                           dt * (a91 * k1[i] + a93 * k3[i] + a94 * k4[i] + a95 * k5[i] +
                            a96 * k6[i] + a97 * k7[i])) # no a92 & a98
    end

    f.f1(k9, duprev, ku, p, t + dt * c8)
    @tight_loop_macros for i in uidx
        @inbounds u[i] = uprev[i] +
                         dt * (duprev[i] +
                          dt *
                          (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i] +
                           b7 * k7[i])) # b1 -- b7, no b2
        @inbounds du[i] = duprev[i] +
                          dt * (bp1 * k1[i] + bp3 * k3[i] + bp4 * k4[i] + bp5 * k5[i] +
                           bp6 * k6[i] + bp7 * k7[i] + bp8 * k8[i]) # bp1 -- bp8, no bp2
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 9
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @tight_loop_macros for i in uidx
            @inbounds uhat[i] = dtsq *
                                (btilde1 * k1[i] + btilde3 * k3[i] + btilde4 * k4[i] +
                                 btilde5 * k5[i] + btilde6 * k6[i] + btilde7 * k7[i])
            @inbounds duhat[i] = dt *
                                 (bptilde1 * k1[i] + bptilde3 * k3[i] + bptilde4 * k4[i] +
                                  bptilde5 * k5[i] + bptilde6 * k6[i] + bptilde7 * k7[i] +
                                  bptilde8 * k8[i] + bptilde9 * k9[i])
        end
        calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::DPRKN12ConstantCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, a21, a31, a32, a41, a42, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a84, a85, a86, a87, a91, a93, a94, a95, a96, a97, a98, a101, a103, a104, a105, a106, a107, a108, a109, a111, a113, a114, a115, a116, a117, a118, a119, a1110, a121, a123, a124, a125, a126, a127, a128, a129, a1210, a1211, a131, a133, a134, a135, a136, a137, a138, a139, a1310, a1311, a1312, a141, a143, a144, a145, a146, a147, a148, a149, a1410, a1411, a1412, a1413, a151, a153, a154, a155, a156, a157, a158, a159, a1510, a1511, a1512, a1513, a1514, a161, a163, a164, a165, a166, a167, a168, a169, a1610, a1611, a1612, a1613, a1614, a1615, a171, a173, a174, a175, a176, a177, a178, a179, a1710, a1711, a1712, a1713, a1714, a1715, b1, b7, b8, b9, b10, b11, b12, b13, b14, b15, bp1, bp7, bp8, bp9, bp10, bp11, bp12, bp13, bp14, bp15, bp16, bp17, btilde1, btilde7, btilde8, btilde9, btilde10, btilde11, btilde12, btilde13, btilde14, btilde15, bptilde1, bptilde7, bptilde8, bptilde9, bptilde10, bptilde11, bptilde12, bptilde13, bptilde14, bptilde15, bptilde16, bptilde17 = cache
    k1 = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c1 * duprev + dt * a21 * k1)

    k2 = f.f1(duprev, ku, p, t + dt * c1)
    ku = uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    k3 = f.f1(duprev, ku, p, t + dt * c2)
    ku = uprev + dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    k4 = f.f1(duprev, ku, p, t + dt * c3)
    ku = uprev + dt * (c4 * duprev + dt * (a51 * k1 + a53 * k3 + a54 * k4)) # no a52

    k5 = f.f1(duprev, ku, p, t + dt * c4)
    ku = uprev + dt * (c5 * duprev + dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)) # no a62

    k6 = f.f1(duprev, ku, p, t + dt * c5)
    ku = uprev +
         dt * (c6 * duprev + dt * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)) # no a72

    k7 = f.f1(duprev, ku, p, t + dt * c6)
    ku = uprev +
         dt * (c7 * duprev + dt * (a81 * k1 + a84 * k4 + a85 * k5 + a86 * k6 + a87 * k7)) # no a82, a83

    k8 = f.f1(duprev, ku, p, t + dt * c7)
    ku = uprev +
         dt * (c8 * duprev +
          dt * (a91 * k1 + a93 * k3 + a94 * k4 + a95 * k5 + a96 * k6 + a97 * k7 + a98 * k8)) # no a92

    k9 = f.f1(duprev, ku, p, t + dt * c8)
    ku = uprev +
         dt * (c9 * duprev +
          dt * (a101 * k1 + a103 * k3 + a104 * k4 + a105 * k5 + a106 * k6 + a107 * k7 +
           a108 * k8 + a109 * k9)) # no a102

    k10 = f.f1(duprev, ku, p, t + dt * c9)
    ku = uprev +
         dt * (c10 * duprev +
          dt * (a111 * k1 + a113 * k3 + a114 * k4 + a115 * k5 + a116 * k6 + a117 * k7 +
           a118 * k8 + a119 * k9 + a1110 * k10)) # no a112

    k11 = f.f1(duprev, ku, p, t + dt * c10)
    ku = uprev +
         dt * (c11 * duprev +
          dt * (a121 * k1 + a123 * k3 + a124 * k4 + a125 * k5 + a126 * k6 + a127 * k7 +
           a128 * k8 + a129 * k9 + a1210 * k10 + a1211 * k11)) # no a122

    k12 = f.f1(duprev, ku, p, t + dt * c11)
    ku = uprev +
         dt * (c12 * duprev +
          dt * (a131 * k1 + a133 * k3 + a134 * k4 + a135 * k5 + a136 * k6 + a137 * k7 +
           a138 * k8 + a139 * k9 + a1310 * k10 + a1311 * k11 + a1312 * k12)) # no a132

    k13 = f.f1(duprev, ku, p, t + dt * c12)
    ku = uprev +
         dt * (c13 * duprev +
          dt * (a141 * k1 + a143 * k3 + a144 * k4 + a145 * k5 + a146 * k6 + a147 * k7 +
           a148 * k8 + a149 * k9 + a1410 * k10 + a1411 * k11 + a1412 * k12 + a1413 * k13)) # no a142

    k14 = f.f1(duprev, ku, p, t + dt * c13)
    ku = uprev +
         dt * (c14 * duprev +
          dt * (a151 * k1 + a153 * k3 + a154 * k4 + a155 * k5 + a156 * k6 + a157 * k7 +
           a158 * k8 + a159 * k9 + a1510 * k10 + a1511 * k11 + a1512 * k12 + a1513 * k13 +
           a1514 * k14)) # no a152

    k15 = f.f1(duprev, ku, p, t + dt * c14)
    ku = uprev +
         dt * (c15 * duprev +
          dt * (a161 * k1 + a163 * k3 + a164 * k4 + a165 * k5 + a166 * k6 + a167 * k7 +
           a168 * k8 + a169 * k9 + a1610 * k10 + a1611 * k11 + a1612 * k12 + a1613 * k13 +
           a1614 * k14 + a1615 * k15)) # no a162

    k16 = f.f1(duprev, ku, p, t + dt * c15)
    ku = uprev +
         dt * (c16 * duprev +
          dt * (a171 * k1 + a173 * k3 + a174 * k4 + a175 * k5 + a176 * k6 + a177 * k7 +
           a178 * k8 + a179 * k9 + a1710 * k10 + a1711 * k11 + a1712 * k12 + a1713 * k13 +
           a1714 * k14 + a1715 * k15)) # no a172, a1716

    k17 = f.f1(duprev, ku, p, t + dt * c16)
    u = uprev +
        dt * (duprev +
         dt * (b1 * k1 + b7 * k7 + b8 * k8 + b9 * k9 + b10 * k10 + b11 * k11 + b12 * k12 +
          b13 * k13 + b14 * k14 + b15 * k15)) # b1 & b7 -- b15
    du = duprev +
         dt *
         (bp1 * k1 + bp7 * k7 + bp8 * k8 + bp9 * k9 + bp10 * k10 + bp11 * k11 + bp12 * k12 +
          bp13 * k13 + bp14 * k14 + bp15 * k15 + bp16 * k16 + bp17 * k17) # bp1 & bp7 -- bp17

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 17
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq *
               (btilde1 * k1 + btilde7 * k7 + btilde8 * k8 + btilde9 * k9 + btilde10 * k10 +
                btilde11 * k11 + btilde12 * k12 + btilde13 * k13 + btilde14 * k14 +
                btilde15 * k15) # btilde1 & btilde7 -- btilde15
        duhat = dt * (bptilde1 * k1 + bptilde7 * k7 + bptilde8 * k8 + bptilde9 * k9 +
                 bptilde10 * k10 + bptilde11 * k11 + bptilde12 * k12 + bptilde13 * k13 +
                 bptilde14 * k14 + bptilde15 * k15 + bptilde16 * k16 + bptilde17 * k17) # bptilde1 & bptilde7 -- bptilde17
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::DPRKN12Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, atmp, fsalfirst, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k, utilde = cache
    @unpack c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, a21, a31, a32, a41, a42, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a84, a85, a86, a87, a91, a93, a94, a95, a96, a97, a98, a101, a103, a104, a105, a106, a107, a108, a109, a111, a113, a114, a115, a116, a117, a118, a119, a1110, a121, a123, a124, a125, a126, a127, a128, a129, a1210, a1211, a131, a133, a134, a135, a136, a137, a138, a139, a1310, a1311, a1312, a141, a143, a144, a145, a146, a147, a148, a149, a1410, a1411, a1412, a1413, a151, a153, a154, a155, a156, a157, a158, a159, a1510, a1511, a1512, a1513, a1514, a161, a163, a164, a165, a166, a167, a168, a169, a1610, a1611, a1612, a1613, a1614, a1615, a171, a173, a174, a175, a176, a177, a178, a179, a1710, a1711, a1712, a1713, a1714, a1715, b1, b7, b8, b9, b10, b11, b12, b13, b14, b15, bp1, bp7, bp8, bp9, bp10, bp11, bp12, bp13, bp14, bp15, bp16, bp17, btilde1, btilde7, btilde8, btilde9, btilde10, btilde11, btilde12, btilde13, btilde14, btilde15, bptilde1, bptilde7, bptilde8, bptilde9, bptilde10, bptilde11, bptilde12, bptilde13, bptilde14, bptilde15, bptilde16, bptilde17 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k1)

    f.f1(k2, duprev, ku, p, t + dt * c1)
    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    f.f1(k3, duprev, ku, p, t + dt * c2)
    @.. broadcast=false ku=uprev +
                           dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    f.f1(k4, duprev, ku, p, t + dt * c3)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt *
                          (c4 * duprev[i] + dt * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i])) # no a52
    end

    f.f1(k5, duprev, ku, p, t + dt * c4)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c5 * duprev[i] +
                           dt * (a61 * k1[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i])) # no a62
    end

    f.f1(k6, duprev, ku, p, t + dt * c5)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c6 * duprev[i] +
                           dt * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] +
                            a76 * k6[i])) # no a72
    end

    f.f1(k7, duprev, ku, p, t + dt * c6)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c7 * duprev[i] +
                           dt * (a81 * k1[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] +
                            a87 * k7[i])) # no a82, a83
    end

    f.f1(k8, duprev, ku, p, t + dt * c7)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c8 * duprev[i] +
                           dt * (a91 * k1[i] + a93 * k3[i] + a94 * k4[i] + a95 * k5[i] +
                            a96 * k6[i] + a97 * k7[i] + a98 * k8[i])) # no a92
    end

    f.f1(k9, duprev, ku, p, t + dt * c8)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c9 * duprev[i] +
                           dt *
                           (a101 * k1[i] + a103 * k3[i] + a104 * k4[i] + a105 * k5[i] +
                            a106 * k6[i] + a107 * k7[i] + a108 * k8[i] + a109 * k9[i])) # no a102
    end

    f.f1(k10, duprev, ku, p, t + dt * c9)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c10 * duprev[i] +
                           dt *
                           (a111 * k1[i] + a113 * k3[i] + a114 * k4[i] + a115 * k5[i] +
                            a116 * k6[i] + a117 * k7[i] + a118 * k8[i] + a119 * k9[i] +
                            a1110 * k10[i])) # no a112
    end

    f.f1(k11, duprev, ku, p, t + dt * c10)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c11 * duprev[i] +
                           dt *
                           (a121 * k1[i] + a123 * k3[i] + a124 * k4[i] + a125 * k5[i] +
                            a126 * k6[i] + a127 * k7[i] + a128 * k8[i] + a129 * k9[i] +
                            a1210 * k10[i] + a1211 * k11[i])) # no a122
    end

    f.f1(k12, duprev, ku, p, t + dt * c11)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c12 * duprev[i] +
                           dt *
                           (a131 * k1[i] + a133 * k3[i] + a134 * k4[i] + a135 * k5[i] +
                            a136 * k6[i] + a137 * k7[i] + a138 * k8[i] + a139 * k9[i] +
                            a1310 * k10[i] + a1311 * k11[i] + a1312 * k12[i])) # no a132
    end

    f.f1(k13, duprev, ku, p, t + dt * c12)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c13 * duprev[i] +
                           dt *
                           (a141 * k1[i] + a143 * k3[i] + a144 * k4[i] + a145 * k5[i] +
                            a146 * k6[i] + a147 * k7[i] + a148 * k8[i] + a149 * k9[i] +
                            a1410 * k10[i] + a1411 * k11[i] + a1412 * k12[i] +
                            a1413 * k13[i])) # no a142
    end

    f.f1(k14, duprev, ku, p, t + dt * c13)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c14 * duprev[i] +
                           dt *
                           (a151 * k1[i] + a153 * k3[i] + a154 * k4[i] + a155 * k5[i] +
                            a156 * k6[i] + a157 * k7[i] + a158 * k8[i] + a159 * k9[i] +
                            a1510 * k10[i] + a1511 * k11[i] + a1512 * k12[i] +
                            a1513 * k13[i] + a1514 * k14[i])) # no a152
    end

    f.f1(k15, duprev, ku, p, t + dt * c14)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c15 * duprev[i] +
                           dt *
                           (a161 * k1[i] + a163 * k3[i] + a164 * k4[i] + a165 * k5[i] +
                            a166 * k6[i] + a167 * k7[i] + a168 * k8[i] + a169 * k9[i] +
                            a1610 * k10[i] + a1611 * k11[i] + a1612 * k12[i] +
                            a1613 * k13[i] + a1614 * k14[i] + a1615 * k15[i])) # no a162
    end

    f.f1(k16, duprev, ku, p, t + dt * c15)
    @tight_loop_macros for i in uidx
        @inbounds ku[i] = uprev[i] +
                          dt * (c16 * duprev[i] +
                           dt *
                           (a171 * k1[i] + a173 * k3[i] + a174 * k4[i] + a175 * k5[i] +
                            a176 * k6[i] + a177 * k7[i] + a178 * k8[i] + a179 * k9[i] +
                            a1710 * k10[i] + a1711 * k11[i] + a1712 * k12[i] +
                            a1713 * k13[i] + a1714 * k14[i] + a1715 * k15[i])) # no a172, a1716
    end

    f.f1(k17, duprev, ku, p, t + dt * c16)
    @tight_loop_macros for i in uidx
        @inbounds u[i] = uprev[i] +
                         dt * (duprev[i] +
                          dt * (b1 * k1[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i] +
                           b10 * k10[i] + b11 * k11[i] + b12 * k12[i] + b13 * k13[i] +
                           b14 * k14[i] + b15 * k15[i])) # b1 & b7 -- b15
        @inbounds du[i] = duprev[i] +
                          dt * (bp1 * k1[i] + bp7 * k7[i] + bp8 * k8[i] + bp9 * k9[i] +
                           bp10 * k10[i] + bp11 * k11[i] + bp12 * k12[i] + bp13 * k13[i] +
                           bp14 * k14[i] + bp15 * k15[i] + bp16 * k16[i] + bp17 * k17[i]) # bp1 & bp7 -- bp17
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 17
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @tight_loop_macros for i in uidx
            @inbounds uhat[i] = dtsq *
                                (btilde1 * k1[i] + btilde7 * k7[i] + btilde8 * k8[i] +
                                 btilde9 * k9[i] + btilde10 * k10[i] + btilde11 * k11[i] +
                                 btilde12 * k12[i] + btilde13 * k13[i] + btilde14 * k14[i] +
                                 btilde15 * k15[i]) # btilde1 & btilde7 -- btilde15
            @inbounds duhat[i] = dt *
                                 (bptilde1 * k1[i] + bptilde7 * k7[i] + bptilde8 * k8[i] +
                                  bptilde9 * k9[i] + bptilde10 * k10[i] +
                                  bptilde11 * k11[i] + bptilde12 * k12[i] +
                                  bptilde13 * k13[i] + bptilde14 * k14[i] +
                                  bptilde15 * k15[i] + bptilde16 * k16[i] +
                                  bptilde17 * k17[i]) # bptilde1 & bptilde7 -- bptilde17
        end
        calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::ERKN4ConstantCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c1, c2, c3, a21, a31, a32, a41, a42, a43, b1, b2, b3, b4, bp1, bp2, bp3, bp4, btilde1, btilde2, btilde3, btilde4, bptilde1, bptilde2, bptilde3, bptilde4 = cache
    k1 = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c1 * duprev + dt * a21 * k1)

    k2 = f.f1(duprev, ku, p, t + dt * c1)
    ku = uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    k3 = f.f1(duprev, ku, p, t + dt * c2)
    ku = uprev + dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    k4 = f.f1(duprev, ku, p, t + dt * c3)
    u = uprev + dt * (duprev + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4))
    du = duprev + dt * (bp1 * k1 + bp2 * k2 + bp3 * k3 + bp4 * k4)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq * (btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4)
        duhat = dt * (bptilde1 * k1 + bptilde2 * k2 + bptilde3 * k3 + bptilde4 * k4)
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::ERKN4Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, atmp, fsalfirst, k2, k3, k4, k, utilde = cache
    @unpack c1, c2, c3, a21, a31, a32, a41, a42, a43, b1, b2, b3, b4, bp1, bp2, bp3, bp4, btilde1, btilde2, btilde3, btilde4, bptilde1, bptilde2, bptilde3, bptilde4 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k1)

    f.f1(k2, duprev, ku, p, t + dt * c1)
    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    f.f1(k3, duprev, ku, p, t + dt * c2)
    @.. broadcast=false ku=uprev +
                           dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    f.f1(k4, duprev, ku, p, t + dt * c3)
    @tight_loop_macros for i in uidx
        @inbounds u[i] = uprev[i] +
                         dt * (duprev[i] +
                          dt * (b1 * k1[i] + b2 * k2[i] + b3 * k3[i] + b4 * k4[i]))
        @inbounds du[i] = duprev[i] +
                          dt * (bp1 * k1[i] + bp2 * k2[i] + bp3 * k3[i] + bp4 * k4[i])
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @tight_loop_macros for i in uidx
            @inbounds uhat[i] = dtsq *
                                (btilde1 * k1[i] + btilde2 * k2[i] + btilde3 * k3[i] +
                                 btilde4 * k4[i])
            @inbounds duhat[i] = dt *
                                 (bptilde1 * k1[i] + bptilde2 * k2[i] + bptilde3 * k3[i] +
                                  bptilde4 * k4[i])
        end
        calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::ERKN5ConstantCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c1, c2, c3, a21, a31, a32, a41, a42, a43, b1, b2, b3, b4, bp1, bp2, bp3, bp4, btilde1, btilde2, btilde3, btilde4 = cache
    k1 = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c1 * duprev + dt * a21 * k1)

    k2 = f.f1(duprev, ku, p, t + dt * c1)
    ku = uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    k3 = f.f1(duprev, ku, p, t + dt * c2)
    ku = uprev + dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    k4 = f.f1(duprev, ku, p, t + dt * c3)
    u = uprev + dt * (duprev + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4))
    du = duprev + dt * (bp1 * k1 + bp2 * k2 + bp3 * k3 + bp4 * k4)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq * (btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4)
        atmp = calculate_residuals(uhat, integrator.uprev.x[2], integrator.u.x[2],
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::ERKN5Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, atmp, fsalfirst, k2, k3, k4, k, utilde = cache
    @unpack c1, c2, c3, a21, a31, a32, a41, a42, a43, b1, b2, b3, b4, bp1, bp2, bp3, bp4, btilde1, btilde2, btilde3, btilde4 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k1)

    f.f1(k2, duprev, ku, p, t + dt * c1)
    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    f.f1(k3, duprev, ku, p, t + dt * c2)
    @.. broadcast=false ku=uprev +
                           dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    f.f1(k4, duprev, ku, p, t + dt * c3)
    @tight_loop_macros for i in uidx
        @inbounds u[i] = uprev[i] +
                         dt * (duprev[i] +
                          dt * (b1 * k1[i] + b2 * k2[i] + b3 * k3[i] + b4 * k4[i]))
        @inbounds du[i] = duprev[i] +
                          dt * (bp1 * k1[i] + bp2 * k2[i] + bp3 * k3[i] + bp4 * k4[i])
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @tight_loop_macros for i in uidx
            @inbounds uhat[i] = dtsq *
                                (btilde1 * k1[i] + btilde2 * k2[i] + btilde3 * k3[i] +
                                 btilde4 * k4[i])
        end
        calculate_residuals!(atmp.x[2], uhat, integrator.uprev.x[2], integrator.u.x[2],
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp.x[2], t)
    end
end

@muladd function perform_step!(integrator, cache::ERKN7ConstantCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    @unpack c1, c2, c3, c4, c5, c6, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a73, a74, a75, a76, b1, b3, b4, b5, b6, bp1, bp3, bp4, bp5, bp6, bp7, btilde1, btilde3, btilde4, btilde5, btilde6, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6, bptilde7 = cache
    k1 = integrator.fsalfirst.x[1]

    ku = uprev + dt * (c1 * duprev + dt * a21 * k1)

    k2 = f.f1(duprev, ku, p, t + dt * c1)
    ku = uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    k3 = f.f1(duprev, ku, p, t + dt * c2)
    ku = uprev + dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    k4 = f.f1(duprev, ku, p, t + dt * c3)
    ku = uprev + dt * (c4 * duprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))

    k5 = f.f1(duprev, ku, p, t + dt * c4)
    ku = uprev +
         dt * (c5 * duprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5))

    k6 = f.f1(duprev, ku, p, t + dt * c5)
    ku = uprev +
         dt * (c6 * duprev + dt * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6))

    k7 = f.f1(duprev, ku, p, t + dt * c6)
    u = uprev + dt * (duprev + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6))
    du = duprev + dt * (bp1 * k1 + bp3 * k3 + bp4 * k4 + bp5 * k5 + bp6 * k6 + bp7 * k7)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq *
               (btilde1 * k1 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5 + btilde6 * k6)
        duhat = dt * (bptilde1 * k1 + bptilde3 * k3 + bptilde4 * k4 + bptilde5 * k5 +
                 bptilde6 * k6 + bptilde7 * k7)
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

@muladd function perform_step!(integrator, cache::ERKN7Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    @unpack tmp, atmp, fsalfirst, k2, k3, k4, k5, k6, k7, k, utilde = cache
    @unpack c1, c2, c3, c4, c5, c6, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a73, a74, a75, a76, b1, b3, b4, b5, b6, bp1, bp3, bp4, bp5, bp6, bp7, btilde1, btilde3, btilde4, btilde5, btilde6, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6, bptilde7 = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast=false ku=uprev + dt * (c1 * duprev + dt * a21 * k1)

    f.f1(k2, duprev, ku, p, t + dt * c1)
    @.. broadcast=false ku=uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    f.f1(k3, duprev, ku, p, t + dt * c2)
    @.. broadcast=false ku=uprev +
                           dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    f.f1(k4, duprev, ku, p, t + dt * c3)
    @.. broadcast=false ku=uprev +
                           dt *
                           (c4 * duprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))

    f.f1(k5, duprev, ku, p, t + dt * c4)
    @.. broadcast=false ku=uprev +
                           dt * (c5 * duprev +
                            dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5))

    f.f1(k6, duprev, ku, p, t + dt * c5)
    @.. broadcast=false ku=uprev +
                           dt * (c6 * duprev +
                            dt * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6))

    f.f1(k7, duprev, ku, p, t + dt * c6)
    @tight_loop_macros for i in uidx
        @inbounds u[i] = uprev[i] +
                         dt * (duprev[i] +
                          dt *
                          (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]))
        @inbounds du[i] = duprev[i] +
                          dt * (bp1 * k1[i] + bp3 * k3[i] + bp4 * k4[i] + bp5 * k5[i] +
                           bp6 * k6[i] + bp7 * k7[i])
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    integrator.stats.nf += 4
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @tight_loop_macros for i in uidx
            @inbounds uhat[i] = dtsq *
                                (btilde1 * k1[i] + btilde3 * k3[i] + btilde4 * k4[i] +
                                 btilde5 * k5[i] + btilde6 * k6[i])
            @inbounds duhat[i] = dt *
                                 (bptilde1 * k1[i] + bptilde3 * k3[i] + bptilde4 * k4[i] +
                                  bptilde5 * k5[i] + bptilde6 * k6[i] + bptilde7 * k7[i])
        end
        calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

function initialize!(integrator, cache::RKN4Cache)
    @unpack fsalfirst, k = cache
    duprev, uprev = integrator.uprev.x
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = cache.k₃
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.k[1].x[1], duprev, uprev, integrator.p, integrator.t)
    integrator.f.f2(integrator.k[1].x[2], duprev, uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
end

@muladd function perform_step!(integrator, cache::RKN4constantCache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    
    #define dt values
    halfdt = dt/2
    dtsq = dt^2
    eightdtsq = dtsq/8
    halfdtsq = dtsq/2
    sixthdtsq = dtsq/6
    sixthdt = dt/6
    ttmp = t + halfdt

    #perform operations to find k values
    k₁ = integrator.fsalfirst.x[1]
    ku = uprev + halfdt * duprev + eightdtsq * k₁
    kdu = duprev + halfdt * k₁

    k₂ = f.f1(kdu, ku, p, ttmp)
    ku = uprev + dt * duprev + halfdtsq * k₂
    kdu = duprev + dt * k₂

    k₃ = f.f1(kdu, ku, p, t + dt)
    ku = uprev + dt * duprev + eightdtsq * k₃
    kdu = duprev + dt * k₃

    #perform final calculations to determine new y and y'.
    u = uprev + sixthdtsq* (1*k₁ + 2*k₂ + 0*k₃) + dt * duprev
    du = duprev + sixthdt * (1*k₁ + 4*k₂ + 1*k₃)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.stats.nf += 3
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::RKN4Cache, repeat_step = false)
    @unpack t, dt, f, p = integrator
    duprev, uprev = integrator.uprev.x
    
    #define dt values
    halfdt = dt/2
    dtsq = dt^2
    eightdtsq = dtsq/8
    halfdtsq = dtsq/2
    sixthdtsq = dtsq/6
    sixthdt = dt/6
    ttmp = t + halfdt

    #perform operations to find k values
    k₁ = integrator.fsalfirst.x[1]
    ku = uprev + halfdt * duprev + eightdtsq * k₁
    kdu = duprev + halfdt * k₁

    f.f1(k₂, kdu, ku, p, ttmp)
    ku = uprev + dt * duprev + halfdtsq * k₂
    kdu = duprev + dt * k₂

    f.f1(k₃, kdu, ku, p, t + dt)
    ku = uprev + dt * duprev + eightdtsq * k₃
    kdu = duprev + dt * k₃

    #perform final calculations to determine new y and y'.
    u = uprev + sixthdtsq* (1*k₁ + 2*k₂ + 0*k₃) + dt * duprev
    du = duprev + sixthdt * (1*k₁ + 4*k₂ + 1*k₃)

    integrator.stats.nf += 3
    integrator.stats.nf2 += 1
end