## y'' = f(t, y, y')
## y(t₀) = y₀; y'(t₀) = y'₀
## kᵢ' = f(t₀+cᵢh, y₀+cᵢhy'₀+h²∑āᵢⱼk'ⱼ, y'₀+h∑aᵢⱼk'ⱼ)
## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
## y'₁ = y'₀ + h∑bᵢk'ᵢ

const NystromCCDefaultInitialization = Union{
    Nystrom4VelocityIndependentConstantCache,
    IRKN3ConstantCache, IRKN4ConstantCache,
}

function initialize!(integrator, cache::NystromCCDefaultInitialization)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    duprev, uprev = integrator.uprev.x
    kdu = integrator.f.f1(duprev, uprev, integrator.p, integrator.t)
    ku = integrator.f.f2(duprev, uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    return integrator.fsalfirst = ArrayPartition((kdu, ku))
end

const NystromDefaultInitialization = Union{
    Nystrom4VelocityIndependentCache,
    IRKN3Cache, IRKN4Cache,
}

function initialize!(integrator, cache::NystromDefaultInitialization)
    (; fsalfirst, k) = cache
    duprev, uprev = integrator.uprev.x
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.k[1].x[1], duprev, uprev, integrator.p, integrator.t)
    integrator.f.f2(integrator.k[1].x[2], duprev, uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return integrator.stats.nf2 += 1
end

## Nystrom4VelocityIndependent perform_step! kept for IRKN3/IRKN4 bootstrap use

@muladd function perform_step!(
        integrator, cache::Nystrom4VelocityIndependentConstantCache,
        repeat_step = false
    )
    (; t, dt, f, p) = integrator
    duprev, uprev = integrator.uprev.x
    k₁ = integrator.fsalfirst.x[1]
    halfdt = dt / 2
    dtsq = dt^2
    eighth_dtsq = dtsq / 8
    half_dtsq = dtsq / 2
    ttmp = t + halfdt

    ku = uprev + halfdt * duprev + eighth_dtsq * k₁

    k₂ = f.f1(duprev, ku, p, ttmp)
    ku = uprev + dt * duprev + half_dtsq * k₂

    k₃ = f.f1(duprev, ku, p, t + dt)
    u = uprev + (dtsq / 6) * (k₁ + 2 * k₂) + dt * duprev
    du = duprev + (dt / 6) * (k₁ + k₃ + 4 * k₂)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(
        integrator, cache::Nystrom4VelocityIndependentCache,
        repeat_step = false
    )
    (; t, dt, f, p) = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    (; tmp, fsalfirst, k₂, k₃, k) = cache
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    k₁ = integrator.fsalfirst.x[1]
    halfdt = dt / 2
    dtsq = dt^2
    eighth_dtsq = dtsq / 8
    half_dtsq = dtsq / 2
    ttmp = t + halfdt

    @.. broadcast = false ku = uprev + halfdt * duprev + eighth_dtsq * k₁

    f.f1(k₂, duprev, ku, p, ttmp)
    @.. broadcast = false ku = uprev + dt * duprev + half_dtsq * k₂

    f.f1(k₃, duprev, ku, p, t + dt)
    @.. broadcast = false u = uprev + (dtsq / 6) * (k₁ + 2 * k₂) + dt * duprev
    @.. broadcast = false du = duprev + (dt / 6) * (k₁ + k₃ + 4 * k₂)

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    integrator.stats.nf2 += 1
end

@muladd function perform_step!(integrator, cache::IRKN3ConstantCache, repeat_step = false)
    (; t, dt, k, tprev, f, p) = integrator
    duprev, uprev = integrator.uprev.x
    duprev2, uprev2 = integrator.uprev2.x
    (; bconst1, bconst2, c1, a21, b1, b2, bbar1, bbar2, k₂) = cache
    k₁ = integrator.fsalfirst
    # if there's a discontinuity or the solver is in the first step
    if integrator.iter < 2 && !integrator.u_modified
        perform_step!(integrator, Nystrom4VelocityIndependentConstantCache())
        k = integrator.fsallast
        k1cache = ArrayPartition((k.x[1], f.f1(duprev, uprev, p, t + c1 * dt)))
        kdu = uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[1])
        k₂.x[1] = f.f1(duprev, kdu, p, t + c1 * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    else
        kdu = uprev2 + dt * (c1 * duprev2 + dt * a21 * k1cache.x[1])
        ku = uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[2])

        k₂x1 = f.f1(duprev, ku, p, t + c1 * dt)
        du = duprev +
            dt * (b1 * k1cache.x[1] + bbar1 * k1cache.x[1] + b2 * (k₂x1 - k₂.x[1]))
        u = uprev + bconst1 * dt * duprev +
            dt * (bconst2 * duprev2 + dt * bbar2 * (k₂x1 - k₂.x[1]))

        integrator.fsallast = ArrayPartition(
            (
                f.f1(du, u, p, t + dt),
                f.f2(du, u, p, t + dt),
            )
        )
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
        integrator.stats.nf2 += 1
        copyto!(k₂.x[1], k₂.x[2])
        k1cache = ArrayPartition((k1cache.x[1], k.x[2]))
    end # end if
end

@muladd function perform_step!(integrator, cache::IRKN3Cache, repeat_step = false)
    (; t, dt, k, tprev, f, p) = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    duprev2, uprev2 = integrator.uprev2.x
    uidx = eachindex(integrator.uprev.x[1])
    (; tmp, fsalfirst, k₂, k) = cache
    (; bconst1, bconst2, c1, a21, b1, b2, bbar1, bbar2) = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    k1cache = cache.tmp2
    k₁ = fsalfirst
    # if there's a discontinuity or the solver is in the first step
    if integrator.iter < 2 && !integrator.u_modified
        perform_step!(integrator, integrator.cache.onestep_cache)
        copyto!(k1cache.x[1], k.x[1])
        f.f1(k1cache.x[2], duprev, uprev, p, t + c1 * dt)
        @.. broadcast = false kdu = uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[2])
        f.f1(k₂.x[1], duprev, kdu, p, t + c1 * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    else
        @.. broadcast = false kdu = uprev2 + dt * (c1 * duprev2 + dt * a21 * k1cache.x[1])
        @.. broadcast = false ku = uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[2])

        f.f1(k₂.x[2], duprev, ku, p, t + c1 * dt)
        @tight_loop_macros for i in uidx
            @inbounds u[i] = uprev[i] + bconst1 * dt * duprev[i] +
                dt *
                (bconst2 * duprev2[i] + dt * bbar2 * (k₂.x[2][i] - k₂.x[1][i]))
            @inbounds du[i] = duprev[i] +
                dt * (
                b1 * k1cache.x[1][i] + bbar1 * k1cache.x[2][i] +
                    b2 * (k₂.x[2][i] - k₂.x[1][i])
            )
        end
        f.f1(k.x[1], du, u, p, t + dt)
        f.f2(k.x[2], du, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
        integrator.stats.nf2 += 1
        copyto!(k₂.x[1], k₂.x[2])
        copyto!(k1cache.x[2], k1cache.x[1])
        copyto!(k1cache.x[1], k.x[1])
    end # end if
end

@muladd function perform_step!(integrator, cache::IRKN4Cache, repeat_step = false)
    (; t, dt, k, tprev, f, p) = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    duprev2, uprev2 = integrator.uprev2.x
    uidx = eachindex(integrator.uprev.x[1])
    (; tmp, tmp2, fsalfirst, k₂, k₃, k) = cache
    (; bconst1, bconst2, c1, c2, a21, a32, b1, b2, b3, bbar1, bbar2, bbar3) = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    k1cache = integrator.cache.tmp2
    k₁ = fsalfirst
    # if there's a discontinuity or the solver is in the first step
    if integrator.iter < 2 && !integrator.u_modified
        perform_step!(integrator, integrator.cache.onestep_cache)
        copyto!(k1cache.x[1], k.x[1])
        f.f1(k1cache.x[2], duprev, uprev, p, t + c1 * dt)
        @.. broadcast = false kdu = uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[1])
        f.f1(k₂.x[1], duprev, kdu, p, t + c1 * dt)
        @.. broadcast = false kdu = uprev + dt * (c2 * duprev + dt * a32 * k1cache.x[2])
        f.f1(k₃.x[1], duprev, kdu, p, t + c1 * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
    else
        @.. broadcast = false ku = uprev + dt * (c1 * duprev + dt * a21 * k1cache.x[1])
        @.. broadcast = false kdu = uprev2 + dt * (c1 * duprev2 + dt * a21 * k1cache.x[2])

        f.f1(k₂.x[2], duprev, ku, p, t + c1 * dt)
        @.. broadcast = false ku = uprev + dt * (c2 * duprev + dt * a32 * k₂.x[2])
        @.. broadcast = false kdu = uprev2 + dt * (c2 * duprev2 + dt * a32 * k₂.x[1])

        f.f1(k₃.x[2], duprev, ku, p, t + c2 * dt)
        @tight_loop_macros for i in uidx
            @inbounds u[i] = uprev[i] + dt * bconst1 * duprev[i] +
                dt * (
                bconst2 * duprev2[i] +
                    dt * (
                    bbar2 * (k₂.x[2][i] - k₂.x[1][i]) +
                        bbar3 * (k₃.x[2][i] - k₃.x[1][i])
                )
            )
            @inbounds du[i] = duprev[i] +
                dt * (
                b1 * k1cache.x[1][i] + bbar1 * k1cache.x[2][i] +
                    b2 * (k₂.x[2][i] - k₂.x[1][i]) +
                    b3 * (k₃.x[2][i] - k₃.x[1][i])
            )
        end
        f.f1(k.x[1], du, u, p, t + dt)
        f.f2(k.x[2], du, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
        integrator.stats.nf2 += 1
        copyto!(k₂.x[1], k₂.x[2])
        copyto!(k₃.x[1], k₃.x[2])
        copyto!(k1cache.x[2], k1cache.x[1])
        copyto!(k1cache.x[1], k.x[1])
    end # end if
end

function initialize!(integrator, cache::DPRKN6ConstantCache)
    duprev, uprev = integrator.uprev.x
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    kdu = integrator.f.f1(duprev, uprev, integrator.p, integrator.t)
    ku = integrator.f.f2(duprev, uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.fsalfirst = ArrayPartition((kdu, ku))
    integrator.fsallast = zero(integrator.fsalfirst)

    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::DPRKN6ConstantCache, repeat_step = false)
    (; t, dt, f, p) = integrator
    duprev, uprev = integrator.uprev.x
    (; c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1, bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4, btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6) = cache
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

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    integrator.k[1] = ArrayPartition(integrator.fsalfirst.x[1], k2)
    integrator.k[2] = ArrayPartition(k3, k4)
    integrator.k[3] = ArrayPartition(k5, k6)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    integrator.stats.nf2 += 1

    if integrator.opts.adaptive
        dtsq = dt^2
        uhat = dtsq *
            (btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4 + btilde5 * k5)
        duhat = dt * (
            bptilde1 * k1 + bptilde3 * k3 + bptilde4 * k4 + bptilde5 * k5 +
                bptilde6 * k6
        )
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(
            utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

function initialize!(integrator, cache::DPRKN6Cache)
    (; fsalfirst, k) = cache
    duprev, uprev = integrator.uprev.x
    integrator.kshortsize = 3
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = ArrayPartition(cache.fsalfirst.x[1], cache.k2)
    integrator.k[2] = ArrayPartition(cache.k3, cache.k4)
    integrator.k[3] = ArrayPartition(cache.k5, cache.k6)
    integrator.f.f1(integrator.fsalfirst.x[1], duprev, uprev, integrator.p, integrator.t)
    integrator.f.f2(integrator.fsalfirst.x[2], duprev, uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return integrator.stats.nf2 += 1
end

@muladd function perform_step!(integrator, cache::DPRKN6Cache, repeat_step = false)
    (; t, dt, f, p) = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    (; tmp, atmp, fsalfirst, k2, k3, k4, k5, k6, k, utilde) = cache
    (; c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1, bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4, btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6) = cache.tab
    kdu, ku = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    uidx = eachindex(integrator.uprev.x[2])
    k1 = integrator.fsalfirst.x[1]

    @.. broadcast = false ku = uprev + dt * (c1 * duprev + dt * a21 * k1)

    f.f1(k2, duprev, ku, p, t + dt * c1)
    @.. broadcast = false ku = uprev + dt * (c2 * duprev + dt * (a31 * k1 + a32 * k2))

    f.f1(k3, duprev, ku, p, t + dt * c2)
    @.. broadcast = false ku = uprev +
        dt * (c3 * duprev + dt * (a41 * k1 + a42 * k2 + a43 * k3))

    f.f1(k4, duprev, ku, p, t + dt * c3)
    @.. broadcast = false ku = uprev +
        dt *
        (c4 * duprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4))

    f.f1(k5, duprev, ku, p, t + dt * c4)
    @.. broadcast = false ku = uprev +
        dt *
        (c5 * duprev + dt * (a61 * k1 + a63 * k3 + a64 * k4 + a65 * k5)) # no a62

    f.f1(k6, duprev, ku, p, t + dt * c5)

    @.. broadcast = false u = uprev +
        dt * (duprev + dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5)) # b1 -- b5, no b2
    @.. broadcast = false du = duprev +
        dt * (bp1 * k1 + bp3 * k3 + bp4 * k4 + bp5 * k5 + bp6 * k6) # bp1 -- bp6, no bp2

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 6)
    integrator.stats.nf2 += 1
    if integrator.opts.adaptive
        duhat, uhat = utilde.x
        dtsq = dt^2
        @.. broadcast = false uhat = dtsq * (
            btilde1 * k1 + btilde2 * k2 + btilde3 * k3 +
                btilde4 * k4 + btilde5 * k5
        )
        @.. broadcast = false duhat = dt * (
            bptilde1 * k1 + bptilde3 * k3 + bptilde4 * k4 +
                bptilde5 * k5 + bptilde6 * k6
        )
        calculate_residuals!(
            atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

