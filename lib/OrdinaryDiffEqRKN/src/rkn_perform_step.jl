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
    if integrator.iter < 2 && !integrator.derivative_discontinuity
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
    if integrator.iter < 2 && !integrator.derivative_discontinuity
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
    if integrator.iter < 2 && !integrator.derivative_discontinuity
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

################################################################################
# Generic Nyström velocity-independent / velocity-dependent perform_step!
################################################################################

## Generic Nyström velocity-independent perform_step!
## Solves: y'' = f(t, y) where f is velocity-independent
## kᵢ = f1(duprev, yᵢ, p, t + cᵢ*dt)   (duprev constant throughout)
## yᵢ = y₀ + cᵢ*h*y'₀ + h²*Σⱼ<ᵢ aᵢⱼ*kⱼ
## y₁ = y₀ + h*y'₀ + h²*Σᵢ bᵢ*kᵢ
## y'₁ = y'₀ + h*Σᵢ bpᵢ*kᵢ

# Dense-output bookkeeping for velocity-independent methods. Standard methods use
# 2-point Hermite (fsalfirst/fsallast); DPRKN6 stores its six stage derivatives for
# its specialized "free" 6th-order interpolant (see interpolants.jl). The required
# `integrator.k` length lives on the tableau as `kshortsize`.

# Constant cache: pack the just-computed stages `ks` (ks[1]=k1 … ks[6]=k6) into integrator.k.
function _store_vi_k!(integrator, ::Any, ks)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return
end
function _store_vi_k!(integrator, ::DPRKN6Tableau, ks)
    integrator.k[1] = ArrayPartition(ks[1], ks[2])
    integrator.k[2] = ArrayPartition(ks[3], ks[4])
    integrator.k[3] = ArrayPartition(ks[5], ks[6])
    return
end

# Mutable cache: point integrator.k at the persistent stage buffers so the interpolant
# sees the current stages without re-storing each step.
function _init_vi_k_mut!(integrator, cache, tab::Any)
    integrator.kshortsize = tab.kshortsize
    resize!(integrator.k, tab.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return
end
function _init_vi_k_mut!(integrator, cache, tab::DPRKN6Tableau)
    integrator.kshortsize = tab.kshortsize
    resize!(integrator.k, tab.kshortsize)
    integrator.k[1] = ArrayPartition(integrator.fsalfirst.x[1], cache.ks[1])
    integrator.k[2] = ArrayPartition(cache.ks[2], cache.ks[3])
    integrator.k[3] = ArrayPartition(cache.ks[4], cache.ks[5])
    return
end

function initialize!(integrator, cache::NystromVIConstantCache)
    integrator.kshortsize = cache.tab.kshortsize
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    duprev, uprev = integrator.uprev.x
    kdu = integrator.f.f1(duprev, uprev, integrator.p, integrator.t)
    ku = integrator.f.f2(duprev, uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.fsalfirst = ArrayPartition((kdu, ku))
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    integrator.k[integrator.kshortsize] = integrator.fsallast
    return
end

function initialize!(integrator, cache::NystromVICache)
    duprev, uprev = integrator.uprev.x
    integrator.f.f1(integrator.fsalfirst.x[1], duprev, uprev, integrator.p, integrator.t)
    integrator.f.f2(integrator.fsalfirst.x[2], duprev, uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    _init_vi_k_mut!(integrator, cache, cache.tab)
    return
end

@muladd function perform_step!(
        integrator, cache::NystromVIConstantCache, repeat_step = false
    )
    (; t, dt, f, p) = integrator
    duprev, uprev = integrator.uprev.x
    (; tab) = cache
    (; a, b, bp, btilde, bptilde, c, pos_only_error) = tab
    k1 = integrator.fsalfirst.x[1]
    nstages = length(b)
    dtsq = dt^2
    ks = Vector{typeof(k1)}(undef, nstages)
    ks[1] = k1

    # ---------------- Stages 2..nstages (explicit ladder, max nstages = 17 for DPRKN12) ----------------
    if nstages >= 2
        ku = uprev + dt * (c[1] * duprev + dt * (a[2, 1] * ks[1]))
        ks[2] = f.f1(duprev, ku, p, t + dt * c[1])
    end
    if nstages >= 3
        ku = uprev + dt * (c[2] * duprev + dt * (a[3, 1] * ks[1] + a[3, 2] * ks[2]))
        ks[3] = f.f1(duprev, ku, p, t + dt * c[2])
    end
    if nstages >= 4
        ku = uprev + dt * (c[3] * duprev + dt * (a[4, 1] * ks[1] + a[4, 2] * ks[2] + a[4, 3] * ks[3]))
        ks[4] = f.f1(duprev, ku, p, t + dt * c[3])
    end
    if nstages >= 5
        ku = uprev + dt * (c[4] * duprev + dt * (a[5, 1] * ks[1] + a[5, 2] * ks[2] + a[5, 3] * ks[3] + a[5, 4] * ks[4]))
        ks[5] = f.f1(duprev, ku, p, t + dt * c[4])
    end
    if nstages >= 6
        ku = uprev + dt * (c[5] * duprev + dt * (a[6, 1] * ks[1] + a[6, 2] * ks[2] + a[6, 3] * ks[3] + a[6, 4] * ks[4] + a[6, 5] * ks[5]))
        ks[6] = f.f1(duprev, ku, p, t + dt * c[5])
    end
    if nstages >= 7
        ku = uprev + dt * (c[6] * duprev + dt * (a[7, 1] * ks[1] + a[7, 2] * ks[2] + a[7, 3] * ks[3] + a[7, 4] * ks[4] + a[7, 5] * ks[5] + a[7, 6] * ks[6]))
        ks[7] = f.f1(duprev, ku, p, t + dt * c[6])
    end
    if nstages >= 8
        ku = uprev + dt * (c[7] * duprev + dt * (a[8, 1] * ks[1] + a[8, 2] * ks[2] + a[8, 3] * ks[3] + a[8, 4] * ks[4] + a[8, 5] * ks[5] + a[8, 6] * ks[6] + a[8, 7] * ks[7]))
        ks[8] = f.f1(duprev, ku, p, t + dt * c[7])
    end
    if nstages >= 9
        ku = uprev + dt * (c[8] * duprev + dt * (a[9, 1] * ks[1] + a[9, 2] * ks[2] + a[9, 3] * ks[3] + a[9, 4] * ks[4] + a[9, 5] * ks[5] + a[9, 6] * ks[6] + a[9, 7] * ks[7] + a[9, 8] * ks[8]))
        ks[9] = f.f1(duprev, ku, p, t + dt * c[8])
    end
    if nstages >= 10
        ku = uprev + dt * (c[9] * duprev + dt * (a[10, 1] * ks[1] + a[10, 2] * ks[2] + a[10, 3] * ks[3] + a[10, 4] * ks[4] + a[10, 5] * ks[5] + a[10, 6] * ks[6] + a[10, 7] * ks[7] + a[10, 8] * ks[8] + a[10, 9] * ks[9]))
        ks[10] = f.f1(duprev, ku, p, t + dt * c[9])
    end
    if nstages >= 11
        ku = uprev + dt * (c[10] * duprev + dt * (a[11, 1] * ks[1] + a[11, 2] * ks[2] + a[11, 3] * ks[3] + a[11, 4] * ks[4] + a[11, 5] * ks[5] + a[11, 6] * ks[6] + a[11, 7] * ks[7] + a[11, 8] * ks[8] + a[11, 9] * ks[9] + a[11, 10] * ks[10]))
        ks[11] = f.f1(duprev, ku, p, t + dt * c[10])
    end
    if nstages >= 12
        ku = uprev + dt * (c[11] * duprev + dt * (a[12, 1] * ks[1] + a[12, 2] * ks[2] + a[12, 3] * ks[3] + a[12, 4] * ks[4] + a[12, 5] * ks[5] + a[12, 6] * ks[6] + a[12, 7] * ks[7] + a[12, 8] * ks[8] + a[12, 9] * ks[9] + a[12, 10] * ks[10] + a[12, 11] * ks[11]))
        ks[12] = f.f1(duprev, ku, p, t + dt * c[11])
    end
    if nstages >= 13
        ku = uprev + dt * (c[12] * duprev + dt * (a[13, 1] * ks[1] + a[13, 2] * ks[2] + a[13, 3] * ks[3] + a[13, 4] * ks[4] + a[13, 5] * ks[5] + a[13, 6] * ks[6] + a[13, 7] * ks[7] + a[13, 8] * ks[8] + a[13, 9] * ks[9] + a[13, 10] * ks[10] + a[13, 11] * ks[11] + a[13, 12] * ks[12]))
        ks[13] = f.f1(duprev, ku, p, t + dt * c[12])
    end
    if nstages >= 14
        ku = uprev + dt * (c[13] * duprev + dt * (a[14, 1] * ks[1] + a[14, 2] * ks[2] + a[14, 3] * ks[3] + a[14, 4] * ks[4] + a[14, 5] * ks[5] + a[14, 6] * ks[6] + a[14, 7] * ks[7] + a[14, 8] * ks[8] + a[14, 9] * ks[9] + a[14, 10] * ks[10] + a[14, 11] * ks[11] + a[14, 12] * ks[12] + a[14, 13] * ks[13]))
        ks[14] = f.f1(duprev, ku, p, t + dt * c[13])
    end
    if nstages >= 15
        ku = uprev + dt * (c[14] * duprev + dt * (a[15, 1] * ks[1] + a[15, 2] * ks[2] + a[15, 3] * ks[3] + a[15, 4] * ks[4] + a[15, 5] * ks[5] + a[15, 6] * ks[6] + a[15, 7] * ks[7] + a[15, 8] * ks[8] + a[15, 9] * ks[9] + a[15, 10] * ks[10] + a[15, 11] * ks[11] + a[15, 12] * ks[12] + a[15, 13] * ks[13] + a[15, 14] * ks[14]))
        ks[15] = f.f1(duprev, ku, p, t + dt * c[14])
    end
    if nstages >= 16
        ku = uprev + dt * (c[15] * duprev + dt * (a[16, 1] * ks[1] + a[16, 2] * ks[2] + a[16, 3] * ks[3] + a[16, 4] * ks[4] + a[16, 5] * ks[5] + a[16, 6] * ks[6] + a[16, 7] * ks[7] + a[16, 8] * ks[8] + a[16, 9] * ks[9] + a[16, 10] * ks[10] + a[16, 11] * ks[11] + a[16, 12] * ks[12] + a[16, 13] * ks[13] + a[16, 14] * ks[14] + a[16, 15] * ks[15]))
        ks[16] = f.f1(duprev, ku, p, t + dt * c[15])
    end
    if nstages >= 17
        ku = uprev + dt * (c[16] * duprev + dt * (a[17, 1] * ks[1] + a[17, 2] * ks[2] + a[17, 3] * ks[3] + a[17, 4] * ks[4] + a[17, 5] * ks[5] + a[17, 6] * ks[6] + a[17, 7] * ks[7] + a[17, 8] * ks[8] + a[17, 9] * ks[9] + a[17, 10] * ks[10] + a[17, 11] * ks[11] + a[17, 12] * ks[12] + a[17, 13] * ks[13] + a[17, 14] * ks[14] + a[17, 15] * ks[15] + a[17, 16] * ks[16]))
        ks[17] = f.f1(duprev, ku, p, t + dt * c[16])
    end

    # ---------------- Final update u, du ----------------
    u = uprev + dt * duprev
    du = duprev
    for i in 1:nstages
        u = u + dtsq * b[i] * ks[i]
        du = du + dt * bp[i] * ks[i]
    end

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, nstages)
    integrator.stats.nf2 += 1
    _store_vi_k!(integrator, tab, ks)

    if integrator.opts.adaptive && !isempty(btilde)
        uhat = zero(uprev)
        for i in 1:nstages
            uhat = uhat + dtsq * btilde[i] * ks[i]
        end
        if pos_only_error
            atmp = calculate_residuals(
                uhat, integrator.uprev.x[2], integrator.u.x[2],
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
        else
            duhat = zero(duprev)
            for i in 1:nstages
                duhat = duhat + dt * bptilde[i] * ks[i]
            end
            utilde = ArrayPartition((duhat, uhat))
            atmp = calculate_residuals(
                utilde, integrator.uprev, integrator.u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
        end
    end
end

@muladd function perform_step!(
        integrator, cache::NystromVICache, repeat_step = false
    )
    (; t, dt, f, p) = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    (; ks, k, utilde, tmp, atmp, tab) = cache
    (; a, b, bp, btilde, bptilde, c, pos_only_error) = tab
    ku = tmp.x[2]
    k1 = integrator.fsalfirst.x[1]
    nstages = length(b)
    dtsq = dt^2

    # ---------------- Stages 2..nstages (explicit ladder; ks[i-1] holds the i-th stage) ----------------
    if nstages >= 2
        @.. broadcast = false ku = uprev + dt * (c[1] * duprev + dt * (a[2, 1] * k1))
        f.f1(ks[1], duprev, ku, p, t + dt * c[1])
    end
    if nstages >= 3
        @.. broadcast = false ku = uprev + dt * (c[2] * duprev + dt * (a[3, 1] * k1 + a[3, 2] * ks[1]))
        f.f1(ks[2], duprev, ku, p, t + dt * c[2])
    end
    if nstages >= 4
        @.. broadcast = false ku = uprev + dt * (c[3] * duprev + dt * (a[4, 1] * k1 + a[4, 2] * ks[1] + a[4, 3] * ks[2]))
        f.f1(ks[3], duprev, ku, p, t + dt * c[3])
    end
    if nstages >= 5
        @.. broadcast = false ku = uprev + dt * (c[4] * duprev + dt * (a[5, 1] * k1 + a[5, 2] * ks[1] + a[5, 3] * ks[2] + a[5, 4] * ks[3]))
        f.f1(ks[4], duprev, ku, p, t + dt * c[4])
    end
    if nstages >= 6
        @.. broadcast = false ku = uprev + dt * (c[5] * duprev + dt * (a[6, 1] * k1 + a[6, 2] * ks[1] + a[6, 3] * ks[2] + a[6, 4] * ks[3] + a[6, 5] * ks[4]))
        f.f1(ks[5], duprev, ku, p, t + dt * c[5])
    end
    if nstages >= 7
        @.. broadcast = false ku = uprev + dt * (c[6] * duprev + dt * (a[7, 1] * k1 + a[7, 2] * ks[1] + a[7, 3] * ks[2] + a[7, 4] * ks[3] + a[7, 5] * ks[4] + a[7, 6] * ks[5]))
        f.f1(ks[6], duprev, ku, p, t + dt * c[6])
    end
    if nstages >= 8
        @.. broadcast = false ku = uprev + dt * (c[7] * duprev + dt * (a[8, 1] * k1 + a[8, 2] * ks[1] + a[8, 3] * ks[2] + a[8, 4] * ks[3] + a[8, 5] * ks[4] + a[8, 6] * ks[5] + a[8, 7] * ks[6]))
        f.f1(ks[7], duprev, ku, p, t + dt * c[7])
    end
    if nstages >= 9
        @.. broadcast = false ku = uprev + dt * (c[8] * duprev + dt * (a[9, 1] * k1 + a[9, 2] * ks[1] + a[9, 3] * ks[2] + a[9, 4] * ks[3] + a[9, 5] * ks[4] + a[9, 6] * ks[5] + a[9, 7] * ks[6] + a[9, 8] * ks[7]))
        f.f1(ks[8], duprev, ku, p, t + dt * c[8])
    end
    if nstages >= 10
        @.. broadcast = false ku = uprev + dt * (c[9] * duprev + dt * (a[10, 1] * k1 + a[10, 2] * ks[1] + a[10, 3] * ks[2] + a[10, 4] * ks[3] + a[10, 5] * ks[4] + a[10, 6] * ks[5] + a[10, 7] * ks[6] + a[10, 8] * ks[7] + a[10, 9] * ks[8]))
        f.f1(ks[9], duprev, ku, p, t + dt * c[9])
    end
    if nstages >= 11
        @.. broadcast = false ku = uprev + dt * (c[10] * duprev + dt * (a[11, 1] * k1 + a[11, 2] * ks[1] + a[11, 3] * ks[2] + a[11, 4] * ks[3] + a[11, 5] * ks[4] + a[11, 6] * ks[5] + a[11, 7] * ks[6] + a[11, 8] * ks[7] + a[11, 9] * ks[8] + a[11, 10] * ks[9]))
        f.f1(ks[10], duprev, ku, p, t + dt * c[10])
    end
    if nstages >= 12
        @.. broadcast = false ku = uprev + dt * (c[11] * duprev + dt * (a[12, 1] * k1 + a[12, 2] * ks[1] + a[12, 3] * ks[2] + a[12, 4] * ks[3] + a[12, 5] * ks[4] + a[12, 6] * ks[5] + a[12, 7] * ks[6] + a[12, 8] * ks[7] + a[12, 9] * ks[8] + a[12, 10] * ks[9] + a[12, 11] * ks[10]))
        f.f1(ks[11], duprev, ku, p, t + dt * c[11])
    end
    if nstages >= 13
        @.. broadcast = false ku = uprev + dt * (c[12] * duprev + dt * (a[13, 1] * k1 + a[13, 2] * ks[1] + a[13, 3] * ks[2] + a[13, 4] * ks[3] + a[13, 5] * ks[4] + a[13, 6] * ks[5] + a[13, 7] * ks[6] + a[13, 8] * ks[7] + a[13, 9] * ks[8] + a[13, 10] * ks[9] + a[13, 11] * ks[10] + a[13, 12] * ks[11]))
        f.f1(ks[12], duprev, ku, p, t + dt * c[12])
    end
    if nstages >= 14
        @.. broadcast = false ku = uprev + dt * (c[13] * duprev + dt * (a[14, 1] * k1 + a[14, 2] * ks[1] + a[14, 3] * ks[2] + a[14, 4] * ks[3] + a[14, 5] * ks[4] + a[14, 6] * ks[5] + a[14, 7] * ks[6] + a[14, 8] * ks[7] + a[14, 9] * ks[8] + a[14, 10] * ks[9] + a[14, 11] * ks[10] + a[14, 12] * ks[11] + a[14, 13] * ks[12]))
        f.f1(ks[13], duprev, ku, p, t + dt * c[13])
    end
    if nstages >= 15
        @.. broadcast = false ku = uprev + dt * (c[14] * duprev + dt * (a[15, 1] * k1 + a[15, 2] * ks[1] + a[15, 3] * ks[2] + a[15, 4] * ks[3] + a[15, 5] * ks[4] + a[15, 6] * ks[5] + a[15, 7] * ks[6] + a[15, 8] * ks[7] + a[15, 9] * ks[8] + a[15, 10] * ks[9] + a[15, 11] * ks[10] + a[15, 12] * ks[11] + a[15, 13] * ks[12] + a[15, 14] * ks[13]))
        f.f1(ks[14], duprev, ku, p, t + dt * c[14])
    end
    if nstages >= 16
        @.. broadcast = false ku = uprev + dt * (c[15] * duprev + dt * (a[16, 1] * k1 + a[16, 2] * ks[1] + a[16, 3] * ks[2] + a[16, 4] * ks[3] + a[16, 5] * ks[4] + a[16, 6] * ks[5] + a[16, 7] * ks[6] + a[16, 8] * ks[7] + a[16, 9] * ks[8] + a[16, 10] * ks[9] + a[16, 11] * ks[10] + a[16, 12] * ks[11] + a[16, 13] * ks[12] + a[16, 14] * ks[13] + a[16, 15] * ks[14]))
        f.f1(ks[15], duprev, ku, p, t + dt * c[15])
    end
    if nstages >= 17
        @.. broadcast = false ku = uprev + dt * (c[16] * duprev + dt * (a[17, 1] * k1 + a[17, 2] * ks[1] + a[17, 3] * ks[2] + a[17, 4] * ks[3] + a[17, 5] * ks[4] + a[17, 6] * ks[5] + a[17, 7] * ks[6] + a[17, 8] * ks[7] + a[17, 9] * ks[8] + a[17, 10] * ks[9] + a[17, 11] * ks[10] + a[17, 12] * ks[11] + a[17, 13] * ks[12] + a[17, 14] * ks[13] + a[17, 15] * ks[14] + a[17, 16] * ks[15]))
        f.f1(ks[16], duprev, ku, p, t + dt * c[16])
    end

    # ---------------- Final update u, du ----------------
    @.. broadcast = false u = uprev + dt * duprev
    @.. broadcast = false du = duprev
    @.. broadcast = false u = u + dtsq * b[1] * k1
    @.. broadcast = false du = du + dt * bp[1] * k1
    for i in 2:nstages
        ki = ks[i - 1]
        @.. broadcast = false u = u + dtsq * b[i] * ki
        @.. broadcast = false du = du + dt * bp[i] * ki
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, nstages)
    integrator.stats.nf2 += 1

    if integrator.opts.adaptive && !isempty(btilde)
        if pos_only_error
            uhat = utilde.x[2]
            @.. broadcast = false uhat = dtsq * btilde[1] * k1
            for i in 2:nstages
                ki = ks[i - 1]
                @.. broadcast = false uhat = uhat + dtsq * btilde[i] * ki
            end
            calculate_residuals!(
                atmp.x[2], uhat, integrator.uprev.x[2], integrator.u.x[2],
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp.x[2], t))
        else
            duhat, uhat = utilde.x
            @.. broadcast = false uhat = dtsq * btilde[1] * k1
            @.. broadcast = false duhat = dt * bptilde[1] * k1
            for i in 2:nstages
                ki = ks[i - 1]
                @.. broadcast = false uhat = uhat + dtsq * btilde[i] * ki
                @.. broadcast = false duhat = duhat + dt * bptilde[i] * ki
            end
            calculate_residuals!(
                atmp, utilde, integrator.uprev, integrator.u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
        end
    end
end

## Generic Nyström velocity-DEPENDENT perform_step!
## Solves: y'' = f(t, y, y') where f depends on both position and velocity
## kᵢ = f1(duprev + dt*Σⱼ abar[i,j]*kⱼ,  uprev + dt*c[i]*duprev + dt²*Σⱼ a[i,j]*kⱼ,  p, t+c[i]*dt)
## y₁ = y₀ + h*y'₀ + h²*Σᵢ bᵢ*kᵢ
## y'₁ = y'₀ + h*Σᵢ bpᵢ*kᵢ

function initialize!(integrator, cache::NystromVDConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    duprev, uprev = integrator.uprev.x
    kdu = integrator.f.f1(duprev, uprev, integrator.p, integrator.t)
    ku = integrator.f.f2(duprev, uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.fsalfirst = ArrayPartition((kdu, ku))
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::NystromVDCache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    duprev, uprev = integrator.uprev.x
    integrator.f.f1(integrator.fsalfirst.x[1], duprev, uprev, integrator.p, integrator.t)
    integrator.f.f2(integrator.fsalfirst.x[2], duprev, uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return integrator.stats.nf2 += 1
end

@muladd function perform_step!(
        integrator, cache::NystromVDConstantCache, repeat_step = false
    )
    (; t, dt, f, p) = integrator
    duprev, uprev = integrator.uprev.x
    (; tab) = cache
    (; a, abar, b, bp, btilde, bptilde, c) = tab
    k1 = integrator.fsalfirst.x[1]
    nstages = length(b)
    dtsq = dt^2
    ks = Vector{typeof(k1)}(undef, nstages)
    ks[1] = k1

    # ---------------- Stages 2..nstages (explicit ladder, max nstages = 17) ----------------
    if nstages >= 2
        ku = uprev + dt * (c[1] * duprev + dt * (a[2, 1] * ks[1]))
        kdu = duprev + dt * abar[2, 1] * ks[1]
        ks[2] = f.f1(kdu, ku, p, t + dt * c[1])
    end
    if nstages >= 3
        ku = uprev + dt * (c[2] * duprev + dt * (a[3, 1] * ks[1] + a[3, 2] * ks[2]))
        kdu = duprev + dt * (abar[3, 1] * ks[1] + abar[3, 2] * ks[2])
        ks[3] = f.f1(kdu, ku, p, t + dt * c[2])
    end
    if nstages >= 4
        ku = uprev + dt * (c[3] * duprev + dt * (a[4, 1] * ks[1] + a[4, 2] * ks[2] + a[4, 3] * ks[3]))
        kdu = duprev +
            dt * (abar[4, 1] * ks[1] + abar[4, 2] * ks[2] + abar[4, 3] * ks[3])
        ks[4] = f.f1(kdu, ku, p, t + dt * c[3])
    end
    if nstages >= 5
        ku = uprev + dt * (c[4] * duprev + dt * (a[5, 1] * ks[1] + a[5, 2] * ks[2] + a[5, 3] * ks[3] + a[5, 4] * ks[4]))
        kdu = duprev +
            dt * (
            abar[5, 1] * ks[1] + abar[5, 2] * ks[2] +
                abar[5, 3] * ks[3] + abar[5, 4] * ks[4]
        )
        ks[5] = f.f1(kdu, ku, p, t + dt * c[4])
    end
    if nstages >= 6
        ku = uprev + dt * (c[5] * duprev + dt * (a[6, 1] * ks[1] + a[6, 2] * ks[2] + a[6, 3] * ks[3] + a[6, 4] * ks[4] + a[6, 5] * ks[5]))
        kdu = duprev +
            dt * (
            abar[6, 1] * ks[1] + abar[6, 2] * ks[2] + abar[6, 3] * ks[3] +
                abar[6, 4] * ks[4] + abar[6, 5] * ks[5]
        )
        ks[6] = f.f1(kdu, ku, p, t + dt * c[5])
    end
    if nstages >= 7
        ku = uprev + dt * (c[6] * duprev + dt * (a[7, 1] * ks[1] + a[7, 2] * ks[2] + a[7, 3] * ks[3] + a[7, 4] * ks[4] + a[7, 5] * ks[5] + a[7, 6] * ks[6]))
        kdu = duprev +
            dt * (
            abar[7, 1] * ks[1] + abar[7, 2] * ks[2] + abar[7, 3] * ks[3] +
                abar[7, 4] * ks[4] + abar[7, 5] * ks[5] + abar[7, 6] * ks[6]
        )
        ks[7] = f.f1(kdu, ku, p, t + dt * c[6])
    end
    if nstages >= 8
        ku = uprev + dt * (c[7] * duprev + dt * (a[8, 1] * ks[1] + a[8, 2] * ks[2] + a[8, 3] * ks[3] + a[8, 4] * ks[4] + a[8, 5] * ks[5] + a[8, 6] * ks[6] + a[8, 7] * ks[7]))
        kdu = duprev +
            dt * (
            abar[8, 1] * ks[1] + abar[8, 2] * ks[2] + abar[8, 3] * ks[3] +
                abar[8, 4] * ks[4] + abar[8, 5] * ks[5] + abar[8, 6] * ks[6] +
                abar[8, 7] * ks[7]
        )
        ks[8] = f.f1(kdu, ku, p, t + dt * c[7])
    end
    if nstages >= 9
        ku = uprev + dt * (c[8] * duprev + dt * (a[9, 1] * ks[1] + a[9, 2] * ks[2] + a[9, 3] * ks[3] + a[9, 4] * ks[4] + a[9, 5] * ks[5] + a[9, 6] * ks[6] + a[9, 7] * ks[7] + a[9, 8] * ks[8]))
        kdu = duprev +
            dt * (
            abar[9, 1] * ks[1] + abar[9, 2] * ks[2] + abar[9, 3] * ks[3] +
                abar[9, 4] * ks[4] + abar[9, 5] * ks[5] + abar[9, 6] * ks[6] +
                abar[9, 7] * ks[7] + abar[9, 8] * ks[8]
        )
        ks[9] = f.f1(kdu, ku, p, t + dt * c[8])
    end
    if nstages >= 10
        ku = uprev + dt * (c[9] * duprev + dt * (a[10, 1] * ks[1] + a[10, 2] * ks[2] + a[10, 3] * ks[3] + a[10, 4] * ks[4] + a[10, 5] * ks[5] + a[10, 6] * ks[6] + a[10, 7] * ks[7] + a[10, 8] * ks[8] + a[10, 9] * ks[9]))
        kdu = duprev +
            dt * (
            abar[10, 1] * ks[1] + abar[10, 2] * ks[2] + abar[10, 3] * ks[3] +
                abar[10, 4] * ks[4] + abar[10, 5] * ks[5] + abar[10, 6] * ks[6] +
                abar[10, 7] * ks[7] + abar[10, 8] * ks[8] + abar[10, 9] * ks[9]
        )
        ks[10] = f.f1(kdu, ku, p, t + dt * c[9])
    end
    if nstages >= 11
        ku = uprev + dt * (c[10] * duprev + dt * (a[11, 1] * ks[1] + a[11, 2] * ks[2] + a[11, 3] * ks[3] + a[11, 4] * ks[4] + a[11, 5] * ks[5] + a[11, 6] * ks[6] + a[11, 7] * ks[7] + a[11, 8] * ks[8] + a[11, 9] * ks[9] + a[11, 10] * ks[10]))
        kdu = duprev +
            dt * (
            abar[11, 1] * ks[1] + abar[11, 2] * ks[2] + abar[11, 3] * ks[3] +
                abar[11, 4] * ks[4] + abar[11, 5] * ks[5] + abar[11, 6] * ks[6] +
                abar[11, 7] * ks[7] + abar[11, 8] * ks[8] + abar[11, 9] * ks[9] +
                abar[11, 10] * ks[10]
        )
        ks[11] = f.f1(kdu, ku, p, t + dt * c[10])
    end
    if nstages >= 12
        ku = uprev + dt * (c[11] * duprev + dt * (a[12, 1] * ks[1] + a[12, 2] * ks[2] + a[12, 3] * ks[3] + a[12, 4] * ks[4] + a[12, 5] * ks[5] + a[12, 6] * ks[6] + a[12, 7] * ks[7] + a[12, 8] * ks[8] + a[12, 9] * ks[9] + a[12, 10] * ks[10] + a[12, 11] * ks[11]))
        kdu = duprev +
            dt * (
            abar[12, 1] * ks[1] + abar[12, 2] * ks[2] + abar[12, 3] * ks[3] +
                abar[12, 4] * ks[4] + abar[12, 5] * ks[5] + abar[12, 6] * ks[6] +
                abar[12, 7] * ks[7] + abar[12, 8] * ks[8] + abar[12, 9] * ks[9] +
                abar[12, 10] * ks[10] + abar[12, 11] * ks[11]
        )
        ks[12] = f.f1(kdu, ku, p, t + dt * c[11])
    end
    if nstages >= 13
        ku = uprev + dt * (c[12] * duprev + dt * (a[13, 1] * ks[1] + a[13, 2] * ks[2] + a[13, 3] * ks[3] + a[13, 4] * ks[4] + a[13, 5] * ks[5] + a[13, 6] * ks[6] + a[13, 7] * ks[7] + a[13, 8] * ks[8] + a[13, 9] * ks[9] + a[13, 10] * ks[10] + a[13, 11] * ks[11] + a[13, 12] * ks[12]))
        kdu = duprev +
            dt * (
            abar[13, 1] * ks[1] + abar[13, 2] * ks[2] + abar[13, 3] * ks[3] +
                abar[13, 4] * ks[4] + abar[13, 5] * ks[5] + abar[13, 6] * ks[6] +
                abar[13, 7] * ks[7] + abar[13, 8] * ks[8] + abar[13, 9] * ks[9] +
                abar[13, 10] * ks[10] + abar[13, 11] * ks[11] + abar[13, 12] * ks[12]
        )
        ks[13] = f.f1(kdu, ku, p, t + dt * c[12])
    end
    if nstages >= 14
        ku = uprev + dt * (c[13] * duprev + dt * (a[14, 1] * ks[1] + a[14, 2] * ks[2] + a[14, 3] * ks[3] + a[14, 4] * ks[4] + a[14, 5] * ks[5] + a[14, 6] * ks[6] + a[14, 7] * ks[7] + a[14, 8] * ks[8] + a[14, 9] * ks[9] + a[14, 10] * ks[10] + a[14, 11] * ks[11] + a[14, 12] * ks[12] + a[14, 13] * ks[13]))
        kdu = duprev +
            dt * (
            abar[14, 1] * ks[1] + abar[14, 2] * ks[2] + abar[14, 3] * ks[3] +
                abar[14, 4] * ks[4] + abar[14, 5] * ks[5] + abar[14, 6] * ks[6] +
                abar[14, 7] * ks[7] + abar[14, 8] * ks[8] + abar[14, 9] * ks[9] +
                abar[14, 10] * ks[10] + abar[14, 11] * ks[11] + abar[14, 12] * ks[12] +
                abar[14, 13] * ks[13]
        )
        ks[14] = f.f1(kdu, ku, p, t + dt * c[13])
    end
    if nstages >= 15
        ku = uprev + dt * (c[14] * duprev + dt * (a[15, 1] * ks[1] + a[15, 2] * ks[2] + a[15, 3] * ks[3] + a[15, 4] * ks[4] + a[15, 5] * ks[5] + a[15, 6] * ks[6] + a[15, 7] * ks[7] + a[15, 8] * ks[8] + a[15, 9] * ks[9] + a[15, 10] * ks[10] + a[15, 11] * ks[11] + a[15, 12] * ks[12] + a[15, 13] * ks[13] + a[15, 14] * ks[14]))
        kdu = duprev +
            dt * (
            abar[15, 1] * ks[1] + abar[15, 2] * ks[2] + abar[15, 3] * ks[3] +
                abar[15, 4] * ks[4] + abar[15, 5] * ks[5] + abar[15, 6] * ks[6] +
                abar[15, 7] * ks[7] + abar[15, 8] * ks[8] + abar[15, 9] * ks[9] +
                abar[15, 10] * ks[10] + abar[15, 11] * ks[11] + abar[15, 12] * ks[12] +
                abar[15, 13] * ks[13] + abar[15, 14] * ks[14]
        )
        ks[15] = f.f1(kdu, ku, p, t + dt * c[14])
    end
    if nstages >= 16
        ku = uprev + dt * (c[15] * duprev + dt * (a[16, 1] * ks[1] + a[16, 2] * ks[2] + a[16, 3] * ks[3] + a[16, 4] * ks[4] + a[16, 5] * ks[5] + a[16, 6] * ks[6] + a[16, 7] * ks[7] + a[16, 8] * ks[8] + a[16, 9] * ks[9] + a[16, 10] * ks[10] + a[16, 11] * ks[11] + a[16, 12] * ks[12] + a[16, 13] * ks[13] + a[16, 14] * ks[14] + a[16, 15] * ks[15]))
        kdu = duprev +
            dt * (
            abar[16, 1] * ks[1] + abar[16, 2] * ks[2] + abar[16, 3] * ks[3] +
                abar[16, 4] * ks[4] + abar[16, 5] * ks[5] + abar[16, 6] * ks[6] +
                abar[16, 7] * ks[7] + abar[16, 8] * ks[8] + abar[16, 9] * ks[9] +
                abar[16, 10] * ks[10] + abar[16, 11] * ks[11] + abar[16, 12] * ks[12] +
                abar[16, 13] * ks[13] + abar[16, 14] * ks[14] + abar[16, 15] * ks[15]
        )
        ks[16] = f.f1(kdu, ku, p, t + dt * c[15])
    end
    if nstages >= 17
        ku = uprev + dt * (c[16] * duprev + dt * (a[17, 1] * ks[1] + a[17, 2] * ks[2] + a[17, 3] * ks[3] + a[17, 4] * ks[4] + a[17, 5] * ks[5] + a[17, 6] * ks[6] + a[17, 7] * ks[7] + a[17, 8] * ks[8] + a[17, 9] * ks[9] + a[17, 10] * ks[10] + a[17, 11] * ks[11] + a[17, 12] * ks[12] + a[17, 13] * ks[13] + a[17, 14] * ks[14] + a[17, 15] * ks[15] + a[17, 16] * ks[16]))
        kdu = duprev +
            dt * (
            abar[17, 1] * ks[1] + abar[17, 2] * ks[2] + abar[17, 3] * ks[3] +
                abar[17, 4] * ks[4] + abar[17, 5] * ks[5] + abar[17, 6] * ks[6] +
                abar[17, 7] * ks[7] + abar[17, 8] * ks[8] + abar[17, 9] * ks[9] +
                abar[17, 10] * ks[10] + abar[17, 11] * ks[11] + abar[17, 12] * ks[12] +
                abar[17, 13] * ks[13] + abar[17, 14] * ks[14] + abar[17, 15] * ks[15] +
                abar[17, 16] * ks[16]
        )
        ks[17] = f.f1(kdu, ku, p, t + dt * c[16])
    end

    # ---------------- Final update u, du ----------------
    u = uprev + dt * duprev
    du = duprev
    for i in 1:nstages
        u = u + dtsq * b[i] * ks[i]
        du = du + dt * bp[i] * ks[i]
    end

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, tab.nf_per_step)
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    if integrator.opts.adaptive && !isempty(btilde)
        uhat = zero(uprev)
        duhat = zero(duprev)
        for i in 1:nstages
            if !iszero(btilde[i])
                uhat = uhat + dtsq * btilde[i] * ks[i]
            end
            if !isempty(bptilde) && !iszero(bptilde[i])
                duhat = duhat + dt * bptilde[i] * ks[i]
            end
        end
        utilde = ArrayPartition((duhat, uhat))
        atmp = calculate_residuals(
            utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

@muladd function perform_step!(
        integrator, cache::NystromVDCache, repeat_step = false
    )
    (; t, dt, f, p) = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    (; ks, k, utilde, tmp, atmp, tab) = cache
    (; a, abar, b, bp, btilde, bptilde, c) = tab
    ku = tmp.x[2]
    kdu = tmp.x[1]
    k1 = integrator.fsalfirst.x[1]
    nstages = length(b)
    dtsq = dt^2

    # ---------------- Stages 2..nstages (explicit ladder; ks[i-1] holds the i-th stage) ----------------
    if nstages >= 2
        @.. broadcast = false ku = uprev + dt * (c[1] * duprev + dt * (a[2, 1] * k1))
        @.. broadcast = false kdu = duprev + dt * abar[2, 1] * k1
        f.f1(ks[1], kdu, ku, p, t + dt * c[1])
    end
    if nstages >= 3
        @.. broadcast = false ku = uprev + dt * (c[2] * duprev + dt * (a[3, 1] * k1 + a[3, 2] * ks[1]))
        @.. broadcast = false kdu = duprev + dt * abar[3, 1] * k1 + dt * abar[3, 2] * ks[1]
        f.f1(ks[2], kdu, ku, p, t + dt * c[2])
    end
    if nstages >= 4
        @.. broadcast = false ku = uprev + dt * (c[3] * duprev + dt * (a[4, 1] * k1 + a[4, 2] * ks[1] + a[4, 3] * ks[2]))
        @.. broadcast = false kdu = duprev +
            dt * abar[4, 1] * k1 + dt * abar[4, 2] * ks[1] + dt * abar[4, 3] * ks[2]
        f.f1(ks[3], kdu, ku, p, t + dt * c[3])
    end
    if nstages >= 5
        @.. broadcast = false ku = uprev + dt * (c[4] * duprev + dt * (a[5, 1] * k1 + a[5, 2] * ks[1] + a[5, 3] * ks[2] + a[5, 4] * ks[3]))
        @.. broadcast = false kdu = duprev +
            dt * abar[5, 1] * k1 + dt * abar[5, 2] * ks[1] +
            dt * abar[5, 3] * ks[2] + dt * abar[5, 4] * ks[3]
        f.f1(ks[4], kdu, ku, p, t + dt * c[4])
    end
    if nstages >= 6
        @.. broadcast = false ku = uprev + dt * (c[5] * duprev + dt * (a[6, 1] * k1 + a[6, 2] * ks[1] + a[6, 3] * ks[2] + a[6, 4] * ks[3] + a[6, 5] * ks[4]))
        @.. broadcast = false kdu = duprev +
            dt * abar[6, 1] * k1 + dt * abar[6, 2] * ks[1] +
            dt * abar[6, 3] * ks[2] + dt * abar[6, 4] * ks[3] +
            dt * abar[6, 5] * ks[4]
        f.f1(ks[5], kdu, ku, p, t + dt * c[5])
    end
    if nstages >= 7
        @.. broadcast = false ku = uprev + dt * (c[6] * duprev + dt * (a[7, 1] * k1 + a[7, 2] * ks[1] + a[7, 3] * ks[2] + a[7, 4] * ks[3] + a[7, 5] * ks[4] + a[7, 6] * ks[5]))
        @.. broadcast = false kdu = duprev +
            dt * abar[7, 1] * k1 + dt * abar[7, 2] * ks[1] +
            dt * abar[7, 3] * ks[2] + dt * abar[7, 4] * ks[3] +
            dt * abar[7, 5] * ks[4] + dt * abar[7, 6] * ks[5]
        f.f1(ks[6], kdu, ku, p, t + dt * c[6])
    end
    if nstages >= 8
        @.. broadcast = false ku = uprev + dt * (c[7] * duprev + dt * (a[8, 1] * k1 + a[8, 2] * ks[1] + a[8, 3] * ks[2] + a[8, 4] * ks[3] + a[8, 5] * ks[4] + a[8, 6] * ks[5] + a[8, 7] * ks[6]))
        @.. broadcast = false kdu = duprev +
            dt * abar[8, 1] * k1 + dt * abar[8, 2] * ks[1] +
            dt * abar[8, 3] * ks[2] + dt * abar[8, 4] * ks[3] +
            dt * abar[8, 5] * ks[4] + dt * abar[8, 6] * ks[5] +
            dt * abar[8, 7] * ks[6]
        f.f1(ks[7], kdu, ku, p, t + dt * c[7])
    end
    if nstages >= 9
        @.. broadcast = false ku = uprev + dt * (c[8] * duprev + dt * (a[9, 1] * k1 + a[9, 2] * ks[1] + a[9, 3] * ks[2] + a[9, 4] * ks[3] + a[9, 5] * ks[4] + a[9, 6] * ks[5] + a[9, 7] * ks[6] + a[9, 8] * ks[7]))
        @.. broadcast = false kdu = duprev +
            dt * abar[9, 1] * k1 + dt * abar[9, 2] * ks[1] +
            dt * abar[9, 3] * ks[2] + dt * abar[9, 4] * ks[3] +
            dt * abar[9, 5] * ks[4] + dt * abar[9, 6] * ks[5] +
            dt * abar[9, 7] * ks[6] + dt * abar[9, 8] * ks[7]
        f.f1(ks[8], kdu, ku, p, t + dt * c[8])
    end
    if nstages >= 10
        @.. broadcast = false ku = uprev + dt * (c[9] * duprev + dt * (a[10, 1] * k1 + a[10, 2] * ks[1] + a[10, 3] * ks[2] + a[10, 4] * ks[3] + a[10, 5] * ks[4] + a[10, 6] * ks[5] + a[10, 7] * ks[6] + a[10, 8] * ks[7] + a[10, 9] * ks[8]))
        @.. broadcast = false kdu = duprev +
            dt * abar[10, 1] * k1 + dt * abar[10, 2] * ks[1] +
            dt * abar[10, 3] * ks[2] + dt * abar[10, 4] * ks[3] +
            dt * abar[10, 5] * ks[4] + dt * abar[10, 6] * ks[5] +
            dt * abar[10, 7] * ks[6] + dt * abar[10, 8] * ks[7] +
            dt * abar[10, 9] * ks[8]
        f.f1(ks[9], kdu, ku, p, t + dt * c[9])
    end
    if nstages >= 11
        @.. broadcast = false ku = uprev + dt * (c[10] * duprev + dt * (a[11, 1] * k1 + a[11, 2] * ks[1] + a[11, 3] * ks[2] + a[11, 4] * ks[3] + a[11, 5] * ks[4] + a[11, 6] * ks[5] + a[11, 7] * ks[6] + a[11, 8] * ks[7] + a[11, 9] * ks[8] + a[11, 10] * ks[9]))
        @.. broadcast = false kdu = duprev +
            dt * abar[11, 1] * k1 + dt * abar[11, 2] * ks[1] +
            dt * abar[11, 3] * ks[2] + dt * abar[11, 4] * ks[3] +
            dt * abar[11, 5] * ks[4] + dt * abar[11, 6] * ks[5] +
            dt * abar[11, 7] * ks[6] + dt * abar[11, 8] * ks[7] +
            dt * abar[11, 9] * ks[8] + dt * abar[11, 10] * ks[9]
        f.f1(ks[10], kdu, ku, p, t + dt * c[10])
    end
    if nstages >= 12
        @.. broadcast = false ku = uprev + dt * (c[11] * duprev + dt * (a[12, 1] * k1 + a[12, 2] * ks[1] + a[12, 3] * ks[2] + a[12, 4] * ks[3] + a[12, 5] * ks[4] + a[12, 6] * ks[5] + a[12, 7] * ks[6] + a[12, 8] * ks[7] + a[12, 9] * ks[8] + a[12, 10] * ks[9] + a[12, 11] * ks[10]))
        @.. broadcast = false kdu = duprev +
            dt * abar[12, 1] * k1 + dt * abar[12, 2] * ks[1] +
            dt * abar[12, 3] * ks[2] + dt * abar[12, 4] * ks[3] +
            dt * abar[12, 5] * ks[4] + dt * abar[12, 6] * ks[5] +
            dt * abar[12, 7] * ks[6] + dt * abar[12, 8] * ks[7] +
            dt * abar[12, 9] * ks[8] + dt * abar[12, 10] * ks[9] +
            dt * abar[12, 11] * ks[10]
        f.f1(ks[11], kdu, ku, p, t + dt * c[11])
    end
    if nstages >= 13
        @.. broadcast = false ku = uprev + dt * (c[12] * duprev + dt * (a[13, 1] * k1 + a[13, 2] * ks[1] + a[13, 3] * ks[2] + a[13, 4] * ks[3] + a[13, 5] * ks[4] + a[13, 6] * ks[5] + a[13, 7] * ks[6] + a[13, 8] * ks[7] + a[13, 9] * ks[8] + a[13, 10] * ks[9] + a[13, 11] * ks[10] + a[13, 12] * ks[11]))
        @.. broadcast = false kdu = duprev +
            dt * abar[13, 1] * k1 + dt * abar[13, 2] * ks[1] +
            dt * abar[13, 3] * ks[2] + dt * abar[13, 4] * ks[3] +
            dt * abar[13, 5] * ks[4] + dt * abar[13, 6] * ks[5] +
            dt * abar[13, 7] * ks[6] + dt * abar[13, 8] * ks[7] +
            dt * abar[13, 9] * ks[8] + dt * abar[13, 10] * ks[9] +
            dt * abar[13, 11] * ks[10] + dt * abar[13, 12] * ks[11]
        f.f1(ks[12], kdu, ku, p, t + dt * c[12])
    end
    if nstages >= 14
        @.. broadcast = false ku = uprev + dt * (c[13] * duprev + dt * (a[14, 1] * k1 + a[14, 2] * ks[1] + a[14, 3] * ks[2] + a[14, 4] * ks[3] + a[14, 5] * ks[4] + a[14, 6] * ks[5] + a[14, 7] * ks[6] + a[14, 8] * ks[7] + a[14, 9] * ks[8] + a[14, 10] * ks[9] + a[14, 11] * ks[10] + a[14, 12] * ks[11] + a[14, 13] * ks[12]))
        @.. broadcast = false kdu = duprev +
            dt * abar[14, 1] * k1 + dt * abar[14, 2] * ks[1] +
            dt * abar[14, 3] * ks[2] + dt * abar[14, 4] * ks[3] +
            dt * abar[14, 5] * ks[4] + dt * abar[14, 6] * ks[5] +
            dt * abar[14, 7] * ks[6] + dt * abar[14, 8] * ks[7] +
            dt * abar[14, 9] * ks[8] + dt * abar[14, 10] * ks[9] +
            dt * abar[14, 11] * ks[10] + dt * abar[14, 12] * ks[11] +
            dt * abar[14, 13] * ks[12]
        f.f1(ks[13], kdu, ku, p, t + dt * c[13])
    end
    if nstages >= 15
        @.. broadcast = false ku = uprev + dt * (c[14] * duprev + dt * (a[15, 1] * k1 + a[15, 2] * ks[1] + a[15, 3] * ks[2] + a[15, 4] * ks[3] + a[15, 5] * ks[4] + a[15, 6] * ks[5] + a[15, 7] * ks[6] + a[15, 8] * ks[7] + a[15, 9] * ks[8] + a[15, 10] * ks[9] + a[15, 11] * ks[10] + a[15, 12] * ks[11] + a[15, 13] * ks[12] + a[15, 14] * ks[13]))
        @.. broadcast = false kdu = duprev +
            dt * abar[15, 1] * k1 + dt * abar[15, 2] * ks[1] +
            dt * abar[15, 3] * ks[2] + dt * abar[15, 4] * ks[3] +
            dt * abar[15, 5] * ks[4] + dt * abar[15, 6] * ks[5] +
            dt * abar[15, 7] * ks[6] + dt * abar[15, 8] * ks[7] +
            dt * abar[15, 9] * ks[8] + dt * abar[15, 10] * ks[9] +
            dt * abar[15, 11] * ks[10] + dt * abar[15, 12] * ks[11] +
            dt * abar[15, 13] * ks[12] + dt * abar[15, 14] * ks[13]
        f.f1(ks[14], kdu, ku, p, t + dt * c[14])
    end
    if nstages >= 16
        @.. broadcast = false ku = uprev + dt * (c[15] * duprev + dt * (a[16, 1] * k1 + a[16, 2] * ks[1] + a[16, 3] * ks[2] + a[16, 4] * ks[3] + a[16, 5] * ks[4] + a[16, 6] * ks[5] + a[16, 7] * ks[6] + a[16, 8] * ks[7] + a[16, 9] * ks[8] + a[16, 10] * ks[9] + a[16, 11] * ks[10] + a[16, 12] * ks[11] + a[16, 13] * ks[12] + a[16, 14] * ks[13] + a[16, 15] * ks[14]))
        @.. broadcast = false kdu = duprev +
            dt * abar[16, 1] * k1 + dt * abar[16, 2] * ks[1] +
            dt * abar[16, 3] * ks[2] + dt * abar[16, 4] * ks[3] +
            dt * abar[16, 5] * ks[4] + dt * abar[16, 6] * ks[5] +
            dt * abar[16, 7] * ks[6] + dt * abar[16, 8] * ks[7] +
            dt * abar[16, 9] * ks[8] + dt * abar[16, 10] * ks[9] +
            dt * abar[16, 11] * ks[10] + dt * abar[16, 12] * ks[11] +
            dt * abar[16, 13] * ks[12] + dt * abar[16, 14] * ks[13] +
            dt * abar[16, 15] * ks[14]
        f.f1(ks[15], kdu, ku, p, t + dt * c[15])
    end
    if nstages >= 17
        @.. broadcast = false ku = uprev + dt * (c[16] * duprev + dt * (a[17, 1] * k1 + a[17, 2] * ks[1] + a[17, 3] * ks[2] + a[17, 4] * ks[3] + a[17, 5] * ks[4] + a[17, 6] * ks[5] + a[17, 7] * ks[6] + a[17, 8] * ks[7] + a[17, 9] * ks[8] + a[17, 10] * ks[9] + a[17, 11] * ks[10] + a[17, 12] * ks[11] + a[17, 13] * ks[12] + a[17, 14] * ks[13] + a[17, 15] * ks[14] + a[17, 16] * ks[15]))
        @.. broadcast = false kdu = duprev +
            dt * abar[17, 1] * k1 + dt * abar[17, 2] * ks[1] +
            dt * abar[17, 3] * ks[2] + dt * abar[17, 4] * ks[3] +
            dt * abar[17, 5] * ks[4] + dt * abar[17, 6] * ks[5] +
            dt * abar[17, 7] * ks[6] + dt * abar[17, 8] * ks[7] +
            dt * abar[17, 9] * ks[8] + dt * abar[17, 10] * ks[9] +
            dt * abar[17, 11] * ks[10] + dt * abar[17, 12] * ks[11] +
            dt * abar[17, 13] * ks[12] + dt * abar[17, 14] * ks[13] +
            dt * abar[17, 15] * ks[14] + dt * abar[17, 16] * ks[15]
        f.f1(ks[16], kdu, ku, p, t + dt * c[16])
    end

    # ---------------- Final update u, du ----------------
    @.. broadcast = false u = uprev + dt * duprev
    @.. broadcast = false du = duprev
    @.. broadcast = false u = u + dtsq * b[1] * k1
    @.. broadcast = false du = du + dt * bp[1] * k1
    for i in 2:nstages
        ki = ks[i - 1]
        @.. broadcast = false u = u + dtsq * b[i] * ki
        @.. broadcast = false du = du + dt * bp[i] * ki
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, tab.nf_per_step)
    integrator.stats.nf2 += 1

    if integrator.opts.adaptive && !isempty(btilde)
        duhat, uhat = utilde.x
        @.. broadcast = false uhat = zero(uhat)
        @.. broadcast = false duhat = zero(duhat)
        for i in 1:nstages
            ki = (i == 1) ? k1 : ks[i - 1]
            if !iszero(btilde[i])
                @.. broadcast = false uhat = uhat + dtsq * btilde[i] * ki
            end
            if !isempty(bptilde) && !iszero(bptilde[i])
                @.. broadcast = false duhat = duhat + dt * bptilde[i] * ki
            end
        end
        calculate_residuals!(
            atmp, utilde, integrator.uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end
