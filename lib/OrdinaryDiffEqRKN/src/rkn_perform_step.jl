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
# its specialized "free" 6th-order interpolant (see interpolants.jl).
_rkn_kshortsize(::Any) = 2
_rkn_kshortsize(::DPRKN6Tableau) = 3

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
function _init_vi_k_mut!(integrator, cache, ::Any)
    integrator.kshortsize = 2
    resize!(integrator.k, 2)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return
end
function _init_vi_k_mut!(integrator, cache, ::DPRKN6Tableau)
    integrator.kshortsize = 3
    resize!(integrator.k, 3)
    integrator.k[1] = ArrayPartition(integrator.fsalfirst.x[1], cache.ks[1])
    integrator.k[2] = ArrayPartition(cache.ks[2], cache.ks[3])
    integrator.k[3] = ArrayPartition(cache.ks[4], cache.ks[5])
    return
end

function initialize!(integrator, cache::NystromVIConstantCache)
    integrator.kshortsize = _rkn_kshortsize(cache.tab)
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
    for i in 2:nstages
        ku = uprev + dt * c[i - 1] * duprev
        for j in 1:(i - 1)
            if !iszero(a[i, j])
                ku = ku + dtsq * a[i, j] * ks[j]
            end
        end
        ks[i] = f.f1(duprev, ku, p, t + dt * c[i - 1])
    end

    u = uprev + dt * duprev
    for i in 1:nstages
        if !iszero(b[i])
            u = u + dtsq * b[i] * ks[i]
        end
    end
    du = duprev
    for i in 1:nstages
        if !iszero(bp[i])
            du = du + dt * bp[i] * ks[i]
        end
    end

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((f.f1(du, u, p, t + dt), f.f2(du, u, p, t + dt)))
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, nstages)
    integrator.stats.nf2 += 1
    _store_vi_k!(integrator, tab, ks)

    if integrator.opts.adaptive && !isempty(btilde)
        uhat = zero(uprev)
        for i in 1:nstages
            if !iszero(btilde[i])
                uhat = uhat + dtsq * btilde[i] * ks[i]
            end
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

    for i in 2:nstages
        @.. broadcast = false ku = uprev + dt * c[i - 1] * duprev
        for j in 1:(i - 1)
            if !iszero(a[i, j])
                kj = (j == 1) ? k1 : ks[j - 1]
                @.. broadcast = false ku = ku + dtsq * a[i, j] * kj
            end
        end
        f.f1(ks[i - 1], duprev, ku, p, t + dt * c[i - 1])
    end

    @.. broadcast = false u = uprev + dt * duprev
    for i in 1:nstages
        if !iszero(b[i])
            ki = (i == 1) ? k1 : ks[i - 1]
            @.. broadcast = false u = u + dtsq * b[i] * ki
        end
    end

    @.. broadcast = false du = duprev
    for i in 1:nstages
        if !iszero(bp[i])
            ki = (i == 1) ? k1 : ks[i - 1]
            @.. broadcast = false du = du + dt * bp[i] * ki
        end
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, nstages)
    integrator.stats.nf2 += 1

    if integrator.opts.adaptive && !isempty(btilde)
        if pos_only_error
            uhat = utilde.x[2]
            @.. broadcast = false uhat = zero(uhat)
            for i in 1:nstages
                if !iszero(btilde[i])
                    ki = (i == 1) ? k1 : ks[i - 1]
                    @.. broadcast = false uhat = uhat + dtsq * btilde[i] * ki
                end
            end
            calculate_residuals!(
                atmp.x[2], uhat, integrator.uprev.x[2], integrator.u.x[2],
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp.x[2], t))
        else
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
    for i in 2:nstages
        ku = uprev + dt * c[i - 1] * duprev
        kdu = duprev
        for j in 1:(i - 1)
            if !iszero(a[i, j])
                ku = ku + dtsq * a[i, j] * ks[j]
            end
            if !iszero(abar[i, j])
                kdu = kdu + dt * abar[i, j] * ks[j]
            end
        end
        ks[i] = f.f1(kdu, ku, p, t + dt * c[i - 1])
    end

    u = uprev + dt * duprev
    for i in 1:nstages
        if !iszero(b[i])
            u = u + dtsq * b[i] * ks[i]
        end
    end
    du = duprev
    for i in 1:nstages
        if !iszero(bp[i])
            du = du + dt * bp[i] * ks[i]
        end
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

    for i in 2:nstages
        @.. broadcast = false ku = uprev + dt * c[i - 1] * duprev
        @.. broadcast = false kdu = duprev
        for j in 1:(i - 1)
            kj = (j == 1) ? k1 : ks[j - 1]
            if !iszero(a[i, j])
                @.. broadcast = false ku = ku + dtsq * a[i, j] * kj
            end
            if !iszero(abar[i, j])
                @.. broadcast = false kdu = kdu + dt * abar[i, j] * kj
            end
        end
        f.f1(ks[i - 1], kdu, ku, p, t + dt * c[i - 1])
    end

    @.. broadcast = false u = uprev + dt * duprev
    for i in 1:nstages
        if !iszero(b[i])
            ki = (i == 1) ? k1 : ks[i - 1]
            @.. broadcast = false u = u + dtsq * b[i] * ki
        end
    end

    @.. broadcast = false du = duprev
    for i in 1:nstages
        if !iszero(bp[i])
            ki = (i == 1) ? k1 : ks[i - 1]
            @.. broadcast = false du = du + dt * bp[i] * ki
        end
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
