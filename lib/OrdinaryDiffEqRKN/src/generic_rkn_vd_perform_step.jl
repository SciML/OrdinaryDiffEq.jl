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

    # Compute intermediate stages k2..knstages
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

    # Position and velocity updates
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
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
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

    # Compute intermediate stages k2..knstages, stored in ks[1..nstages-1]
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

    # Position update: u = uprev + dt*duprev + dt^2 * sum(b[i]*ki)
    @.. broadcast = false u = uprev + dt * duprev
    for i in 1:nstages
        if !iszero(b[i])
            ki = (i == 1) ? k1 : ks[i - 1]
            @.. broadcast = false u = u + dtsq * b[i] * ki
        end
    end

    # Velocity update: du = duprev + dt * sum(bp[i]*ki)
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
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end
