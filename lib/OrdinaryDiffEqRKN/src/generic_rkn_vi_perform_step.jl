## Generic Nyström velocity-independent perform_step!
## Solves: y'' = f(t, y) where f is velocity-independent
## kᵢ = f1(duprev, yᵢ, p, t + cᵢ*dt)   (duprev constant throughout)
## yᵢ = y₀ + cᵢ*h*y'₀ + h²*Σⱼ<ᵢ aᵢⱼ*kⱼ
## y₁ = y₀ + h*y'₀ + h²*Σᵢ bᵢ*kᵢ
## y'₁ = y'₀ + h*Σᵢ bpᵢ*kᵢ

function initialize!(integrator, cache::NystromVIConstantCache)
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

function initialize!(integrator, cache::NystromVICache)
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
        integrator, cache::NystromVIConstantCache, repeat_step = false)
    (; t, dt, f, p) = integrator
    duprev, uprev = integrator.uprev.x
    (; tab) = cache
    (; a, b, bp, btilde, bptilde, c, pos_only_error) = tab
    k1 = integrator.fsalfirst.x[1]
    nstages = length(b)
    dtsq = dt^2

    # Compute intermediate stages
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
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, nstages)
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    if integrator.opts.adaptive && !isempty(btilde)
        uhat = zero(uprev)
        for i in 1:nstages
            if !iszero(btilde[i])
                uhat = uhat + dtsq * btilde[i] * ks[i]
            end
        end
        if pos_only_error
            atmp = calculate_residuals(uhat, integrator.uprev.x[2], integrator.u.x[2],
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        else
            duhat = zero(duprev)
            for i in 1:nstages
                if !isempty(bptilde) && !iszero(bptilde[i])
                    duhat = duhat + dt * bptilde[i] * ks[i]
                end
            end
            utilde = ArrayPartition((duhat, uhat))
            atmp = calculate_residuals(utilde, integrator.uprev, integrator.u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        end
    end
end

@muladd function perform_step!(
        integrator, cache::NystromVICache, repeat_step = false)
    (; t, dt, f, p) = integrator
    du, u = integrator.u.x
    duprev, uprev = integrator.uprev.x
    (; ks, k, utilde, tmp, atmp, tab) = cache
    (; a, b, bp, btilde, bptilde, c, pos_only_error) = tab
    ku = tmp.x[2]
    k1 = integrator.fsalfirst.x[1]
    nstages = length(b)
    dtsq = dt^2

    # Compute intermediate stages k2..knstages, stored in ks[1..nstages-1]
    for i in 2:nstages
        @.. broadcast=false ku = uprev + dt * c[i - 1] * duprev
        for j in 1:(i - 1)
            if !iszero(a[i, j])
                kj = (j == 1) ? k1 : ks[j - 1]
                @.. broadcast=false ku = ku + dtsq * a[i, j] * kj
            end
        end
        f.f1(ks[i - 1], duprev, ku, p, t + dt * c[i - 1])
    end

    # Position update: u = uprev + dt*duprev + dt^2 * sum(b[i]*ki)
    @.. broadcast=false u = uprev + dt * duprev
    for i in 1:nstages
        if !iszero(b[i])
            ki = (i == 1) ? k1 : ks[i - 1]
            @.. broadcast=false u = u + dtsq * b[i] * ki
        end
    end

    # Velocity update: du = duprev + dt * sum(bp[i]*ki)
    @.. broadcast=false du = duprev
    for i in 1:nstages
        if !iszero(bp[i])
            ki = (i == 1) ? k1 : ks[i - 1]
            @.. broadcast=false du = du + dt * bp[i] * ki
        end
    end

    f.f1(k.x[1], du, u, p, t + dt)
    f.f2(k.x[2], du, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, nstages)
    integrator.stats.nf2 += 1

    if integrator.opts.adaptive && !isempty(btilde)
        if pos_only_error
            uhat = utilde.x[2]
            @.. broadcast=false uhat = zero(uhat)
            for i in 1:nstages
                if !iszero(btilde[i])
                    ki = (i == 1) ? k1 : ks[i - 1]
                    @.. broadcast=false uhat = uhat + dtsq * btilde[i] * ki
                end
            end
            calculate_residuals!(atmp.x[2], uhat, integrator.uprev.x[2], integrator.u.x[2],
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(atmp.x[2], t)
        else
            duhat, uhat = utilde.x
            @.. broadcast=false uhat = zero(uhat)
            @.. broadcast=false duhat = zero(duhat)
            for i in 1:nstages
                ki = (i == 1) ? k1 : ks[i - 1]
                if !iszero(btilde[i])
                    @.. broadcast=false uhat = uhat + dtsq * btilde[i] * ki
                end
                if !isempty(bptilde) && !iszero(bptilde[i])
                    @.. broadcast=false duhat = duhat + dt * bptilde[i] * ki
                end
            end
            calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        end
    end
end
