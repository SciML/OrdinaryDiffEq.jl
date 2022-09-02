function initialize!(integrator, cache::ExplicitRKConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.destats.nf += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ExplicitRKConstantCache,
                               repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack A, c, α, αEEst, stages = cache
    @unpack kk = cache

    # Calc First
    kk[1] = integrator.fsalfirst

    # Calc Middle
    for i in 2:(stages - 1)
        utilde = zero(kk[1])
        for j in 1:(i - 1)
            utilde = utilde + A[j, i] * kk[j]
        end
        kk[i] = f(uprev + dt * utilde, p, t + c[i] * dt)
        integrator.destats.nf += 1
    end

    #Calc Last
    utilde = zero(kk[1])
    for j in 1:(stages - 1)
        utilde = utilde + A[j, end] * kk[j]
    end
    kk[end] = f(uprev + dt * utilde, p, t + c[end] * dt)
    integrator.destats.nf += 1
    integrator.fsallast = kk[end] # Uses fsallast as temp even if not fsal

    # Accumulate Result
    utilde = α[1] * kk[1]
    for i in 2:stages
        utilde = utilde + α[i] * kk[i]
    end
    u = uprev + dt * utilde

    if integrator.opts.adaptive
        utilde = (α[1] - αEEst[1]) .* kk[1]
        for i in 2:stages
            utilde = utilde + (α[i] - αEEst[i]) * kk[i]
        end
        atmp = calculate_residuals(dt * utilde, uprev, u, integrator.opts.abstol,
                                   integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if !isfsal(alg.tableau)
        integrator.fsallast = f(u, p, t + dt)
        integrator.destats.nf += 1
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::ExplicitRKCache)
    integrator.kshortsize = 2
    integrator.fsallast = cache.fsallast
    integrator.fsalfirst = cache.kk[1]
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::ExplicitRKCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    alg = unwrap_alg(integrator, nothing)
    @unpack A, c, α, αEEst, stages = cache.tab
    @unpack kk, utilde, tmp, atmp = cache

    # Middle
    for i in 2:(stages - 1)
        @.. broadcast=false utilde=zero(kk[1][1])
        for j in 1:(i - 1)
            @.. broadcast=false utilde=utilde + A[j, i] * kk[j]
        end
        @.. broadcast=false tmp=uprev + dt * utilde
        f(kk[i], tmp, p, t + c[i] * dt)
        integrator.destats.nf += 1
    end

    #Last
    @.. broadcast=false utilde=zero(kk[1][1])
    for j in 1:(stages - 1)
        @.. broadcast=false utilde=utilde + A[j, end] * kk[j]
    end
    @.. broadcast=false u=uprev + dt * utilde
    f(kk[end], u, p, t + c[end] * dt) #fsallast is tmp even if not fsal
    integrator.destats.nf += 1

    #Accumulate
    if !isfsal(alg.tableau)
        @.. broadcast=false utilde=α[1] * kk[1]
        for i in 2:stages
            @.. broadcast=false utilde=utilde + α[i] * kk[i]
        end
        @.. broadcast=false u=uprev + dt * utilde
    end

    if integrator.opts.adaptive
        @.. broadcast=false utilde=(α[1] - αEEst[1]) * kk[1]
        for i in 2:stages
            @.. broadcast=false utilde=utilde + (α[i] - αEEst[i]) * kk[i]
        end
        @.. broadcast=false tmp=dt * utilde
        calculate_residuals!(atmp, tmp, uprev, u,
                             integrator.opts.abstol, integrator.opts.reltol,
                             integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if !isfsal(alg.tableau)
        f(integrator.fsallast, u, p, t + dt)
        integrator.destats.nf += 1
    end
end
