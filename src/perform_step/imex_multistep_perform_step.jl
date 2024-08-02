# CNAB2

function initialize!(integrator, cache::CNAB2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f.f1(integrator.uprev, integrator.p, integrator.t) +
                           integrator.f.f2(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::CNAB2ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack k2, nlsolver = cache
    cnt = integrator.iter
    f1 = integrator.f.f1
    f2 = integrator.f.f2
    du₁ = f1(uprev, p, t)
    integrator.stats.nf += 1
    k1 = integrator.fsalfirst - du₁
    # Explicit part
    if cnt == 1
        tmp = uprev + dt * k1
    else
        tmp = uprev + dt * (3 // 2 * k1 - 1 // 2 * k2)
    end
    nlsolver.tmp = tmp
    # Implicit part
    # precalculations
    γ = 1 // 2
    γdt = γ * dt

    # initial guess
    zprev = dt * du₁
    nlsolver.z = z = zprev # Constant extrapolation

    nlsolver.tmp += γ * zprev
    markfirststage!(nlsolver)
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = nlsolver.tmp + 1 // 2 * z

    cache.k2 = k1
    integrator.fsallast = f1(u, p, t + dt) + f2(u, p, t + dt)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::CNAB2Cache)
    integrator.kshortsize = 2
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
end

function perform_step!(integrator, cache::CNAB2Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack k1, k2, du₁, nlsolver = cache
    @unpack z, tmp = nlsolver
    @unpack f1 = f
    cnt = integrator.iter

    f1(du₁, uprev, p, t)
    integrator.stats.nf += 1
    @.. broadcast=false k1=integrator.fsalfirst - du₁
    # Explicit part
    if cnt == 1
        @.. broadcast=false tmp=uprev + dt * k1
    else
        @.. broadcast=false tmp=uprev + dt * (3 // 2 * k1 - 1 // 2 * k2)
    end
    # Implicit part
    # precalculations
    γ = 1 // 2
    γdt = γ * dt

    # initial guess
    @.. broadcast=false z=dt * du₁
    @.. broadcast=false tmp+=γ * z
    markfirststage!(nlsolver)
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast=false u=tmp + 1 // 2 * z

    cache.k2 .= k1
    f(integrator.fsallast, u, p, t + dt)
    integrator.stats.nf += 1
end

# CNLF2

function initialize!(integrator, cache::CNLF2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f.f1(integrator.uprev, integrator.p, integrator.t) +
                           integrator.f.f2(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::CNLF2ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack k2, uprev2, nlsolver = cache
    cnt = integrator.iter
    f1 = integrator.f.f1
    f2 = integrator.f.f2
    du₁ = f1(uprev, p, t)
    integrator.stats.nf += 1
    # Explicit part
    if cnt == 1
        tmp = uprev + dt * (integrator.fsalfirst - du₁)
    else
        tmp = uprev2 + 2 // 1 * dt * (integrator.fsalfirst - du₁)
    end
    # Implicit part
    # precalculations
    γ = 1 // 1
    if cnt != 1
        tmp += γ * dt * k2
    end
    γdt = γ * dt
    nlsolver.tmp = tmp

    # initial guess
    zprev = dt * du₁
    nlsolver.z = z = zprev # Constant extrapolation

    markfirststage!(nlsolver)
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = nlsolver.tmp + γ * z

    cache.uprev2 = uprev
    cache.k2 = du₁
    integrator.fsallast = f1(u, p, t + dt) + f2(u, p, t + dt)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::CNLF2Cache)
    integrator.kshortsize = 2
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
end

function perform_step!(integrator, cache::CNLF2Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack uprev2, k2, du₁, nlsolver = cache
    @unpack z, tmp = nlsolver
    @unpack f1 = f
    cnt = integrator.iter

    f1(du₁, uprev, p, t)
    integrator.stats.nf += 1
    # Explicit part
    if cnt == 1
        @.. broadcast=false tmp=uprev + dt * (integrator.fsalfirst - du₁)
    else
        @.. broadcast=false tmp=uprev2 + 2 // 1 * dt * (integrator.fsalfirst - du₁)
    end
    # Implicit part
    # precalculations
    γ = 1 // 1
    if cnt != 1
        @.. broadcast=false tmp+=γ * dt * k2
    end
    γdt = γ * dt

    # initial guess
    @.. broadcast=false z=dt * du₁
    markfirststage!(nlsolver)
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast=false u=tmp + γ * z

    cache.uprev2 .= uprev
    cache.k2 .= du₁
    f(integrator.fsallast, u, p, t + dt)
    integrator.stats.nf += 1
end