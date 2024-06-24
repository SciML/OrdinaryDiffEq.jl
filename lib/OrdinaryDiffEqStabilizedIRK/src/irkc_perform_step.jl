function initialize!(integrator, cache::IRKCConstantCache)
    @unpack uprev, p, t = integrator
    @unpack f1, f2 = integrator.f
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    cache.du₁ = f1(uprev, p, t)
    cache.du₂ = f2(uprev, p, t)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
    integrator.fsalfirst = cache.du₁ + cache.du₂

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::IRKCConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p, fsalfirst = integrator
    @unpack minm, du₁, du₂, nlsolver = cache
    @unpack f1, f2 = integrator.f
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)

    # The the number of degree for Chebyshev polynomial
    #maxm = max(2,Int(floor(sqrt(integrator.opts.internalnorm(integrator.opts.reltol,t)/(10 *eps(integrator.opts.internalnorm(uprev,t)))))))
    maxm = 50
    mdeg = 1 + floor(Int, sqrt(1.54 * abs(dt) * integrator.eigen_est + 1))
    mdeg = min(maxm, max(minm, mdeg))

    ω₀ = 1 + 2 / (13 * (mdeg^2))
    temp₁ = ω₀^2 - 1
    temp₂ = sqrt(temp₁)
    θ = mdeg * log(ω₀ + temp₂)
    ω₁ = (sinh(θ) * temp₁) / (cosh(θ) * mdeg * temp₂ - ω₀ * sinh(θ))
    Bⱼ₋₂ = 1 / (4 * ω₀^2)
    Bⱼ₋₁ = 1 / ω₀

    #stage-1
    f1ⱼ₋₂ = du₁
    gprev2 = copy(uprev)
    μs = ω₁ * Bⱼ₋₁
    μs₁ = μs

    # initial guess for implicit part
    # if alg.extrapolant == :linear
    #   nlsolver.z = dt*du₁
    # else # :constant
    #   nlsolver.z = zero(u)
    # end

    nlsolver.z = dt * du₁

    nlsolver.tmp = uprev + dt * μs₁ * du₂
    nlsolver.γ = μs₁
    nlsolver.c = μs
    markfirststage!(nlsolver)
    z = nlsolve!(nlsolver, integrator, cache, false)
    # nlsolvefail(nlsolver) && return
    gprev = nlsolver.tmp + μs₁ * z

    Cⱼ₋₂ = zero(eltype(u))
    Cⱼ₋₁ = μs
    Tⱼ₋₁ = ω₀
    Tⱼ₋₂ = one(eltype(u))
    Tⱼ₋₁′ = one(eltype(u))
    Tⱼ₋₂′ = zero(eltype(u))
    Tⱼ₋₁″ = zero(eltype(u))
    Tⱼ₋₂″ = zero(eltype(u))

    #stage- 2...mdeg
    for iter in 2:mdeg
        Tⱼ = 2 * ω₀ * Tⱼ₋₁ - Tⱼ₋₂
        Tⱼ′ = 2 * ω₀ * Tⱼ₋₁′ + 2 * Tⱼ₋₁ - Tⱼ₋₂′
        Tⱼ″ = 2 * ω₀ * Tⱼ₋₁″ + 4 * Tⱼ₋₁′ - Tⱼ₋₂″
        Bⱼ = Tⱼ″ / (Tⱼ′^2)
        μ = (2 * ω₀ * Bⱼ) / Bⱼ₋₁
        ν = -Bⱼ / Bⱼ₋₂
        μs = (μ * ω₁) / ω₀
        νs = -(1 - Tⱼ₋₁ * Bⱼ₋₁) * μs
        Cⱼ = μ * Cⱼ₋₁ + ν * Cⱼ₋₂ + μs + νs

        f1ⱼ₋₁ = f1(gprev, p, t + Cⱼ₋₁ * dt)
        f2ⱼ₋₁ = f2(gprev, p, t + Cⱼ₋₁ * dt)
        integrator.stats.nf += 1
        integrator.stats.nf2 += 1
        nlsolver.tmp = (1 - μ - ν) * uprev + μ * gprev + ν * gprev2 + dt * μs * f2ⱼ₋₁ +
                       dt * νs * du₂ + (νs - (1 - μ - ν) * μs₁) * dt * du₁ -
                       ν * μs₁ * dt * f1ⱼ₋₂
        nlsolver.z = dt * f1ⱼ₋₁
        nlsolver.c = Cⱼ
        z = nlsolve!(nlsolver, integrator, cache, false)
        # ignoring newton method's convergence failure
        # nlsolvefail(nlsolver) && return
        u = nlsolver.tmp + μs₁ * z
        if (iter < mdeg)
            f1ⱼ₋₂ = f1ⱼ₋₁
            gprev2 = gprev
            gprev = u
            Cⱼ₋₂ = Cⱼ₋₁
            Cⱼ₋₁ = Cⱼ
            Bⱼ₋₂ = Bⱼ₋₁
            Bⱼ₋₁ = Bⱼ
            Tⱼ₋₂ = Tⱼ₋₁
            Tⱼ₋₁ = Tⱼ
            Tⱼ₋₂′ = Tⱼ₋₁′
            Tⱼ₋₁′ = Tⱼ′
            Tⱼ₋₂″ = Tⱼ₋₁″
            Tⱼ₋₁″ = Tⱼ″
        end
    end

    cache.du₁ = f1(u, p, t + dt)
    cache.du₂ = f2(u, p, t + dt)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
    # error estimate
    if isnewton(nlsolver) && integrator.opts.adaptive
        update_W!(integrator, cache, dt, false)
        tmp = dt * (0.5 * (cache.du₂ - du₂) + (0.5 - μs₁) * (cache.du₁ - du₁))
        tmp = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.fsallast = cache.du₁ + cache.du₂
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::IRKCCache)
    @unpack uprev, p, t = integrator
    @unpack f1, f2 = integrator.f
    integrator.kshortsize = 2
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = du_alias_or_new(cache.nlsolver, integrator.fsalfirst)
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    f1(cache.du₁, uprev, p, t)
    f2(cache.du₂, uprev, p, t)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
    @.. broadcast=false integrator.fsalfirst=cache.du₁ + cache.du₂
end

function perform_step!(integrator, cache::IRKCCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack gprev, gprev2, f1ⱼ₋₁, f1ⱼ₋₂, f2ⱼ₋₁, du₁, du₂, atmp, nlsolver = cache
    @unpack tmp, z = nlsolver
    @unpack minm = cache.constantcache
    @unpack f1, f2 = integrator.f

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    # The the number of degree for Chebyshev polynomial
    #maxm = max(2,int(floor(sqrt(integrator.opts.internalnorm(integrator.opts.reltol,t)/(10 *eps(integrator.opts.internalnorm(uprev,t)))))))
    maxm = 50
    mdeg = 1 + Int(floor(sqrt(1.54 * abs(dt) * integrator.eigen_est + 1)))
    mdeg = (mdeg < minm) ? minm : mdeg
    mdeg = (mdeg >= maxm) ? maxm : mdeg

    ω₀ = 1 + 2 / (13 * (mdeg^2))
    temp₁ = ω₀^2 - 1
    temp₂ = sqrt(temp₁)
    θ = mdeg * log(ω₀ + temp₂)
    ω₁ = (sinh(θ) * temp₁) / (cosh(θ) * mdeg * temp₂ - ω₀ * sinh(θ))
    Bⱼ₋₂ = 1 / (4 * ω₀^2)
    Bⱼ₋₁ = 1 / ω₀

    #stage-1
    f1ⱼ₋₂ = du₁
    @.. broadcast=false gprev2=uprev
    μs = ω₁ * Bⱼ₋₁
    μs₁ = μs

    # initial guess
    # if alg.extrapolant == :linear
    #   @.. broadcast=false z = dt*du₁
    # else # :constant
    #   @.. broadcast=false z = zero(eltype(u))
    # end
    @.. broadcast=false nlsolver.z=dt * du₁

    @.. broadcast=false nlsolver.tmp=uprev + dt * μs₁ * du₂
    nlsolver.γ = μs₁
    nlsolver.c = μs
    markfirststage!(nlsolver)
    z = nlsolve!(nlsolver, integrator, cache, false)
    # ignoring newton method's convergence failure
    # nlsolvefail(nlsolver) && return
    @.. broadcast=false gprev=nlsolver.tmp + μs₁ * nlsolver.z

    Cⱼ₋₂ = zero(eltype(u))
    Cⱼ₋₁ = μs
    Tⱼ₋₁ = ω₀
    Tⱼ₋₂ = one(eltype(u))
    Tⱼ₋₁′ = one(eltype(u))
    Tⱼ₋₂′ = zero(eltype(u))
    Tⱼ₋₁″ = zero(eltype(u))
    Tⱼ₋₂″ = zero(eltype(u))

    #stage- 2...mdeg
    for iter in 2:mdeg
        Tⱼ = 2 * ω₀ * Tⱼ₋₁ - Tⱼ₋₂
        Tⱼ′ = 2 * ω₀ * Tⱼ₋₁′ + 2 * Tⱼ₋₁ - Tⱼ₋₂′
        Tⱼ″ = 2 * ω₀ * Tⱼ₋₁″ + 4 * Tⱼ₋₁′ - Tⱼ₋₂″
        Bⱼ = Tⱼ″ / (Tⱼ′^2)
        μ = (2 * ω₀ * Bⱼ) / Bⱼ₋₁
        ν = -Bⱼ / Bⱼ₋₂
        μs = (μ * ω₁) / ω₀
        νs = -(1 - Tⱼ₋₁ * Bⱼ₋₁) * μs
        Cⱼ = μ * Cⱼ₋₁ + ν * Cⱼ₋₂ + μs + νs

        f1(f1ⱼ₋₁, gprev, p, t + Cⱼ₋₁ * dt)
        f2(f2ⱼ₋₁, gprev, p, t + Cⱼ₋₁ * dt)
        integrator.stats.nf += 1
        integrator.stats.nf2 += 1
        @.. broadcast=false nlsolver.tmp=(1 - μ - ν) * uprev + μ * gprev + ν * gprev2 +
                                         dt * μs * f2ⱼ₋₁ + dt * νs * du₂ +
                                         (νs - (1 - μ - ν) * μs₁) * dt * du₁ -
                                         ν * μs₁ * dt * f1ⱼ₋₂
        @.. broadcast=false nlsolver.z=dt * f1ⱼ₋₁
        nlsolver.c = Cⱼ

        z = nlsolve!(nlsolver, integrator, cache, false)
        # nlsolvefail(nlsolver) && return
        @.. broadcast=false u=nlsolver.tmp + μs₁ * nlsolver.z
        if (iter < mdeg)
            @.. broadcast=false f1ⱼ₋₂=f1ⱼ₋₁
            @.. broadcast=false gprev2=gprev
            @.. broadcast=false gprev=u
            Cⱼ₋₂ = Cⱼ₋₁
            Cⱼ₋₁ = Cⱼ
            Bⱼ₋₂ = Bⱼ₋₁
            Bⱼ₋₁ = Bⱼ
            Tⱼ₋₂ = Tⱼ₋₁
            Tⱼ₋₁ = Tⱼ
            Tⱼ₋₂′ = Tⱼ₋₁′
            Tⱼ₋₁′ = Tⱼ′
            Tⱼ₋₂″ = Tⱼ₋₁″
            Tⱼ₋₁″ = Tⱼ″
        end
    end

    @.. broadcast=false f1ⱼ₋₁=du₁
    @.. broadcast=false f2ⱼ₋₁=du₂
    f1(du₁, u, p, t + dt)
    f2(du₂, u, p, t + dt)
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
    # error estimate
    if isnewton(nlsolver) && integrator.opts.adaptive
        update_W!(integrator, cache, dt, false)
        @.. broadcast=false gprev=dt * 0.5 * (du₂ - f2ⱼ₋₁) +
                                  dt * (0.5 - μs₁) * (du₁ - f1ⱼ₋₁)

        linsolve = nlsolver.cache.linsolve
        linres = dolinsolve(integrator, linsolve; b = _vec(gprev), linu = _vec(tmp))

        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    @.. broadcast=false integrator.fsallast=du₁ + du₂
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end