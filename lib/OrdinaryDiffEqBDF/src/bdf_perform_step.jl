function initialize!(integrator, cache::ABDF2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::ABDF2ConstantCache, repeat_step = false)
    (; t, f, p) = integrator
    (; dtₙ₋₁, nlsolver) = cache
    alg = unwrap_alg(integrator, true)
    dtₙ, uₙ, uₙ₋₁, uₙ₋₂ = integrator.dt, integrator.u, integrator.uprev, integrator.uprev2

    # TODO: this doesn't look right
    if integrator.iter == 1 && !integrator.u_modified
        cache.dtₙ₋₁ = dtₙ
        cache.eulercache.nlsolver.method = DIRK
        perform_step!(integrator, cache.eulercache, repeat_step)
        cache.fsalfirstprev = integrator.fsalfirst
        return
    end

    fₙ₋₁ = integrator.fsalfirst
    ρ = dtₙ / dtₙ₋₁
    β₀ = Int64(2) // 3
    β₁ = -(ρ - 1) / 3
    α₀ = 1
    α̂ = ρ^2 / 3
    α₁ = 1 + α̂
    α₂ = -α̂

    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        u = @.. broadcast = false uₙ₋₁ + dtₙ * fₙ₋₁
    else # :constant
        u = uₙ₋₁
    end
    nlsolver.z = u

    mass_matrix = f.mass_matrix

    if mass_matrix === I
        nlsolver.tmp = @.. (
            (dtₙ * β₁) * fₙ₋₁ +
                (α₁ * uₙ₋₁ + α₂ * uₙ₋₂)
        ) / (dtₙ * β₀)
    else
        _tmp = mass_matrix * @.. (α₁ * uₙ₋₁ + α₂ * uₙ₋₂)
        nlsolver.tmp = @.. ((dtₙ * β₁) * fₙ₋₁ + _tmp) / (dtₙ * β₀)
    end
    nlsolver.γ = β₀
    nlsolver.α = α₀
    nlsolver.method = COEFFICIENT_MULTISTEP
    uₙ = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    integrator.fsallast = f(uₙ, p, t + dtₙ)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        tmp = integrator.fsallast - (1 + dtₙ / dtₙ₋₁) * integrator.fsalfirst +
            (dtₙ / dtₙ₋₁) * cache.fsalfirstprev
        est = (dtₙ₋₁ + dtₙ) / 6 * tmp
        atmp = calculate_residuals(
            est, uₙ₋₁, uₙ, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    ################################### Finalize

    if integrator.EEst < one(integrator.EEst)
        cache.fsalfirstprev = integrator.fsalfirst
        cache.dtₙ₋₁ = dtₙ
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = uₙ
    return
end

function initialize!(integrator, cache::ABDF2Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::ABDF2Cache, repeat_step = false)
    (; t, dt, f, p) = integrator
    #TODO: remove zₙ₋₁ from the cache
    (; atmp, dtₙ₋₁, zₙ₋₁, nlsolver, step_limiter!) = cache
    (; z, tmp, ztmp) = nlsolver
    alg = unwrap_alg(integrator, true)
    uₙ, uₙ₋₁, uₙ₋₂, dtₙ = integrator.u, integrator.uprev, integrator.uprev2, integrator.dt

    if integrator.iter == 1 && !integrator.u_modified
        cache.dtₙ₋₁ = dtₙ
        cache.eulercache.nlsolver.method = DIRK
        perform_step!(integrator, cache.eulercache, repeat_step)
        cache.fsalfirstprev .= integrator.fsalfirst
        nlsolver.tmp = tmp
        return
    end

    fₙ₋₁ = integrator.fsalfirst
    ρ = dtₙ / dtₙ₋₁
    β₀ = Int64(2) // 3
    β₁ = -(ρ - 1) / 3
    α₀ = 1
    α̂ = ρ^2 / 3
    α₁ = 1 + α̂
    α₂ = -α̂

    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        @.. broadcast = false z = uₙ₋₁ + dtₙ * fₙ₋₁
    else # :constant
        @.. broadcast = false z = uₙ₋₁
    end

    mass_matrix = f.mass_matrix

    if mass_matrix === I
        @.. broadcast = false tmp = ((dtₙ * β₁) * fₙ₋₁ + (α₁ * uₙ₋₁ + α₂ * uₙ₋₂)) / (dtₙ * β₀)
    else
        dz = nlsolver.cache.dz
        @.. broadcast = false dz = α₁ * uₙ₋₁ + α₂ * uₙ₋₂
        mul!(ztmp, mass_matrix, dz)
        @.. broadcast = false tmp = ((dtₙ * β₁) * fₙ₋₁ + ztmp) / (dtₙ * β₀)
    end
    nlsolver.γ = β₀
    nlsolver.α = α₀
    nlsolver.method = COEFFICIENT_MULTISTEP
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast = false uₙ = z

    step_limiter!(uₙ, integrator, p, t + dtₙ)

    f(integrator.fsallast, uₙ, p, t + dtₙ)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    if integrator.opts.adaptive
        btilde0 = (dtₙ₋₁ + dtₙ) * Int64(1) // 6
        btilde1 = 1 + dtₙ / dtₙ₋₁
        btilde2 = dtₙ / dtₙ₋₁
        @.. broadcast = false tmp = btilde0 *
            (
            integrator.fsallast - btilde1 * integrator.fsalfirst +
                btilde2 * cache.fsalfirstprev
        )
        calculate_residuals!(
            atmp, tmp, uₙ₋₁, uₙ, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    ################################### Finalize

    if integrator.EEst < one(integrator.EEst)
        @.. broadcast = false cache.fsalfirstprev = integrator.fsalfirst
        cache.dtₙ₋₁ = dtₙ
    end
    return
end

# SBDF

function initialize!(integrator, cache::SBDFConstantCache)
    (; uprev, p, t) = integrator
    (; f1, f2) = integrator.f
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    cache.du₁ = f1(uprev, p, t)
    cache.du₂ = f2(uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.fsalfirst = cache.du₁ + cache.du₂

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::SBDFConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    alg = unwrap_alg(integrator, true)
    (; uprev2, uprev3, uprev4, du₁, du₂, k₁, k₂, k₃, nlsolver) = cache
    (; f1, f2) = integrator.f
    cnt = cache.cnt = min(alg.order, integrator.iter + 1)
    integrator.iter == 1 && !integrator.u_modified && (cnt = cache.cnt = 1)
    nlsolver.γ = γ = inv(γₖ[cnt])
    if cache.ark
        # Additive Runge-Kutta Method
        du₂ = f2(uprev + dt * du₁, p, t)
        integrator.stats.nf2 += 1
    end
    if cnt == 1
        tmp = uprev + dt * du₂
    elseif cnt == 2
        tmp = γ * (2 * uprev - Int64(1) // 2 * uprev2 + dt * (2 * du₂ - k₁))
    elseif cnt == 3
        tmp = γ *
            (3 * uprev - Int64(3) // 2 * uprev2 + Int64(1) // 3 * uprev3 + dt * (3 * (du₂ - k₁) + k₂))
    else
        tmp = γ * (
            4 * uprev - 3 * uprev2 + Int64(4) // 3 * uprev3 - Int64(1) // 4 * uprev4 +
                dt * (4 * du₂ - 6 * k₁ + 4 * k₂ - k₃)
        )
    end
    nlsolver.tmp = tmp

    # Implicit part
    # precalculations
    γdt = γ * dt
    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        z = dt * du₁
    else # :constant
        z = zero(u)
    end
    nlsolver.z = z

    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = nlsolver.tmp + γ * z

    cnt == 4 && (
        cache.uprev4 = uprev3;
        cache.k₃ = k₂
    )
    cnt >= 3 && (
        cache.uprev3 = uprev2;
        cache.k₂ = k₁
    )
    (
        cache.uprev2 = uprev;
        cache.k₁ = du₂
    )
    cache.du₁ = f1(u, p, t + dt)
    cache.du₂ = f2(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.fsallast = cache.du₁ + cache.du₂
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function initialize!(integrator, cache::SBDFCache)
    (; uprev, p, t) = integrator
    (; f1, f2) = integrator.f
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    f1(cache.du₁, uprev, p, t)
    f2(cache.du₂, uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    return @.. broadcast = false integrator.fsalfirst = cache.du₁ + cache.du₂
end

function perform_step!(integrator, cache::SBDFCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    alg = unwrap_alg(integrator, true)
    (; uprev2, uprev3, uprev4, k₁, k₂, k₃, du₁, du₂, nlsolver) = cache
    (; tmp, z) = nlsolver
    (; f1, f2) = integrator.f
    cnt = cache.cnt = min(alg.order, integrator.iter + 1)
    integrator.iter == 1 && !integrator.u_modified && (cnt = cache.cnt = 1)
    nlsolver.γ = γ = inv(γₖ[cnt])
    # Explicit part
    if cache.ark
        # Additive Runge-Kutta Method
        f2(du₂, uprev + dt * du₁, p, t)
        integrator.stats.nf2 += 1
    end
    if cnt == 1
        @.. broadcast = false tmp = uprev + dt * du₂
    elseif cnt == 2
        @.. broadcast = false tmp = γ * (2 * uprev - 1 // 2 * uprev2 + dt * (2 * du₂ - k₁))
    elseif cnt == 3
        @.. broadcast = false tmp = γ * (
            3 * uprev - 3 // 2 * uprev2 + 1 // 3 * uprev3 +
                dt * (3 * (du₂ - k₁) + k₂)
        )
    else
        @.. broadcast = false tmp = γ * (
            4 * uprev - 3 * uprev2 + 4 // 3 * uprev3 -
                1 // 4 * uprev4 + dt * (4 * du₂ - 6 * k₁ + 4 * k₂ - k₃)
        )
    end
    # Implicit part
    # precalculations
    γdt = γ * dt
    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        @.. broadcast = false z = dt * du₁
    else # :constant
        @.. broadcast = false z = zero(eltype(u))
    end

    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast = false u = tmp + γ * z

    cnt == 4 && (
        cache.uprev4 .= uprev3;
        cache.k₃ .= k₂
    )
    cnt >= 3 && (
        cache.uprev3 .= uprev2;
        cache.k₂ .= k₁
    )
    (
        cache.uprev2 .= uprev;
        cache.k₁ .= du₂
    )
    f1(du₁, u, p, t + dt)
    f2(du₂, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    return @.. broadcast = false integrator.fsallast = du₁ + du₂
end

# QNDF1

function initialize!(integrator, cache::QNDF1ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::QNDF1ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; uprev2, D, D2, R, U, dtₙ₋₁, nlsolver) = cache
    alg = unwrap_alg(integrator, true)
    κ = alg.kappa
    cnt = integrator.iter
    k = 1
    if cnt > 1
        ρ = dt / dtₙ₋₁
        D[1] = uprev - uprev2   # backward diff
        if ρ != 1
            R!(k, ρ, cache)
            R .= R * U
            D[1] = D[1] * R[1, 1]
        end
    else
        κ = zero(alg.kappa)
    end

    #Changing uprev2 after D Array has changed with step-size
    uprev2 = uprev - D[1]

    β₀ = 1
    α₀ = 1 - κ
    α₁ = 1 - 2 * κ
    α₂ = κ

    markfirststage!(nlsolver)

    # initial guess
    nlsolver.z = uprev + sum(D)

    mass_matrix = f.mass_matrix

    if mass_matrix === I
        nlsolver.tmp = @.. broadcast = false (α₁ * uprev + α₂ * uprev2) / (dt * β₀)
    else
        _tmp = mass_matrix * @.. broadcast = false (α₁ * uprev + α₂ * uprev2)
        nlsolver.tmp = @.. broadcast = false _tmp / (dt * β₀)
    end

    nlsolver.γ = β₀
    nlsolver.α = α₀
    nlsolver.method = COEFFICIENT_MULTISTEP

    u = nlsolve!(nlsolver, integrator, cache, repeat_step)

    nlsolvefail(nlsolver) && return
    if integrator.opts.adaptive
        if integrator.success_iter == 0
            integrator.EEst = one(integrator.EEst)
        else
            D2[1] = u - uprev
            D2[2] = D2[1] - D[1]
            utilde = (κ + inv(k + 1)) * D2[2]
            atmp = calculate_residuals(
                utilde, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t
            )
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        end
    end
    if integrator.EEst > one(integrator.EEst)
        return
    end
    cache.dtₙ₋₁ = dt
    cache.uprev2 = uprev
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return
end

function initialize!(integrator, cache::QNDF1Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::QNDF1Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; uprev2, D, D2, R, U, dtₙ₋₁, utilde, atmp, nlsolver, step_limiter!) = cache
    (; z, tmp, ztmp) = nlsolver
    alg = unwrap_alg(integrator, true)
    κ = alg.kappa
    cnt = integrator.iter
    k = 1
    if cnt > 1
        ρ = dt / dtₙ₋₁
        @.. broadcast = false D[1] = uprev - uprev2 # backward diff
        if ρ != 1
            R!(k, ρ, cache)
            R .= R * U
            @.. broadcast = false D[1] = D[1] * R[1, 1]
        end
    else
        κ = zero(alg.kappa)
    end

    #Changing uprev2 after D Array has changed with step-size
    uprev2 = uprev - D[1]

    β₀ = 1
    α₀ = 1 - κ
    α₁ = 1 - 2 * κ
    α₂ = κ

    markfirststage!(nlsolver)

    # initial guess
    nlsolver.z = uprev + sum(D)

    mass_matrix = f.mass_matrix

    if mass_matrix === I
        @.. broadcast = false tmp = (α₁ * uprev + α₂ * uprev2) / (dt * β₀)
    else
        dz = nlsolver.cache.dz
        @.. broadcast = false dz = α₁ * uprev + α₂ * uprev2
        mul!(ztmp, mass_matrix, dz)
        @.. broadcast = false tmp = ztmp / (dt * β₀)
    end

    nlsolver.γ = β₀
    nlsolver.α = α₀
    nlsolver.method = COEFFICIENT_MULTISTEP
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast = false u = z

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        if integrator.success_iter == 0
            integrator.EEst = one(integrator.EEst)
        else
            @.. broadcast = false D2[1] = u - uprev
            @.. broadcast = false D2[2] = D2[1] - D[1]
            @.. broadcast = false utilde = (κ + inv(k + 1)) * D2[2]
            calculate_residuals!(
                atmp, utilde, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        end
    end
    if integrator.EEst > one(integrator.EEst)
        return
    end
    cache.uprev2 .= uprev
    cache.dtₙ₋₁ = dt
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return
end

function initialize!(integrator, cache::QNDF2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::QNDF2ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; uprev2, uprev3, dtₙ₋₁, dtₙ₋₂, D, D2, R, U, nlsolver) = cache
    alg = unwrap_alg(integrator, true)
    cnt = integrator.iter
    k = 2
    if cnt == 1 || cnt == 2
        κ = zero(alg.kappa)
        γ₁ = Int64(1) // 1
        γ₂ = Int64(1) // 1
    elseif dtₙ₋₁ != dtₙ₋₂
        κ = alg.kappa
        γ₁ = Int64(1) // 1
        γ₂ = Int64(1) // 1 + Int64(1) // 2
        ρ₁ = dt / dtₙ₋₁
        ρ₂ = dt / dtₙ₋₂
        D[1] = uprev - uprev2
        D[1] = D[1] * ρ₁
        D[2] = D[1] - ((uprev2 - uprev3) * ρ₂)
    else
        κ = alg.kappa
        γ₁ = Int64(1) // 1
        γ₂ = Int64(1) // 1 + Int64(1) // 2
        ρ = dt / dtₙ₋₁
        # backward diff
        D[1] = uprev - uprev2
        D[2] = D[1] - (uprev2 - uprev3)
        if ρ != 1
            R!(k, ρ, cache)
            R .= R * U
            D[1] = D[1] * R[1, 1] + D[2] * R[2, 1]
            D[2] = D[1] * R[1, 2] + D[2] * R[2, 2]
        end
    end

    β₀ = inv((1 - κ) * γ₂)
    α₀ = 1

    u₀ = uprev + D[1] + D[2]
    ϕ = (γ₁ * D[1] + γ₂ * D[2]) * β₀

    markfirststage!(nlsolver)

    # initial guess
    nlsolver.z = uprev + sum(D)

    mass_matrix = f.mass_matrix

    if mass_matrix === I
        nlsolver.tmp = @.. broadcast = false (u₀ - ϕ) / (dt * β₀)
    else
        _tmp = mass_matrix * @.. broadcast = false (u₀ - ϕ)
        nlsolver.tmp = @.. broadcast = false _tmp / (dt * β₀)
    end

    nlsolver.γ = β₀
    nlsolver.α = α₀
    nlsolver.method = COEFFICIENT_MULTISTEP

    u = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    if integrator.opts.adaptive
        if integrator.success_iter == 0
            integrator.EEst = one(integrator.EEst)
        elseif integrator.success_iter == 1
            utilde = (u - uprev) - ((uprev - uprev2) * dt / dtₙ₋₁)
            atmp = calculate_residuals(
                utilde, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t
            )
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        else
            D2[1] = u - uprev
            D2[2] = D2[1] - D[1]
            D2[3] = D2[2] - D[2]
            utilde = (κ * γ₂ + inv(k + 1)) * D2[3]
            atmp = calculate_residuals(
                utilde, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t
            )
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        end
    end
    if integrator.EEst > one(integrator.EEst)
        return
    end

    cache.uprev3 = uprev2
    cache.uprev2 = uprev
    cache.dtₙ₋₂ = dtₙ₋₁
    cache.dtₙ₋₁ = dt
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return
end

function initialize!(integrator, cache::QNDF2Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::QNDF2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (;
        uprev2, uprev3, dtₙ₋₁, dtₙ₋₂, D, D2, R, U, utilde, atmp, nlsolver,
        step_limiter!,
    ) = cache
    (; z, tmp, ztmp) = nlsolver
    alg = unwrap_alg(integrator, true)
    cnt = integrator.iter
    k = 2
    if cnt == 1 || cnt == 2
        κ = zero(alg.kappa)
        γ₁ = Int64(1) // 1
        γ₂ = Int64(1) // 1
    elseif dtₙ₋₁ != dtₙ₋₂
        κ = alg.kappa
        γ₁ = Int64(1) // 1
        γ₂ = Int64(1) // 1 + Int64(1) // 2
        ρ₁ = dt / dtₙ₋₁
        ρ₂ = dt / dtₙ₋₂
        @.. broadcast = false D[1] = uprev - uprev2
        @.. broadcast = false D[1] = D[1] * ρ₁
        @.. broadcast = false D[2] = D[1] - ((uprev2 - uprev3) * ρ₂)
    else
        κ = alg.kappa
        γ₁ = Int64(1) // 1
        γ₂ = Int64(1) // 1 + Int64(1) // 2
        ρ = dt / dtₙ₋₁
        # backward diff
        @.. broadcast = false D[1] = uprev - uprev2
        @.. broadcast = false D[2] = D[1] - (uprev2 - uprev3)
        if ρ != 1
            R!(k, ρ, cache)
            R .= R * U
            @.. broadcast = false D[1] = D[1] * R[1, 1] + D[2] * R[2, 1]
            @.. broadcast = false D[2] = D[1] * R[1, 2] + D[2] * R[2, 2]
        end
    end

    β₀ = inv((1 - κ) * γ₂)
    α₀ = 1

    u₀ = uprev + D[1] + D[2]
    ϕ = (γ₁ * D[1] + γ₂ * D[2]) * β₀

    markfirststage!(nlsolver)

    # initial guess
    nlsolver.z = uprev + sum(D)

    mass_matrix = f.mass_matrix

    if mass_matrix === I
        @.. broadcast = false tmp = (u₀ - ϕ) / (dt * β₀)
    else
        dz = nlsolver.cache.dz
        @.. broadcast = false dz = u₀ - ϕ
        mul!(ztmp, mass_matrix, dz)
        @.. broadcast = false tmp = ztmp / (dt * β₀)
    end

    nlsolver.γ = β₀
    nlsolver.α = α₀
    nlsolver.method = COEFFICIENT_MULTISTEP

    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast = false u = z

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        if integrator.success_iter == 0
            integrator.EEst = one(integrator.EEst)
        elseif integrator.success_iter == 1
            @.. broadcast = false utilde = (u - uprev) - ((uprev - uprev2) * dt / dtₙ₋₁)
            calculate_residuals!(
                atmp, utilde, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        else
            @.. broadcast = false D2[1] = u - uprev
            @.. broadcast = false D2[2] = D2[1] - D[1]
            @.. broadcast = false D2[3] = D2[2] - D[2]
            @.. broadcast = false utilde = (κ * γ₂ + inv(k + 1)) * D2[3]
            calculate_residuals!(
                atmp, utilde, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        end
    end
    if integrator.EEst > one(integrator.EEst)
        return
    end

    cache.uprev3 .= uprev2
    cache.uprev2 .= uprev
    cache.dtₙ₋₂ = dtₙ₋₁
    cache.dtₙ₋₁ = dt
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return
end

function initialize!(integrator, cache::QNDFConstantCache{max_order}) where {max_order}
    integrator.kshortsize = max_order
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    for i in 1:max_order
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    return nothing
end

function perform_step!(
        integrator, cache::QNDFConstantCache{max_order},
        repeat_step = false
    ) where {max_order}
    (; t, dt, uprev, u, f, p) = integrator
    (; dtprev, order, D, U, nlsolver, γₖ) = cache
    alg = unwrap_alg(integrator, true)

    if integrator.u_modified
        dtprev = one(dt)
        order = 1
        cache.nconsteps = 0
        cache.consfailcnt = 0
        for i in eachindex(D)
            D[i] = zero(D[i])
        end
        for i in eachindex(cache.prevD)
            cache.prevD[i] = zero(cache.prevD[i])
        end
    end

    k = order
    κlist = alg.kappa
    κ = κlist[k]
    if cache.consfailcnt > 0
        # Deep copy to avoid aliasing: D[i] and prevD[i] must not share arrays
        for i in eachindex(D)
            D[i] = copy(cache.prevD[i])
        end
    end
    if dt != dtprev || cache.prevorder != k
        ρ = dt / dtprev
        integrator.cache.nconsteps = 0
        (; U) = cache
        R = calc_R(ρ, k, Val(max_order))
        RU = R * U
        D_new = [sum(D[i] * RU[i, j] for i in 1:k) for j in 1:k]
        for j in 1:k
            D[j] = D_new[j]
        end
    end

    α₀ = 1
    β₀ = inv((1 - κ) * γₖ[k])
    if u isa Number
        u₀ = sum(D[i] for i in 1:k) + uprev
        ϕ = zero(u)
        for i in 1:k
            ϕ += γₖ[i] * D[i]
        end
    else
        u₀ = sum(D[i] for i in 1:k) .+ uprev
        ϕ = zero(u)
        for i in 1:k
            ϕ = @.. ϕ + γₖ[i] * D[i]
        end
    end
    markfirststage!(nlsolver)
    nlsolver.z = u₀
    mass_matrix = f.mass_matrix

    if mass_matrix === I
        nlsolver.tmp = @.. (u₀ / β₀ - ϕ) / dt
    else
        nlsolver.tmp = mass_matrix * @.. (u₀ / β₀ - ϕ) / dt
    end

    nlsolver.γ = β₀
    nlsolver.α = α₀
    nlsolver.method = COEFFICIENT_MULTISTEP

    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = z
    dd = u - u₀
    update_D!(D, dd, k)

    if integrator.opts.adaptive
        (; abstol, reltol, internalnorm) = integrator.opts
        if cache.consfailcnt > 1 && mass_matrix !== I
            # if we get repeated failure and mass_matrix !== I it's likely that
            # there's a discontinuity on the algebraic equations
            atmp = calculate_residuals(
                mass_matrix * dd, uprev, u, abstol, reltol,
                internalnorm, t
            )
        else
            atmp = calculate_residuals(dd, uprev, u, abstol, reltol, internalnorm, t)
        end
        integrator.EEst = error_constant(integrator, k) * internalnorm(atmp, t)
        if k > 1
            atmpm1 = calculate_residuals(
                D[k],
                uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            cache.EEst1 = error_constant(integrator, k - 1) * internalnorm(atmpm1, t)
        end
        if k < max_order
            atmpp1 = calculate_residuals(
                D[k + 2],
                uprev, u, abstol, reltol, internalnorm, t
            )
            cache.EEst2 = error_constant(integrator, k + 1) * internalnorm(atmpp1, t)
        end
    end
    if integrator.EEst <= one(integrator.EEst)
        # Deep copy to avoid aliasing: prevD[i] must not share arrays with D[i]
        for i in eachindex(D)
            cache.prevD[i] = copy(D[i])
        end
        cache.dtprev = dt
        cache.prevorder = k
        if integrator.opts.dense
            integrator.fsallast = f(u, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
    end
    if integrator.opts.calck
        for j in 1:max_order
            if j <= k
                # Deep copy: k[j] must not alias D[j] since update_D! mutates in-place
                integrator.k[j] = copy(D[j])
            else
                integrator.k[j] = zero(u)
            end
        end
    end
    return integrator.u = u
end

function initialize!(integrator, cache::QNDFCache{max_order}) where {max_order}
    integrator.kshortsize = max_order

    resize!(integrator.k, integrator.kshortsize)
    for i in 1:max_order
        integrator.k[i] = cache.dense[i]
    end
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(
        integrator, cache::QNDFCache{max_order},
        repeat_step = false
    ) where {max_order}
    (; t, dt, uprev, u, f, p) = integrator
    (;
        dtprev, order, D, nlsolver, γₖ, dd, atmp, atmpm1, atmpp1,
        utilde, utildem1, utildep1, ϕ, u₀, step_limiter!,
    ) = cache
    alg = unwrap_alg(integrator, true)

    if integrator.u_modified
        dtprev = one(dt)
        order = 1
        cache.nconsteps = 0
        cache.consfailcnt = 0
        for d in D
            recursivefill!(d, false)
        end
        for d in cache.prevD
            recursivefill!(d, false)
        end
    end

    k = order
    κlist = alg.kappa
    κ = κlist[k]
    if cache.consfailcnt > 0
        for i in eachindex(D)
            copyto!(D[i], cache.prevD[i])
        end
    end
    if dt != dtprev || cache.prevorder != k
        ρ = dt / dtprev
        cache.nconsteps = 0
        (; RU, U, Dtmp) = cache
        R = calc_R(ρ, k, Val(max_order))
        copyto!(RU, R * U)
        for j in 1:k
            fill!(Dtmp[j], zero(eltype(Dtmp[j])))
            for i in 1:k
                @.. broadcast = false Dtmp[j] += D[i] * RU[i, j]
            end
        end
        D, Dtmp = Dtmp, D
        cache.D = D
        cache.Dtmp = Dtmp
    end

    α₀ = 1
    β₀ = inv((1 - κ) * γₖ[k])
    # D[j] contains scaled j-th derivative approximation.
    # Thus, it’s likely that ||D[j+1]|| <= ||D[j]|| holds,
    # we want to sum small numbers first to minimize the accumulation error.
    # Hence, we sum it backwards.
    @.. broadcast = false u₀ = D[k]
    for j in (k - 1):-1:1
        @.. broadcast = false u₀ += D[j]
    end
    @.. broadcast = false u₀ += uprev
    @.. broadcast = false ϕ = zero(u)
    for i in 1:k
        @.. broadcast = false ϕ += γₖ[i] * D[i]
    end
    markfirststage!(nlsolver)
    @.. broadcast = false nlsolver.z = u₀
    mass_matrix = f.mass_matrix

    if mass_matrix === I
        @.. broadcast = false nlsolver.tmp = (u₀ / β₀ - ϕ) / dt
    else
        (; tmp2) = cache
        @.. broadcast = false tmp2 = (u₀ / β₀ - ϕ) / dt
        mul!(nlsolver.tmp, mass_matrix, tmp2)
    end

    nlsolver.γ = β₀
    nlsolver.α = α₀
    nlsolver.method = COEFFICIENT_MULTISTEP

    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast = false u = z
    @.. broadcast = false dd = u - u₀
    update_D!(D, dd, k)

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        (; abstol, reltol, internalnorm) = integrator.opts
        if cache.consfailcnt > 1 && mass_matrix !== I
            # if we get repeated failure and mass_matrix !== I it's likely that
            # there's a discontinuity on the algebraic equations
            calculate_residuals!(
                atmp, mul!(nlsolver.tmp, mass_matrix, dd), uprev, u,
                abstol, reltol, internalnorm, t
            )
        else
            calculate_residuals!(atmp, dd, uprev, u, abstol, reltol, internalnorm, t)
        end
        integrator.EEst = error_constant(integrator, k) * internalnorm(atmp, t)
        if k > 1
            calculate_residuals!(
                atmpm1, D[k], uprev, u, abstol,
                reltol, internalnorm, t
            )
            cache.EEst1 = error_constant(integrator, k - 1) * internalnorm(atmpm1, t)
        end
        if k < max_order
            calculate_residuals!(
                atmpp1, D[k + 2], uprev, u, abstol,
                reltol, internalnorm, t
            )
            cache.EEst2 = error_constant(integrator, k + 1) * internalnorm(atmpp1, t)
        end
    end
    if integrator.EEst <= one(integrator.EEst)
        for i in eachindex(D)
            copyto!(cache.prevD[i], D[i])
        end
        cache.dtprev = dt
        cache.prevorder = k
        if integrator.opts.dense
            f(integrator.fsallast, u, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
    end
    if integrator.opts.calck
        for j in 1:length(integrator.k)
            if j <= k
                copyto!(integrator.k[j], D[j])
            else
                fill!(integrator.k[j], zero(eltype(u)))
            end
        end
    end
    return nothing
end

### MEBDF2
function initialize!(integrator, cache::MEBDF2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::MEBDF2ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    nlsolver = cache.nlsolver
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        nlsolver.z = dt * integrator.fsalfirst
    else # :constant
        nlsolver.z = zero(u)
    end

    ### STEP 1
    nlsolver.tmp = uprev
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    z₁ = nlsolver.tmp + z
    ### STEP 2
    nlsolver.tmp = z₁
    nlsolver.c = 2
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    z₂ = z₁ + z
    ### STEP 3
    tmp2 = 0.5uprev + z₁ - 0.5z₂
    nlsolver.tmp = tmp2
    nlsolver.c = 1
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = tmp2 + z

    ### finalize
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::MEBDF2Cache)
    integrator.kshortsize = 2

    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::MEBDF2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; z₁, z₂, tmp2, nlsolver) = cache
    z = nlsolver.z
    mass_matrix = integrator.f.mass_matrix
    alg = unwrap_alg(integrator, true)
    markfirststage!(nlsolver)

    # initial guess
    if alg.extrapolant == :linear
        @.. broadcast = false z = dt * integrator.fsalfirst
    else # :constant
        z .= zero(eltype(u))
    end

    ### STEP 1
    nlsolver.tmp = uprev
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast = false z₁ = uprev + z
    ### STEP 2
    nlsolver.tmp = z₁
    nlsolver.c = 2
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast = false z₂ = z₁ + z
    ### STEP 3
    # z .= zero(eltype(u))
    @.. broadcast = false tmp2 = 0.5uprev + z₁ - 0.5z₂
    nlsolver.tmp = tmp2
    nlsolver.c = 1
    isnewton(nlsolver) && set_new_W!(nlsolver, false)
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast = false u = tmp2 + z

    ### finalize
    f(integrator.fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator, cache::FBDFConstantCache{max_order}) where {max_order}
    integrator.kshortsize = max_order + 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # k[1..max_order+1] = solution values at fixed Chebyshev reference nodes
    integrator.fsallast = zero(integrator.fsalfirst)
    for i in 1:(max_order + 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end

    u_modified = integrator.u_modified
    integrator.u_modified = true
    reinitFBDF!(integrator, cache)
    return integrator.u_modified = u_modified
end

function perform_step!(
        integrator, cache::FBDFConstantCache{max_order},
        repeat_step = false
    ) where {max_order}
    (;
        ts, u_history, order, u_corrector, bdf_coeffs, r, nlsolver,
        ts_tmp, iters_from_event, nconsteps,
    ) = cache
    (; t, dt, u, f, p, uprev) = integrator

    tdt = t + dt
    k = order
    reinitFBDF!(integrator, cache)

    # Predictor: evaluate Lagrange interpolant through u_history at Θ=1
    # using actual (variable) theta nodes. No need to fill integrator.k.
    n_pred = k + 1
    pred_thetas = Vector{typeof(t)}(undef, n_pred)
    u₀ = zero(u)
    if iters_from_event >= 1
        for j in 1:n_pred
            pred_thetas[j] = (ts[j] - t) / dt
        end
        u₀ = _eval_lagrange_oop(one(t), pred_thetas, u_history, n_pred)
    else
        u₀ = u
    end
    markfirststage!(nlsolver)

    nlsolver.z = u₀
    mass_matrix = f.mass_matrix

    # Corrector: evaluate Lagrange interpolant at equidistant past points
    if u isa Number
        fill!(u_corrector, zero(eltype(u)))
    else
        for i in eachindex(u_corrector)
            u_corrector[i] = zero(u_corrector[i])
        end
    end
    if iters_from_event >= 1
        for i in 1:(k - 1)
            u_corrector[i] = _eval_lagrange_oop(
                oftype(t, -i), pred_thetas, u_history, n_pred
            )
        end
    end
    if u isa Number
        tmp = -uprev * bdf_coeffs[k, 2]
        for i in 1:(k - 1)
            tmp -= u_corrector[i] * bdf_coeffs[k, i + 2]
        end
    else
        tmp = -uprev * bdf_coeffs[k, 2]
        for i in 1:(k - 1)
            tmp = @.. tmp - u_corrector[i] * bdf_coeffs[k, i + 2]
        end
    end

    if mass_matrix === I
        nlsolver.tmp = tmp / dt
    else
        nlsolver.tmp = mass_matrix * tmp / dt
    end
    β₀ = inv(bdf_coeffs[k, 1])
    α₀ = 1 #bdf_coeffs[k,1]
    nlsolver.γ = β₀
    nlsolver.α = α₀

    nlsolver.method = COEFFICIENT_MULTISTEP
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    u = z

    for j in 2:k
        r[j] = (1 - j)
        for i in 2:(k + 1)
            r[j] *= ((tdt - j * dt) - ts[i]) / (i * dt)
        end
    end

    terkp1 = (u - u₀)
    for j in 1:(k + 1)
        terkp1 *= j * dt / (tdt - ts[j])
    end

    lte = -1 / (1 + k)
    for j in 2:k
        lte -= bdf_coeffs[k, j] * r[j]
    end
    lte *= terkp1

    if integrator.opts.adaptive
        for i in 1:(k + 1)
            ts_tmp[i + 1] = ts[i]
        end
        ts_tmp[1] = tdt
        atmp = calculate_residuals(
            lte, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        terk = estimate_terk(integrator, cache, k + 1, Val(max_order), u)
        fd_weights = calc_finite_difference_weights(ts_tmp, tdt, k, Val(max_order))
        terk = @.. broadcast = false fd_weights[1, k + 1] * u
        if u isa Number
            for i in 2:(k + 1)
                terk += fd_weights[i, k + 1] * u_history[i - 1]
            end
            terk *= abs(dt^(k))
        else
            for i in 2:(k + 1)
                terk = @.. terk + fd_weights[i, k + 1] * u_history[i - 1]
            end
            terk *= abs(dt^(k))
        end

        atmp = calculate_residuals(
            terk, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        cache.terk = integrator.opts.internalnorm(atmp, t)

        if k > 1
            terkm1 = estimate_terk(integrator, cache, k, Val(max_order), u)
            atmp = calculate_residuals(
                terkm1, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            cache.terkm1 = integrator.opts.internalnorm(atmp, t)
        end
        if k > 2
            terkm2 = estimate_terk(integrator, cache, k - 1, Val(max_order), u)
            atmp = calculate_residuals(
                terkm2, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            cache.terkm2 = integrator.opts.internalnorm(atmp, t)
        end
        if cache.qwait == 0 && k < max_order
            atmp = calculate_residuals(
                terkp1, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            cache.terkp1 = integrator.opts.internalnorm(atmp, t)
        else
            cache.terkp1 = zero(cache.terk)
        end
    end

    integrator.fsallast = f(u, p, tdt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    if integrator.opts.calck
        # Store dense output: resample Lagrange interpolant at Chebyshev nodes
        n = min(k + 1, max_order + 1)
        calck_thetas = Vector{typeof(t)}(undef, n)
        calck_thetas[1] = one(t)
        for j in 1:min(k, max_order)
            calck_thetas[1 + j] = (ts[j] - t) / dt
        end
        calck_values = Vector{typeof(u)}(undef, n)
        calck_values[1] = u isa Number ? u : copy(u)
        for j in 1:min(k, max_order)
            calck_values[1 + j] = u isa Number ? u_history[j] : copy(u_history[j])
        end
        _resample_at_chebyshev!(integrator.k, calck_values, calck_thetas, n)
        for j in (n + 1):(max_order + 1)
            integrator.k[j] = zero(u)
        end
    end
    return integrator.u = u
end

function initialize!(integrator, cache::FBDFCache{max_order}) where {max_order}
    integrator.kshortsize = max_order + 1

    resize!(integrator.k, integrator.kshortsize)
    for i in 1:(max_order + 1)
        integrator.k[i] = cache.dense[i]
    end
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    u_modified = integrator.u_modified
    integrator.u_modified = true
    reinitFBDF!(integrator, cache)
    return integrator.u_modified = u_modified
end

function perform_step!(
        integrator, cache::FBDFCache{max_order},
        repeat_step = false
    ) where {max_order}
    (;
        ts, u_history, order, u_corrector, bdf_coeffs, r, nlsolver, terk_tmp,
        terkp1_tmp, atmp, tmp, u₀, ts_tmp, equi_ts, dense, step_limiter!,
    ) = cache
    (; t, dt, u, f, p, uprev) = integrator

    reinitFBDF!(integrator, cache)
    k = order
    tdt = t + dt

    # Predictor: evaluate Lagrange interpolant through u_history at Θ=1
    # using actual (variable) theta nodes. No need to fill integrator.k.
    n_pred = k + 1
    @.. broadcast = false u₀ = zero(u)
    if cache.iters_from_event >= 1
        for j in 1:n_pred
            equi_ts[j] = (ts[j] - t) / dt
        end
        _eval_lagrange_iip!(u₀, one(t), equi_ts, u_history, n_pred)
    else
        @.. broadcast = false u₀ = u
    end
    markfirststage!(nlsolver)
    @.. broadcast = false nlsolver.z = u₀
    mass_matrix = f.mass_matrix

    # Corrector: evaluate Lagrange interpolant at equidistant past points
    for h in u_corrector
        fill!(h, zero(eltype(h)))
    end
    if cache.iters_from_event >= 1
        for i in 1:(k - 1)
            _eval_lagrange_iip!(u_corrector[i], oftype(t, -i), equi_ts, u_history, n_pred)
        end
    end
    @.. broadcast = false tmp = -uprev * bdf_coeffs[k, 2]
    for i in 1:(k - 1)
        @.. broadcast = false tmp -= u_corrector[i] * bdf_coeffs[k, i + 2]
    end

    invdt = inv(dt)
    if mass_matrix === I
        @.. broadcast = false nlsolver.tmp = tmp * invdt
    else
        @.. broadcast = false terkp1_tmp = tmp * invdt
        mul!(nlsolver.tmp, mass_matrix, terkp1_tmp)
    end

    β₀ = inv(bdf_coeffs[k, 1])
    α₀ = 1
    nlsolver.γ = β₀
    nlsolver.α = α₀
    nlsolver.method = COEFFICIENT_MULTISTEP
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast = false u = z

    step_limiter!(u, integrator, p, t + dt)

    #This is to correct the interpolation error of error estimation.
    for j in 2:k
        r[j] = (1 - j)
        for i in 2:(k + 1)
            r[j] *= ((tdt - j * dt) - ts[i]) / (i * dt)
        end
    end

    #for terkp1, we could use corrector and predictor to make an estimation.
    @.. broadcast = false terkp1_tmp = (u - u₀)
    for j in 1:(k + 1)
        @.. broadcast = false terkp1_tmp *= j * dt / (tdt - ts[j])
    end

    lte = -1 / (1 + k)
    for j in 2:k
        lte -= bdf_coeffs[k, j] * r[j]
    end
    @.. broadcast = false terk_tmp = lte * terkp1_tmp
    if integrator.opts.adaptive
        (; abstol, reltol, internalnorm) = integrator.opts
        for i in 1:(k + 1)
            ts_tmp[i + 1] = ts[i]
        end
        ts_tmp[1] = tdt
        calculate_residuals!(
            atmp, terk_tmp, uprev, u, abstol, reltol,
            internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
        estimate_terk!(integrator, cache, k + 1, Val(max_order))
        calculate_residuals!(
            atmp, terk_tmp, uprev, u, abstol, reltol,
            internalnorm, t
        )
        cache.terk = integrator.opts.internalnorm(atmp, t)

        if k > 1
            estimate_terk!(integrator, cache, k, Val(max_order))
            calculate_residuals!(
                atmp, terk_tmp, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            cache.terkm1 = integrator.opts.internalnorm(atmp, t)
        end
        if k > 2
            estimate_terk!(integrator, cache, k - 1, Val(max_order))
            calculate_residuals!(
                atmp, terk_tmp, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            cache.terkm2 = integrator.opts.internalnorm(atmp, t)
        end
        if cache.qwait == 0 && k < max_order
            calculate_residuals!(
                atmp, terkp1_tmp, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            cache.terkp1 = integrator.opts.internalnorm(atmp, t)
        else
            cache.terkp1 = zero(cache.terkp1)
        end
    end

    # Compute fsallast algebraically from NL solver convergence condition:
    #   f(u) = α*invγdt*u - nlsolver.tmp           (mass_matrix === I)
    #   f(u) = α*invγdt*(M*u) - nlsolver.tmp       (mass_matrix !== I)
    # This avoids calling f() through FunctionWrapper dispatch which can allocate.
    invγdt = inv(nlsolver.γ * dt)
    if mass_matrix === I
        @.. integrator.fsallast =
            nlsolver.α * invγdt * u - nlsolver.tmp
    else
        mul!(terkp1_tmp, mass_matrix, u)
        @.. integrator.fsallast =
            nlsolver.α * invγdt * terkp1_tmp - nlsolver.tmp
    end
    if integrator.opts.calck
        # Store dense output: resample Lagrange interpolant at Chebyshev nodes.
        # Use _resample_at_chebyshev_direct_iip! to read from u and u_history
        # directly, avoiding scratch buffer type mismatches during AD.
        n = min(k + 1, max_order + 1)
        equi_ts[1] = one(eltype(equi_ts))
        for j in 1:min(k, max_order)
            equi_ts[1 + j] = (ts[j] - t) / dt
        end
        _resample_at_chebyshev_direct_iip!(integrator.k, u, u_history, equi_ts, n)
        for j in (n + 1):(max_order + 1)
            fill!(integrator.k[j], zero(eltype(u)))
        end
    end
    return nothing
end
