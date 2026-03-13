function initialize!(
        integrator, cache::Union{
            Rosenbrock23Cache,
            Rosenbrock32Cache,
        }
    )
    integrator.kshortsize = 2
    (; k₁, k₂, fsalfirst, fsallast) = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = k₁
    integrator.k[2] = k₂
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(
        integrator,
        cache::Union{
            Rosenbrock23ConstantCache,
            Rosenbrock32ConstantCache,
        }
    )
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = zero(integrator.fsalfirst)
    return integrator.k[2] = zero(integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock23Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, opts) = integrator
    (; k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W, tmp, uf, tf, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter!) = cache
    (; c₃₂, d) = cache.tab

    # Assignments
    sizeu = size(u)
    mass_matrix = integrator.f.mass_matrix

    # Precalculations
    dtγ = dt * d
    neginvdtγ = -inv(dtγ)
    dto2 = dt / 2
    dto6 = dt / 6

    if repeat_step
        f(integrator.fsalfirst, uprev, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    calc_rosenbrock_differentiation!(integrator, cache, dtγ, dtγ, repeat_step)

    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtγ)
        )
    else
        linres = dolinsolve(
            integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtγ)
        )
    end

    vecu = _vec(linres.u)
    veck₁ = _vec(k₁)

    @.. veck₁ = vecu * neginvdtγ
    integrator.stats.nsolve += 1

    @.. u = uprev + dto2 * k₁
    stage_limiter!(u, integrator, p, t + dto2)
    f(f₁, u, p, t + dto2)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        copyto!(tmp, k₁)
    else
        mul!(_vec(tmp), mass_matrix, _vec(k₁))
    end

    @.. linsolve_tmp = f₁ - tmp

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck₂ = _vec(k₂)

    @.. veck₂ = vecu * neginvdtγ + veck₁
    integrator.stats.nsolve += 1

    @.. u = uprev + dt * k₂
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        f(fsallast, u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        if mass_matrix === I
            @.. broadcast = false linsolve_tmp = fsallast - c₃₂ * (k₂ - f₁) -
                2(k₁ - fsalfirst) + dt * dT
        else
            @.. broadcast = false du2 = c₃₂ * k₂ + 2k₁
            mul!(_vec(du1), mass_matrix, _vec(du2))
            @.. broadcast = false linsolve_tmp = fsallast - du1 + c₃₂ * f₁ + 2fsalfirst +
                dt * dT
        end

        linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
        vecu = _vec(linres.u)
        veck3 = _vec(k₃)
        @.. veck3 = vecu * neginvdtγ

        integrator.stats.nsolve += 1

        if mass_matrix === I
            @.. broadcast = false tmp = dto6 * (k₁ - 2 * k₂ + k₃)
        else
            veck₁ = _vec(k₁)
            veck₂ = _vec(k₂)
            veck₃ = _vec(k₃)
            vectmp = _vec(tmp)
            @.. broadcast = false vectmp = ifelse(
                cache.algebraic_vars,
                false, dto6 * (veck₁ - 2 * veck₂ + veck₃)
            )
        end
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            algvar = reshape(cache.algebraic_vars, size(u))
            invatol = inv(integrator.opts.abstol)
            @.. atmp = ifelse(algvar, fsallast, false) * invatol
            integrator.EEst += integrator.opts.internalnorm(atmp, t)
        end
    end
    cache.linsolve = linres.cache
end

@muladd function perform_step!(integrator, cache::Rosenbrock32Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, opts) = integrator
    (; k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W, tmp, uf, tf, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter!) = cache
    (; c₃₂, d) = cache.tab

    # Assignments
    sizeu = size(u)
    mass_matrix = integrator.f.mass_matrix

    # Precalculations
    dtγ = dt * d
    neginvdtγ = -inv(dtγ)
    dto2 = dt / 2
    dto6 = dt / 6

    if repeat_step
        f(integrator.fsalfirst, uprev, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    calc_rosenbrock_differentiation!(integrator, cache, dtγ, dtγ, repeat_step)

    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtγ)
        )
    else
        linres = dolinsolve(
            integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtγ)
        )
    end

    vecu = _vec(linres.u)
    veck₁ = _vec(k₁)

    @.. veck₁ = vecu * neginvdtγ
    integrator.stats.nsolve += 1

    @.. broadcast = false u = uprev + dto2 * k₁
    stage_limiter!(u, integrator, p, t + dto2)
    f(f₁, u, p, t + dto2)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        tmp .= k₁
    else
        mul!(_vec(tmp), mass_matrix, _vec(k₁))
    end

    @.. broadcast = false linsolve_tmp = f₁ - tmp

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck₂ = _vec(k₂)

    @.. veck₂ = vecu * neginvdtγ + veck₁
    integrator.stats.nsolve += 1

    @.. tmp = uprev + dt * k₂
    stage_limiter!(u, integrator, p, t + dt)
    f(fsallast, tmp, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast = false linsolve_tmp = fsallast - c₃₂ * (k₂ - f₁) - 2(k₁ - fsalfirst) +
            dt * dT
    else
        @.. broadcast = false du2 = c₃₂ * k₂ + 2k₁
        mul!(_vec(du1), mass_matrix, _vec(du2))
        @.. broadcast = false linsolve_tmp = fsallast - du1 + c₃₂ * f₁ + 2fsalfirst + dt * dT
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck3 = _vec(k₃)

    @.. veck3 = vecu * neginvdtγ
    integrator.stats.nsolve += 1

    @.. broadcast = false u = uprev + dto6 * (k₁ + 4k₂ + k₃)

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        @.. broadcast = false tmp = dto6 * (k₁ - 2 * k₂ + k₃)
        calculate_residuals!(
            atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            invatol = inv(integrator.opts.abstol)
            @.. atmp = ifelse(cache.algebraic_vars, fsallast, false) * invatol
            integrator.EEst += integrator.opts.internalnorm(atmp, t)
        end
    end
    cache.linsolve = linres.cache
end

@muladd function perform_step!(
        integrator, cache::Rosenbrock23ConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; c₃₂, d, tf, uf) = cache

    # Precalculations
    dtγ = dt * d
    neginvdtγ = -inv(dtγ)
    dto2 = dt / 2
    dto6 = dt / 6

    if repeat_step
        integrator.fsalfirst = f(uprev, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtγ, repeat_step)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    k₁ = _reshape(W \ _vec((integrator.fsalfirst + dtγ * dT)), axes(uprev)) * neginvdtγ
    integrator.stats.nsolve += 1
    tmp = @.. uprev + dto2 * k₁
    f₁ = f(tmp, p, t + dto2)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        k₂ = _reshape(W \ _vec(f₁ - k₁), axes(uprev))
    else
        k₂ = _reshape(W \ _vec(f₁ - mass_matrix * k₁), axes(uprev))
    end
    k₂ = @.. k₂ * neginvdtγ + k₁
    integrator.stats.nsolve += 1
    u = uprev + dt * k₂

    if integrator.opts.adaptive
        integrator.fsallast = f(u, p, t + dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        if mass_matrix === I
            linsolve_tmp = @.. (
                integrator.fsallast - c₃₂ * (k₂ - f₁) -
                    2 * (k₁ - integrator.fsalfirst) + dt * dT
            )
        else
            linsolve_tmp = mass_matrix * (@.. c₃₂ * k₂ + 2 * k₁)
            linsolve_tmp = @.. (
                integrator.fsallast - linsolve_tmp +
                    c₃₂ * f₁ + 2 * integrator.fsalfirst + dt * dT
            )
        end
        k₃ = _reshape(W \ _vec(linsolve_tmp), axes(uprev)) * neginvdtγ
        integrator.stats.nsolve += 1

        if u isa Number
            utilde = dto6 * f.mass_matrix[1, 1] * (k₁ - 2 * k₂ + k₃)
        else
            utilde = f.mass_matrix * (@.. dto6 * (k₁ - 2 * k₂ + k₃))
        end
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            invatol = inv(integrator.opts.abstol)
            atmp = @. ifelse(integrator.differential_vars, false, integrator.fsallast) *
                invatol
            integrator.EEst += integrator.opts.internalnorm(atmp, t)
        end
    end
    integrator.k[1] = k₁
    integrator.k[2] = k₂
    integrator.u = u
    return nothing
end

@muladd function perform_step!(
        integrator, cache::Rosenbrock32ConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; c₃₂, d, tf, uf) = cache

    # Precalculations
    dtγ = dt * d
    neginvdtγ = -inv(dtγ)
    dto2 = dt / 2
    dto6 = dt / 6

    mass_matrix = integrator.f.mass_matrix

    if repeat_step
        integrator.fsalfirst = f(uprev, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    # Time derivative
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtγ, repeat_step)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    k₁ = _reshape(W \ -_vec((integrator.fsalfirst + dtγ * dT)), axes(uprev)) / dtγ
    integrator.stats.nsolve += 1
    tmp = @.. uprev + dto2 * k₁
    f₁ = f(tmp, p, t + dto2)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        k₂ = _reshape(W \ _vec(f₁ - k₁), axes(uprev))
    else
        linsolve_tmp = f₁ - mass_matrix * k₁
        k₂ = _reshape(W \ _vec(linsolve_tmp), axes(uprev))
    end
    k₂ = @.. k₂ * neginvdtγ + k₁

    integrator.stats.nsolve += 1
    tmp = @.. uprev + dt * k₂
    integrator.fsallast = f(tmp, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        linsolve_tmp = @.. (
            integrator.fsallast - c₃₂ * (k₂ - f₁) -
                2(k₁ - integrator.fsalfirst) + dt * dT
        )
    else
        linsolve_tmp = mass_matrix * (@.. c₃₂ * k₂ + 2 * k₁)
        linsolve_tmp = @.. (
            integrator.fsallast - linsolve_tmp +
                c₃₂ * f₁ + 2 * integrator.fsalfirst + dt * dT
        )
    end
    k₃ = _reshape(W \ _vec(linsolve_tmp), axes(uprev)) * neginvdtγ
    integrator.stats.nsolve += 1
    u = @.. uprev + dto6 * (k₁ + 4k₂ + k₃)

    if integrator.opts.adaptive
        utilde = @.. dto6 * (k₁ - 2k₂ + k₃)
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            invatol = inv(integrator.opts.abstol)
            atmp = ifelse(integrator.differential_vars, false, integrator.fsallast) .* invatol
            integrator.EEst += integrator.opts.internalnorm(atmp, t)
        end
    end

    integrator.k[1] = k₁
    integrator.k[2] = k₂
    integrator.u = u
    return nothing
end

function initialize!(
        integrator,
        cache::Union{
            Rosenbrock33ConstantCache,
            Rosenbrock34ConstantCache,
        }
    )
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function initialize!(
        integrator,
        cache::Union{
            Rosenbrock33Cache,
            Rosenbrock34Cache,
        }
    )
    integrator.kshortsize = 2
    (; fsalfirst, fsallast) = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = fsalfirst
    integrator.k[2] = fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(
        integrator, cache::Rosenbrock33ConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; tf, uf) = cache
    (; a21, a31, a32, C21, C31, C32, b1, b2, b3, btilde1, btilde2, btilde3, gamma, c2, c3, d1, d2, d3) = cache.tab

    # Precalculations
    dtC21 = C21 / dt
    dtC31 = C31 / dt
    dtC32 = C32 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtgamma = dt * gamma

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    linsolve_tmp = integrator.fsalfirst + dtd1 * dT

    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a21 * k1
    du = f(u, p, t + c2 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        linsolve_tmp = du + dtd2 * dT + dtC21 * k1
    else
        linsolve_tmp = du + dtd2 * dT + mass_matrix * (dtC21 * k1)
    end

    k2 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a31 * k1 + a32 * k2
    du = f(u, p, t + c3 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        linsolve_tmp = du + dtd3 * dT + dtC31 * k1 + dtC32 * k2
    else
        linsolve_tmp = du + dtd3 * dT + mass_matrix * (dtC31 * k1 + dtC32 * k2)
    end

    k3 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + b1 * k1 + b2 * k2 + b3 * k3
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        utilde = btilde1 * k1 + btilde2 * k2 + btilde3 * k3
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return nothing
end

@muladd function perform_step!(integrator, cache::Rosenbrock33Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; du, du1, du2, fsalfirst, fsallast, k1, k2, k3, dT, J, W, uf, tf, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter!) = cache
    (; a21, a31, a32, C21, C31, C32, b1, b2, b3, btilde1, btilde2, btilde3, gamma, c2, c3, d1, d2, d3) = cache.tab

    # Assignments
    mass_matrix = integrator.f.mass_matrix
    sizeu = size(u)
    utilde = du

    # Precalculations
    dtC21 = C21 / dt
    dtC31 = C31 / dt
    dtC32 = C32 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtgamma = dt * gamma

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step)

    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma)
        )
    else
        linres = dolinsolve(
            integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma)
        )
    end

    vecu = _vec(linres.u)
    veck1 = _vec(k1)

    @.. broadcast = false veck1 = -vecu
    integrator.stats.nsolve += 1

    @.. broadcast = false u = uprev + a21 * k1
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(du, u, p, t + c2 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast = false linsolve_tmp = du + dtd2 * dT + dtC21 * k1
    else
        @.. broadcast = false du1 = dtC21 * k1
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast = false linsolve_tmp = du + dtd2 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck2 = _vec(k2)

    @.. broadcast = false veck2 = -vecu

    integrator.stats.nsolve += 1

    @.. broadcast = false u = uprev + a31 * k1 + a32 * k2
    stage_limiter!(u, integrator, p, t + c3 * dt)
    f(du, u, p, t + c3 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast = false linsolve_tmp = du + dtd3 * dT + dtC31 * k1 + dtC32 * k2
    else
        @.. broadcast = false du1 = dtC31 * k1 + dtC32 * k2
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast = false linsolve_tmp = du + dtd3 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck3 = _vec(k3)

    @.. broadcast = false veck3 = -vecu

    integrator.stats.nsolve += 1

    @.. broadcast = false u = uprev + b1 * k1 + b2 * k2 + b3 * k3

    step_limiter!(u, integrator, p, t + dt)

    f(fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        @.. broadcast = false utilde = btilde1 * k1 + btilde2 * k2 + btilde3 * k3
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    cache.linsolve = linres.cache
end

################################################################################

@muladd function perform_step!(
        integrator, cache::Rosenbrock34ConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; tf, uf) = cache
    (; a21, a31, a32, a41, a42, a43, C21, C31, C32, C41, C42, C43, b1, b2, b3, b4, btilde1, btilde2, btilde3, btilde4, gamma, c2, c3, d1, d2, d3, d4) = cache.tab

    # Precalculations
    dtC21 = C21 / dt
    dtC31 = C31 / dt
    dtC32 = C32 / dt
    dtC41 = C41 / dt
    dtC42 = C42 / dt
    dtC43 = C43 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtd4 = dt * d4
    dtgamma = dt * gamma

    mass_matrix = integrator.f.mass_matrix
    # Time derivative
    tf.u = uprev
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    linsolve_tmp = integrator.fsalfirst + dtd1 * dT

    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev # +a21*k1 a21 == 0
    # du = f(u, p, t+c2*dt) c2 == 0 and a21 == 0 => du = f(uprev, p, t) == fsalfirst

    if mass_matrix === I
        linsolve_tmp = integrator.fsalfirst + dtd2 * dT + dtC21 * k1
    else
        linsolve_tmp = integrator.fsalfirst + dtd2 * dT + mass_matrix * (dtC21 * k1)
    end

    k2 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a31 * k1 + a32 * k2
    du = f(u, p, t + c3 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        linsolve_tmp = du + dtd3 * dT + dtC31 * k1 + dtC32 * k2
    else
        linsolve_tmp = du + dtd3 * dT + mass_matrix * (dtC31 * k1 + dtC32 * k2)
    end

    k3 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a41 * k1 + a42 * k2 + a43 * k3
    du = f(u, p, t + dt) #-- c4 = 1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        linsolve_tmp = du + dtd4 * dT + dtC41 * k1 + dtC42 * k2 + dtC43 * k3
    else
        linsolve_tmp = du + dtd4 * dT + mass_matrix * (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    end

    k4 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        utilde = btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return nothing
end

@muladd function perform_step!(integrator, cache::Rosenbrock34Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; du, du1, du2, fsalfirst, fsallast, k1, k2, k3, k4, dT, J, W, uf, tf, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter!) = cache
    (; a21, a31, a32, a41, a42, a43, C21, C31, C32, C41, C42, C43, b1, b2, b3, b4, btilde1, btilde2, btilde3, btilde4, gamma, c2, c3, d1, d2, d3, d4) = cache.tab

    # Assignments
    uidx = eachindex(integrator.uprev)
    sizeu = size(u)
    mass_matrix = integrator.f.mass_matrix
    utilde = du

    # Precalculations
    dtC21 = C21 / dt
    dtC31 = C31 / dt
    dtC32 = C32 / dt
    dtC41 = C41 / dt
    dtC42 = C42 / dt
    dtC43 = C43 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtd4 = dt * d4
    dtgamma = dt * gamma

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step)

    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma)
        )
    else
        linres = dolinsolve(
            integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma)
        )
    end

    vecu = _vec(linres.u)
    veck1 = _vec(k1)

    @.. broadcast = false veck1 = -vecu
    integrator.stats.nsolve += 1

    #=
    a21 == 0 and c2 == 0
    so du = integrator.fsalfirst!
    @.. broadcast=false u = uprev + a21*k1

    f(du, u, p, t+c2*dt)
    =#

    if mass_matrix === I
        @.. broadcast = false linsolve_tmp = fsalfirst + dtd2 * dT + dtC21 * k1
    else
        @.. broadcast = false du1 = dtC21 * k1
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast = false linsolve_tmp = fsalfirst + dtd2 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck2 = _vec(k2)
    @.. broadcast = false veck2 = -vecu
    integrator.stats.nsolve += 1

    @.. broadcast = false u = uprev + a31 * k1 + a32 * k2
    stage_limiter!(u, integrator, p, t + c3 * dt)
    f(du, u, p, t + c3 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast = false linsolve_tmp = du + dtd3 * dT + dtC31 * k1 + dtC32 * k2
    else
        @.. broadcast = false du1 = dtC31 * k1 + dtC32 * k2
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast = false linsolve_tmp = du + dtd3 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck3 = _vec(k3)
    @.. broadcast = false veck3 = -vecu
    integrator.stats.nsolve += 1
    @.. broadcast = false u = uprev + a41 * k1 + a42 * k2 + a43 * k3
    stage_limiter!(u, integrator, p, t + dt)
    f(du, u, p, t + dt) #-- c4 = 1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast = false linsolve_tmp = du + dtd4 * dT + dtC41 * k1 + dtC42 * k2 +
            dtC43 * k3
    else
        @.. broadcast = false du1 = dtC41 * k1 + dtC42 * k2 + dtC43 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast = false linsolve_tmp = du + dtd4 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck4 = _vec(k4)
    @.. broadcast = false veck4 = -vecu
    integrator.stats.nsolve += 1

    @.. broadcast = false u = uprev + b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4

    step_limiter!(u, integrator, p, t + dt)

    f(fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        @.. broadcast = false utilde = btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    cache.linsolve = linres.cache
end

################################################################################

#### Shared stage loop for all Rosenbrock methods (generic + Rodas)

@muladd function _rosenbrock_stage_loop!(integrator, cache, tab_a, dtC, tab_c, dtd, ks,
        mass_matrix, linres, _stage_limiter!)
    (; t, dt, uprev, u, f, p) = integrator
    (; du, du1, du2, dT, linsolve_tmp) = cache
    s = length(ks)
    for stage in 2:s
        u .= uprev
        for j in 1:(stage - 1)
            @.. u += tab_a[stage, j] * ks[j]
        end
        _stage_limiter!(u, integrator, p, t + tab_c[stage] * dt)
        f(du, u, p, t + tab_c[stage] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        du1 .= 0
        if mass_matrix === I
            for j in 1:(stage - 1)
                @.. du1 += dtC[stage, j] * ks[j]
            end
        else
            for j in 1:(stage - 1)
                @.. du1 += dtC[stage, j] * ks[j]
            end
            mul!(_vec(du2), mass_matrix, _vec(du1))
            du1 .= du2
        end
        @.. linsolve_tmp = du + dtd[stage] * dT + du1
        linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
        @.. $(_vec(ks[stage])) = -linres.u
        integrator.stats.nsolve += 1
    end
    return linres
end

################################################################################

#### Generic Rosenbrock initialize! (replaces generated code from macros)

function initialize!(integrator, cache::GenericRosenbrockConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::GenericRosenbrockMutableCache)
    integrator.kshortsize = 2
    (; fsalfirst, fsallast) = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = fsallast
    resize!(integrator.k, integrator.kshortsize)
    integrator.k .= [fsalfirst, fsallast]
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

################################################################################

#### Generic Rosenbrock perform_step! — inplace (mutable cache)

@muladd function perform_step!(integrator, cache::GenericRosenbrockMutableCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; du, du1, du2, fsallast, dT, J, W, uf, tf, linsolve_tmp, jac_config, atmp, weight) = cache
    tab = cache.tab
    ks = _ks(cache)
    s = length(ks)
    mass_matrix = integrator.f.mass_matrix

    # Precalculations — compute dtC and dtd from the full-array tableau
    dtC = tab.C ./ dt
    dtd = dt .* tab.d
    dtgamma = dt * tab.gamma

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)

    calc_rosenbrock_differentiation!(integrator, cache, dtd[1], dtgamma, repeat_step)

    linsolve = cache.linsolve

    # Stage 1 solve
    linres = dolinsolve(
        integrator, linsolve; A = !repeat_step ? W : nothing, b = _vec(linsolve_tmp))
    @.. $(_vec(ks[1])) = -linres.u
    integrator.stats.nsolve += 1

    # Stages 2..s (shared loop)
    linres = _rosenbrock_stage_loop!(integrator, cache, tab.a, dtC, tab.c, dtd, ks,
        mass_matrix, linres, Returns(nothing))

    # Final solution: u = uprev + sum(b[i]*k[i])
    u .= uprev
    for i in 1:s
        @.. u += tab.b[i] * ks[i]
    end

    f(fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Adaptive error estimation
    if tab isa RosenbrockAdaptiveTableau && integrator.opts.adaptive
        utilde = du
        @.. utilde = zero(eltype(u))
        for i in 1:s
            if !iszero(tab.btilde[i])
                @.. utilde += tab.btilde[i] * ks[i]
            end
        end
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    cache.linsolve = linres.cache
end

################################################################################

#### Generic Rosenbrock perform_step! — non-inplace (constant cache)

@muladd function perform_step!(integrator, cache::GenericRosenbrockConstantCache,
        repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tf, uf) = cache
    tab = cache.tab

    # Precalculations
    dtC = tab.C ./ dt
    dtd = dt .* tab.d
    dtgamma = dt * tab.gamma

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    tf.u = uprev
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step)
    linsolve_tmp = integrator.fsalfirst + dtd[1] * dT

    # Stage 1
    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1

    s = length(tab.b)
    ks = ntuple(Returns(k1), Val(6))  # max 6 stages for generic methods

    # Stages 2..s
    for stage in 2:s
        u = uprev
        for j in 1:(stage - 1)
            u = @.. u + tab.a[stage, j] * ks[j]
        end

        du = f(u, p, t + tab.c[stage] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        linsolve_tmp = zero(du)
        if mass_matrix === I
            for j in 1:(stage - 1)
                linsolve_tmp = @.. linsolve_tmp + dtC[stage, j] * ks[j]
            end
        else
            for j in 1:(stage - 1)
                linsolve_tmp = @.. linsolve_tmp + dtC[stage, j] * ks[j]
            end
            linsolve_tmp = mass_matrix * linsolve_tmp
        end
        linsolve_tmp = @.. du + dtd[stage] * dT + linsolve_tmp

        ks = Base.setindex(ks, _reshape(W \ -_vec(linsolve_tmp), axes(uprev)), stage)
        integrator.stats.nsolve += 1
    end

    # Final solution: u = uprev + sum(b[i]*k[i])
    u = uprev
    for i in 1:s
        u = @.. u + tab.b[i] * ks[i]
    end

    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u

    # Adaptive error estimation
    if tab isa RosenbrockAdaptiveTableau && integrator.opts.adaptive
        utilde = zero(u)
        for i in 1:s
            if !iszero(tab.btilde[i])
                utilde = @.. utilde + tab.btilde[i] * ks[i]
            end
        end
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
end

################################################################################

#### Rodas3P type method

function initialize!(integrator, cache::Union{Rodas23WConstantCache, Rodas3PConstantCache})
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    # Avoid undefined entries if k is an array of arrays
    integrator.k[1] = zero(integrator.u)
    integrator.k[2] = zero(integrator.u)
    return integrator.k[3] = zero(integrator.u)
end

@muladd function perform_step!(
        integrator, cache::Union{Rodas23WConstantCache, Rodas3PConstantCache},
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; tf, uf) = cache
    (; a21, a41, a42, a43, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, gamma, c2, c3, d1, d2, d3) = cache.tab

    # Precalculations
    dtC21 = C21 / dt
    dtC31 = C31 / dt
    dtC32 = C32 / dt
    dtC41 = C41 / dt
    dtC42 = C42 / dt
    dtC43 = C43 / dt
    dtC51 = C51 / dt
    dtC52 = C52 / dt
    dtC53 = C53 / dt
    dtC54 = C54 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtgamma = dt * gamma

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    tf.u = uprev
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    du = f(uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    k3 = copy(du)  #-- save for stage 3

    linsolve_tmp = du + dtd1 * dT

    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a21 * k1
    du = f(u, p, t + c2 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        linsolve_tmp = du + dtd2 * dT + dtC21 * k1
    else
        linsolve_tmp = du + dtd2 * dT + mass_matrix * (dtC21 * k1)
    end

    k2 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1

    if mass_matrix === I
        linsolve_tmp = k3 + dtd3 * dT + (dtC31 * k1 + dtC32 * k2)
    else
        linsolve_tmp = k3 + dtd3 * dT + mass_matrix * (dtC31 * k1 + dtC32 * k2)
    end

    k3 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a41 * k1 + a42 * k2 + a43 * k3
    du = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        linsolve_tmp = du + (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    else
        linsolve_tmp = du + mass_matrix * (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    end

    k4 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1

    if mass_matrix === I
        linsolve_tmp = du + (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    else
        linsolve_tmp = du +
            mass_matrix * (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    end

    k5 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    du = u + k4 #-- solution p=2
    u = u + k5 #-- solution p=3

    EEst = 0.0
    if integrator.opts.calck
        (; h21, h22, h23, h24, h25, h31, h32, h33, h34, h35, h2_21, h2_22, h2_23, h2_24, h2_25) = cache.tab
        integrator.k[1] = h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 + h25 * k5
        integrator.k[2] = h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 + h35 * k5
        integrator.k[3] = h2_21 * k1 + h2_22 * k2 + h2_23 * k3 + h2_24 * k4 + h2_25 * k5
        if integrator.opts.adaptive
            if isa(linsolve_tmp, AbstractFloat)
                u_int, u_diff = calculate_interpoldiff(
                    uprev, du, u, integrator.k[1], integrator.k[2], integrator.k[3]
                )
            else
                u_int = linsolve_tmp
                u_diff = linsolve_tmp .+ 0
                calculate_interpoldiff!(
                    u_int, u_diff, uprev, du, u, integrator.k[1],
                    integrator.k[2], integrator.k[3]
                )
            end
            atmp = calculate_residuals(
                u_diff, uprev, u_int, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            EEst = max(EEst, integrator.opts.internalnorm(atmp, t))  #-- role of t unclear
        end
    end

    if (integrator.alg isa Rodas23W)
        k1 = u .+ 0
        u = du .+ 0
        du = k1 .+ 0
        if integrator.opts.calck
            integrator.k[1] = integrator.k[3] .+ 0
            integrator.k[2] = 0 * integrator.k[2]
        end
    end

    if integrator.opts.adaptive
        atmp = calculate_residuals(
            u - du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = max(EEst, integrator.opts.internalnorm(atmp, t))
    end

    integrator.u = u
    return nothing
end

function initialize!(integrator, cache::Union{Rodas23WCache, Rodas3PCache})
    integrator.kshortsize = 3
    (; dense1, dense2, dense3) = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = dense1
    integrator.k[2] = dense2
    return integrator.k[3] = dense3
end

@muladd function perform_step!(
        integrator, cache::Union{Rodas23WCache, Rodas3PCache}, repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; du, du1, du2, dT, J, W, uf, tf, k1, k2, k3, k4, k5, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter!) = cache
    (; a21, a41, a42, a43, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, gamma, c2, c3, d1, d2, d3) = cache.tab

    # Assignments
    sizeu = size(u)
    uidx = eachindex(integrator.uprev)
    mass_matrix = integrator.f.mass_matrix

    # Precalculations
    dtC21 = C21 / dt
    dtC31 = C31 / dt
    dtC32 = C32 / dt
    dtC41 = C41 / dt
    dtC42 = C42 / dt
    dtC43 = C43 / dt
    dtC51 = C51 / dt
    dtC52 = C52 / dt
    dtC53 = C53 / dt
    dtC54 = C54 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtgamma = dt * gamma

    f(cache.fsalfirst, uprev, p, t) # used in calc_rosenbrock_differentiation!
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step)

    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = cache.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma)
        )
    else
        linres = dolinsolve(
            integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = cache.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma)
        )
    end

    @.. broadcast = false $(_vec(k1)) = -linres.u

    integrator.stats.nsolve += 1

    @.. broadcast = false u = uprev + a21 * k1
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(du, u, p, t + c2 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast = false linsolve_tmp = du + dtd2 * dT + dtC21 * k1
    else
        @.. broadcast = false du1 = dtC21 * k1
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast = false linsolve_tmp = du + dtd2 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast = false $(_vec(k2)) = -linres.u
    integrator.stats.nsolve += 1

    if mass_matrix === I
        @.. broadcast = false linsolve_tmp = cache.fsalfirst + dtd3 * dT +
            (dtC31 * k1 + dtC32 * k2)
    else
        @.. broadcast = false du1 = dtC31 * k1 + dtC32 * k2
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast = false linsolve_tmp = cache.fsalfirst + dtd3 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast = false $(_vec(k3)) = -linres.u
    integrator.stats.nsolve += 1

    @.. broadcast = false u = uprev + a41 * k1 + a42 * k2 + a43 * k3
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(du, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast = false linsolve_tmp = du +
            (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    else
        @.. broadcast = false du1 = dtC41 * k1 + dtC42 * k2 + dtC43 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast = false linsolve_tmp = du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast = false $(_vec(k4)) = -linres.u
    integrator.stats.nsolve += 1

    if mass_matrix === I
        @.. broadcast = false linsolve_tmp = du +
            (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    else
        @.. broadcast = false du1 = dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast = false linsolve_tmp = du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast = false $(_vec(k5)) = -linres.u
    integrator.stats.nsolve += 1

    du = u + k4 #-- p=2 solution
    u .+= k5

    step_limiter!(u, integrator, p, t + dt)

    EEst = 0.0
    if integrator.opts.calck
        (; h21, h22, h23, h24, h25, h31, h32, h33, h34, h35, h2_21, h2_22, h2_23, h2_24, h2_25) = cache.tab
        @.. broadcast = false integrator.k[1] = h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 +
            h25 * k5
        @.. broadcast = false integrator.k[2] = h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 +
            h35 * k5
        @.. broadcast = false integrator.k[3] = h2_21 * k1 + h2_22 * k2 + h2_23 * k3 +
            h2_24 * k4 + h2_25 * k5
        if integrator.opts.adaptive
            calculate_interpoldiff!(
                du1, du2, uprev, du, u, integrator.k[1], integrator.k[2], integrator.k[3]
            )
            calculate_residuals!(
                atmp, du2, uprev, du1, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            EEst = max(EEst, integrator.opts.internalnorm(atmp, t))  #-- role of t unclear
        end
    end

    if (integrator.alg isa Rodas23W)
        du1[:] = u[:]
        u[:] = du[:]
        du[:] = du1[:]
        if integrator.opts.calck
            integrator.k[1][:] = integrator.k[3][:]
            integrator.k[2][:] .= 0.0
        end
    end

    if integrator.opts.adaptive
        calculate_residuals!(
            atmp, u - du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = max(EEst, integrator.opts.internalnorm(atmp, t))
    end
    cache.linsolve = linres.cache
end

function calculate_interpoldiff(uprev, up2, up3, c_koeff, d_koeff, c2_koeff)
    u_int = 0.0
    u_diff = 0.0
    a1 = up3 + c_koeff - up2 - c2_koeff
    a2 = d_koeff - c_koeff + c2_koeff
    a3 = -d_koeff
    dis = a2^2 - 3 * a1 * a3
    u_int = up3
    u_diff = 0.0
    if dis > 0.0 #-- Min/Max occurs
        tau1 = (-a2 - sqrt(dis)) / (3 * a3)
        tau2 = (-a2 + sqrt(dis)) / (3 * a3)
        if tau1 > tau2
            tau1, tau2 = tau2, tau1
        end
        for tau in (tau1, tau2)
            if (tau > 0.0) && (tau < 1.0)
                y_tau = (1 - tau) * uprev +
                    tau * (up3 + (1 - tau) * (c_koeff + tau * d_koeff))
                dy_tau = ((a3 * tau + a2) * tau + a1) * tau
                if abs(dy_tau) > abs(u_diff)
                    u_diff = dy_tau
                    u_int = y_tau
                end
            end
        end
    end
    return u_int, u_diff
end

function calculate_interpoldiff!(u_int, u_diff, uprev, up2, up3, c_koeff, d_koeff, c2_koeff)
    for i in eachindex(up2)
        a1 = up3[i] + c_koeff[i] - up2[i] - c2_koeff[i]
        a2 = d_koeff[i] - c_koeff[i] + c2_koeff[i]
        a3 = -d_koeff[i]
        dis = a2^2 - 3 * a1 * a3
        u_int[i] = up3[i]
        u_diff[i] = 0.0
        if dis > 0.0 #-- Min/Max occurs
            tau1 = (-a2 - sqrt(dis)) / (3 * a3)
            tau2 = (-a2 + sqrt(dis)) / (3 * a3)
            if tau1 > tau2
                tau1, tau2 = tau2, tau1
            end
            for tau in (tau1, tau2)
                if (tau > 0.0) && (tau < 1.0)
                    y_tau = (1 - tau) * uprev[i] +
                        tau * (up3[i] + (1 - tau) * (c_koeff[i] + tau * d_koeff[i]))
                    dy_tau = ((a3 * tau + a2) * tau + a1) * tau
                    if abs(dy_tau) > abs(u_diff[i])
                        u_diff[i] = dy_tau
                        u_int[i] = y_tau
                    end
                end
            end
        end
    end
    return
end

#### Rodas4 type method

function initialize!(integrator, cache::RosenbrockCombinedConstantCache)
    integrator.kshortsize = size(cache.tab.H, 1)
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    # Avoid undefined entries if k is an array of arrays
    for i in 1:(integrator.kshortsize)
        integrator.k[i] = zero(integrator.u)
    end
    return
end

@muladd function perform_step!(
        integrator, cache::RosenbrockCombinedConstantCache, repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; tf, uf) = cache
    (; A, C, gamma, c, d, H) = cache.tab

    # Precalculations
    dtC = C ./ dt
    dtd = dt .* d
    dtgamma = dt * gamma

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    tf.u = uprev
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    # Initialize ks
    num_stages = size(A, 1)
    du = f(uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    linsolve_tmp = @.. du + dtd[1] * dT
    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    # constant number for type stability make sure this is greater than num_stages
    ks = ntuple(Returns(k1), Val(20))
    # Loop for stages
    for stage in 2:num_stages
        u = uprev
        for i in 1:(stage - 1)
            u = @.. u + A[stage, i] * ks[i]
        end

        du = f(u, p, t + c[stage] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        # Compute linsolve_tmp for current stage
        linsolve_tmp = zero(du)
        if mass_matrix === I
            for i in 1:(stage - 1)
                linsolve_tmp = @.. linsolve_tmp + dtC[stage, i] * ks[i]
            end
        else
            for i in 1:(stage - 1)
                linsolve_tmp = @.. linsolve_tmp + dtC[stage, i] * ks[i]
            end
            linsolve_tmp = mass_matrix * linsolve_tmp
        end
        linsolve_tmp = @.. du + dtd[stage] * dT + linsolve_tmp

        ks = Base.setindex(ks, _reshape(W \ -_vec(linsolve_tmp), axes(uprev)), stage)
        integrator.stats.nsolve += 1
    end
    if (integrator.alg isa Rodas6P)
        du = ks[16]
        u = uprev
        for i in 1:15
            u = @.. u + A[16, i] * ks[i]
        end
        u = u .+ ks[16]
    else
        du = ks[num_stages]
        u = u .+ ks[num_stages]
    end

    if integrator.opts.adaptive
        if (integrator.alg isa Rodas5Pe)
            du = 0.2606326497975715 * ks[1] - 0.005158627295444251 * ks[2] +
                1.3038988631109731 * ks[3] + 1.235000722062074 * ks[4] +
                -0.7931985603795049 * ks[5] - 1.005448461135913 * ks[6] -
                0.18044626132120234 * ks[7] + 0.17051519239113755 * ks[8]
        end
        atmp = calculate_residuals(
            du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        for j in eachindex(integrator.k)
            integrator.k[j] = zero(integrator.k[1])
        end
        for i in 1:num_stages
            for j in eachindex(integrator.k)
                integrator.k[j] = @.. integrator.k[j] + H[j, i] * ks[i]
            end
        end
        if (integrator.alg isa Rodas5Pr) && integrator.opts.adaptive &&
                (integrator.EEst < 1.0)
            k2 = 0.5 * (
                uprev + u +
                    0.5 * (integrator.k[1] + 0.5 * (integrator.k[2] + 0.5 * integrator.k[3]))
            )
            du1 = (0.25 * (integrator.k[2] + integrator.k[3]) - uprev + u) / dt
            du = f(k2, p, t + dt / 2)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            if mass_matrix === I
                du2 = du1 - du
            else
                du2 = mass_matrix * du1 - du
            end
            EEst = norm(du2) / norm(integrator.opts.abstol .+ integrator.opts.reltol .* k2)
            integrator.EEst = max(EEst, integrator.EEst)
        end
    end

    integrator.u = u
    return nothing
end

function initialize!(integrator, cache::RosenbrockCache)
    integrator.kshortsize = size(cache.tab.H, 1)
    resize!(integrator.k, integrator.kshortsize)
    for i in 1:(integrator.kshortsize)
        integrator.k[i] = cache.dense[i]
    end
    return
end

@muladd function perform_step!(integrator, cache::RosenbrockCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; du, du1, du2, dT, dtC, dtd, J, W, uf, tf, ks, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter!) = cache
    (; A, C, gamma, c, d, H) = cache.tab

    # Assignments
    sizeu = size(u)
    uidx = eachindex(integrator.uprev)
    mass_matrix = integrator.f.mass_matrix

    # Precalculations
    @. dtC = C * inv(dt)
    @. dtd = dt * d
    dtgamma = dt * gamma

    f(cache.fsalfirst, uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    calc_rosenbrock_differentiation!(integrator, cache, dtd[1], dtgamma, repeat_step)

    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = cache.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma)
        )
    else
        linres = dolinsolve(
            integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = cache.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma)
        )
    end

    @.. $(_vec(ks[1])) = -linres.u
    integrator.stats.nsolve += 1

    linres = _rosenbrock_stage_loop!(integrator, cache, A, dtC, c, dtd, ks,
        mass_matrix, linres, stage_limiter!)
    if (integrator.alg isa Rodas6P)
        du .= ks[16]
        u .= uprev
        for i in 1:15
            @.. u += A[16, i] * ks[i]
        end
        u .+= ks[16]
    else
        du .= ks[end]
        u .+= ks[end]
    end

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        if (integrator.alg isa Rodas5Pe)
            @.. du = 0.2606326497975715 * ks[1] - 0.005158627295444251 * ks[2] +
                1.3038988631109731 * ks[3] + 1.235000722062074 * ks[4] +
                -0.7931985603795049 * ks[5] - 1.005448461135913 * ks[6] -
                0.18044626132120234 * ks[7] + 0.17051519239113755 * ks[8]
        end
        calculate_residuals!(
            atmp, du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        for j in eachindex(integrator.k)
            integrator.k[j] .= 0
        end
        for i in eachindex(ks)
            for j in eachindex(integrator.k)
                @.. integrator.k[j] += H[j, i] * ks[i]
            end
        end
        if (integrator.alg isa Rodas5Pr) && integrator.opts.adaptive &&
                (integrator.EEst < 1.0)
            ks[2] = 0.5 * (
                uprev + u +
                    0.5 *
                    (integrator.k[1] + 0.5 * (integrator.k[2] + 0.5 * integrator.k[3]))
            )
            du1 = (0.25 * (integrator.k[2] + integrator.k[3]) - uprev + u) / dt
            f(du, ks[2], p, t + dt / 2)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            if mass_matrix === I
                @.. du2 = du1 - du
            else
                mul!(_vec(du2), mass_matrix, _vec(du1))
                @.. du2 -= du
            end
            EEst = norm(du2) /
                norm(integrator.opts.abstol .+ integrator.opts.reltol .* ks[2])
            integrator.EEst = max(EEst, integrator.EEst)
        end
    end
    cache.linsolve = linres.cache
end

################################################################################
# Tsit5DA - hybrid explicit/linear-implicit method for DAEs
################################################################################

function initialize!(integrator, cache::HybridExplicitImplicitConstantCache)
    integrator.kshortsize = size(cache.tab.H, 1)
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    for i in 1:(integrator.kshortsize)
        integrator.k[i] = zero(integrator.u)
    end
    return
end

@muladd function perform_step!(
        integrator, cache::HybridExplicitImplicitConstantCache, repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    (; tf, uf, tab) = cache
    (; A, C, gamma, b, bhat, c, d, H) = tab

    mass_matrix = integrator.f.mass_matrix
    num_stages = size(A, 1)

    # Detect algebraic variables
    if mass_matrix === I
        has_alg = false
        alg_vars = Int[]
        diff_vars = Int[]
    else
        n = length(uprev)
        diff_vars = findall(i -> mass_matrix[i, i] != 0, 1:n)
        alg_vars = findall(i -> mass_matrix[i, i] == 0, 1:n)
        has_alg = !isempty(alg_vars)
    end

    # Compute Jacobian and time derivative for DAE case
    g_z = g_y = W_z_factored = dT = nothing
    if has_alg
        tf.u = uprev
        dT = calc_tderivative(integrator, cache)

        uf.t = t
        autodiff_alg = cache.autodiff
        if autodiff_alg isa AutoFiniteDiff
            autodiff_alg = SciMLBase.@set autodiff_alg.dir = sign(dt)
        end
        J = DI.jacobian(uf, autodiff_alg, uprev)

        n_g = length(alg_vars)
        n_f = length(diff_vars)
        g_z = J[alg_vars, alg_vars]
        g_y = J[alg_vars, diff_vars]
        W_z_factored = lu(-gamma * g_z)
    end

    # Stage loop
    du = f(uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    ks1 = zero(uprev)
    if mass_matrix === I
        ks1 = dt .* du
    else
        for iv in diff_vars
            ks1[iv] = dt * du[iv]
        end
        if has_alg
            # rhs = g(U_1) + g_y * γ * l_1 + h * d_1 * g_t
            rhs_z = du[alg_vars] .+
                g_y * (C[1, 1] .* ks1[diff_vars]) .+
                dt * d[1] .* dT[alg_vars]
            ks1_alg = W_z_factored \ rhs_z
            for (idx, iv) in enumerate(alg_vars)
                ks1[iv] = ks1_alg[idx]
            end
        end
    end
    ks = ntuple(Returns(ks1), Val(12))

    for stage in 2:num_stages
        # Assemble stage value
        u_stage = copy(uprev)
        for j in 1:(stage - 1)
            u_stage = @.. u_stage + A[stage, j] * ks[j]
        end

        du = f(u_stage, p, t + c[stage] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        ks_new = zero(uprev)
        if mass_matrix === I
            ks_new = dt .* du
        else
            for iv in diff_vars
                ks_new[iv] = dt * du[iv]
            end
            if has_alg
                # Build coupling: g_y * (sum_{j<i} C[i,j]*ks[j][diff] + C[i,i]*ks_new[diff])
                gy_coupling = C[stage, stage] .* ks_new[diff_vars]
                for j in 1:(stage - 1)
                    gy_coupling = @.. gy_coupling + C[stage, j] * ks[j][diff_vars]
                end
                gy_term = g_y * gy_coupling

                # Build coupling: g_z * sum_{j<i} C[i,j]*ks[j][alg]
                gz_coupling = zeros(eltype(uprev), length(alg_vars))
                for j in 1:(stage - 1)
                    gz_coupling = @.. gz_coupling + C[stage, j] * ks[j][alg_vars]
                end
                gz_term = g_z * gz_coupling

                rhs_z = du[alg_vars] .+ gy_term .+ gz_term .+ dt * d[stage] .* dT[alg_vars]
                ks_alg = W_z_factored \ rhs_z
                for (idx, iv) in enumerate(alg_vars)
                    ks_new[iv] = ks_alg[idx]
                end
            end
        end
        ks = Base.setindex(ks, ks_new, stage)
    end

    # Solution update
    u = copy(uprev)
    for i in 1:num_stages
        u = @.. u + b[i] * ks[i]
    end

    # Error estimation
    if integrator.opts.adaptive
        err_vec = zero(uprev)
        for i in 1:num_stages
            err_vec = @.. err_vec + (b[i] - bhat[i]) * ks[i]
        end
        atmp = calculate_residuals(
            err_vec, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    # Dense output
    if integrator.opts.calck
        for j in eachindex(integrator.k)
            integrator.k[j] = zero(integrator.k[1])
        end
        for i in 1:num_stages
            for j in eachindex(integrator.k)
                integrator.k[j] = @.. integrator.k[j] + H[j, i] * ks[i]
            end
        end
    end

    integrator.u = u
    return nothing
end

function initialize!(integrator, cache::HybridExplicitImplicitCache)
    integrator.kshortsize = size(cache.tab.H, 1)
    resize!(integrator.k, integrator.kshortsize)
    for i in 1:(integrator.kshortsize)
        integrator.k[i] = cache.dense[i]
    end
    return
end

@muladd function perform_step!(integrator, cache::HybridExplicitImplicitCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (;
        du, du1, du2, dT, J, W, uf, tf, ks, linsolve_tmp, jac_config, atmp, weight,
        stage_limiter!, step_limiter!, diff_vars, alg_vars,
        g_z, g_y, linsolve_tmp_z,
    ) = cache
    W_z = cache.W_z
    (; A, C, gamma, b, bhat, c, d, H) = cache.tab

    mass_matrix = integrator.f.mass_matrix
    num_stages = size(A, 1)
    has_alg = !isempty(alg_vars)
    n_g = length(alg_vars)
    n_f = length(diff_vars)

    # Compute Jacobian for DAE case
    if has_alg && !repeat_step
        # Use existing Rosenbrock differentiation infrastructure
        dtgamma = dt * gamma
        calc_rosenbrock_differentiation!(integrator, cache, dt * d[1], dtgamma, repeat_step)

        # Extract algebraic Jacobian blocks from full J
        for (gi, ai) in enumerate(alg_vars)
            for (gj, aj) in enumerate(alg_vars)
                g_z[gi, gj] = J[ai, aj]
            end
            for (fj, dj) in enumerate(diff_vars)
                g_y[gi, fj] = J[ai, dj]
            end
        end

        # Form W_z = -gamma * g_z
        for gi in 1:n_g
            for gj in 1:n_g
                W_z[gi, gj] = -gamma * g_z[gi, gj]
            end
        end
    end

    # Evaluate f at initial point
    f(cache.fsalfirst, uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Stage 1
    for iv in eachindex(ks[1])
        ks[1][iv] = zero(eltype(u))
    end
    if mass_matrix === I
        @.. ks[1] = dt * cache.fsalfirst
    else
        for iv in diff_vars
            ks[1][iv] = dt * cache.fsalfirst[iv]
        end
        if has_alg
            # rhs = g(U_1) + g_y * γ * l_1 + h * d_1 * g_t
            for gi in 1:n_g
                val = cache.fsalfirst[alg_vars[gi]] + dt * d[1] * dT[alg_vars[gi]]
                for fj in 1:n_f
                    val += g_y[gi, fj] * C[1, 1] * ks[1][diff_vars[fj]]
                end
                linsolve_tmp_z[gi] = val
            end
            # Solve W_z * k_alg = rhs  (W_z = -γ*g_z)
            sol_z = W_z \ linsolve_tmp_z
            for (gi, ai) in enumerate(alg_vars)
                ks[1][ai] = sol_z[gi]
            end
        end
    end

    # Stages 2 through num_stages
    for stage in 2:num_stages
        # Assemble u_stage
        u .= uprev
        for j in 1:(stage - 1)
            @.. u += A[stage, j] * ks[j]
        end

        stage_limiter!(u, integrator, p, t + c[stage] * dt)
        f(du, u, p, t + c[stage] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        # Set ks[stage] differential part
        for iv in eachindex(ks[stage])
            ks[stage][iv] = zero(eltype(u))
        end
        if mass_matrix === I
            @.. ks[stage] = dt * du
        else
            for iv in diff_vars
                ks[stage][iv] = dt * du[iv]
            end
            if has_alg
                # rhs = g(U_i) + g_z*Σ_{j<i} γ_{ij}*k_j + g_y*Σ_{j≤i} γ_{ij}*l_j + h*d_i*g_t
                for gi in 1:n_g
                    val = du[alg_vars[gi]] + dt * d[stage] * dT[alg_vars[gi]]

                    # g_y * Σ_{j≤i} γ_{ij} * l_j
                    for fj in 1:n_f
                        gy_coupling_fj = C[stage, stage] * ks[stage][diff_vars[fj]]
                        for j in 1:(stage - 1)
                            gy_coupling_fj += C[stage, j] * ks[j][diff_vars[fj]]
                        end
                        val += g_y[gi, fj] * gy_coupling_fj
                    end

                    # g_z * Σ_{j<i} γ_{ij} * k_j
                    for gj in 1:n_g
                        gz_coupling_gj = zero(eltype(u))
                        for j in 1:(stage - 1)
                            gz_coupling_gj += C[stage, j] * ks[j][alg_vars[gj]]
                        end
                        val += g_z[gi, gj] * gz_coupling_gj
                    end

                    linsolve_tmp_z[gi] = val
                end

                # Solve W_z * k_alg = rhs  (W_z = -γ*g_z)
                sol_z = W_z \ linsolve_tmp_z
                for (gi, ai) in enumerate(alg_vars)
                    ks[stage][ai] = sol_z[gi]
                end
            end
        end
    end

    # Solution update: u = uprev + sum_i b[i] * ks[i]
    u .= uprev
    for i in 1:num_stages
        @.. u += b[i] * ks[i]
    end

    step_limiter!(u, integrator, p, t + dt)

    # Error estimation
    if integrator.opts.adaptive
        du .= 0
        for i in 1:num_stages
            @.. du += (b[i] - bhat[i]) * ks[i]
        end
        calculate_residuals!(
            atmp, du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    # Dense output
    if integrator.opts.calck
        for j in eachindex(integrator.k)
            integrator.k[j] .= 0
        end
        for i in eachindex(ks)
            for j in eachindex(integrator.k)
                @.. integrator.k[j] += H[j, i] * ks[i]
            end
        end
    end
    return nothing
end
