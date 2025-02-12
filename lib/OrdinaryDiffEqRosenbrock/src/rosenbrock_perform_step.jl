function initialize!(integrator, cache::Union{Rosenbrock23Cache,
        Rosenbrock32Cache})
    integrator.kshortsize = 2
    @unpack k₁, k₂, fsalfirst, fsallast = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = k₁
    integrator.k[2] = k₂
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function initialize!(integrator,
        cache::Union{Rosenbrock23ConstantCache,
            Rosenbrock32ConstantCache})
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = zero(integrator.fsalfirst)
    integrator.k[2] = zero(integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock23Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p, opts = integrator
    @unpack k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W, tmp, uf, tf, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter! = cache
    @unpack c₃₂, d = cache.tab

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

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtγ))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtγ))
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
            @.. broadcast=false linsolve_tmp=fsallast - c₃₂ * (k₂ - f₁) -
                                             2(k₁ - fsalfirst) + dt * dT
        else
            @.. broadcast=false du2=c₃₂ * k₂ + 2k₁
            mul!(_vec(du1), mass_matrix, _vec(du2))
            @.. broadcast=false linsolve_tmp=fsallast - du1 + c₃₂ * f₁ + 2fsalfirst +
                                             dt * dT
        end

        linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
        vecu = _vec(linres.u)
        veck3 = _vec(k₃)
        @.. veck3 = vecu * neginvdtγ

        integrator.stats.nsolve += 1

        if mass_matrix === I
            @.. broadcast=false tmp=dto6 * (k₁ - 2 * k₂ + k₃)
        else
            veck₁ = _vec(k₁)
            veck₂ = _vec(k₂)
            veck₃ = _vec(k₃)
            vectmp = _vec(tmp)
            @.. broadcast=false vectmp=ifelse(cache.algebraic_vars,
                false, dto6 * (veck₁ - 2 * veck₂ + veck₃))
        end
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
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
    @unpack t, dt, uprev, u, f, p, opts = integrator
    @unpack k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W, tmp, uf, tf, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter! = cache
    @unpack c₃₂, d = cache.tab

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

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtγ))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtγ))
    end

    vecu = _vec(linres.u)
    veck₁ = _vec(k₁)

    @.. veck₁ = vecu * neginvdtγ
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + dto2 * k₁
    stage_limiter!(u, integrator, p, t + dto2)
    f(f₁, u, p, t + dto2)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        tmp .= k₁
    else
        mul!(_vec(tmp), mass_matrix, _vec(k₁))
    end

    @.. broadcast=false linsolve_tmp=f₁ - tmp

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
        @.. broadcast=false linsolve_tmp=fsallast - c₃₂ * (k₂ - f₁) - 2(k₁ - fsalfirst) +
                                         dt * dT
    else
        @.. broadcast=false du2=c₃₂ * k₂ + 2k₁
        mul!(_vec(du1), mass_matrix, _vec(du2))
        @.. broadcast=false linsolve_tmp=fsallast - du1 + c₃₂ * f₁ + 2fsalfirst + dt * dT
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck3 = _vec(k₃)

    @.. veck3 = vecu * neginvdtγ
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + dto6 * (k₁ + 4k₂ + k₃)

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        @.. broadcast=false tmp=dto6 * (k₁ - 2 * k₂ + k₃)
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            invatol = inv(integrator.opts.abstol)
            @.. atmp = ifelse(cache.algebraic_vars, fsallast, false) * invatol
            integrator.EEst += integrator.opts.internalnorm(atmp, t)
        end
    end
    cache.linsolve = linres.cache
end

@muladd function perform_step!(integrator, cache::Rosenbrock23ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack c₃₂, d, tf, uf = cache

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
            linsolve_tmp = @.. (integrator.fsallast - c₃₂ * (k₂ - f₁) -
                                2 * (k₁ - integrator.fsalfirst) + dt * dT)
        else
            linsolve_tmp = mass_matrix * (@.. c₃₂ * k₂ + 2 * k₁)
            linsolve_tmp = @.. (integrator.fsallast - linsolve_tmp +
                                c₃₂ * f₁ + 2 * integrator.fsalfirst + dt * dT)
        end
        k₃ = _reshape(W \ _vec(linsolve_tmp), axes(uprev)) * neginvdtγ
        integrator.stats.nsolve += 1

        if u isa Number
            utilde = dto6 * f.mass_matrix[1, 1] * (k₁ - 2 * k₂ + k₃)
        else
            utilde = f.mass_matrix * (@.. dto6 * (k₁ - 2 * k₂ + k₃))
        end
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
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

@muladd function perform_step!(integrator, cache::Rosenbrock32ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack c₃₂, d, tf, uf = cache

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
        linsolve_tmp = @.. (integrator.fsallast - c₃₂ * (k₂ - f₁) -
                            2(k₁ - integrator.fsalfirst) + dt * dT)
    else
        linsolve_tmp = mass_matrix * (@.. c₃₂ * k₂ + 2 * k₁)
        linsolve_tmp = @.. (integrator.fsallast - linsolve_tmp +
                            c₃₂ * f₁ + 2 * integrator.fsalfirst + dt * dT)
    end
    k₃ = _reshape(W \ _vec(linsolve_tmp), axes(uprev)) * neginvdtγ
    integrator.stats.nsolve += 1
    u = @.. uprev + dto6 * (k₁ + 4k₂ + k₃)

    if integrator.opts.adaptive
        utilde = @.. dto6 * (k₁ - 2k₂ + k₃)
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            invatol = inv(integrator.opts.abstol)
            atmp = ifelse(integrator.differential_vars, false, integrator.fsallast) .*
                   invatol
            integrator.EEst += integrator.opts.internalnorm(atmp, t)
        end
    end

    integrator.k[1] = k₁
    integrator.k[2] = k₂
    integrator.u = u
    return nothing
end

function initialize!(integrator,
        cache::Union{Rosenbrock33ConstantCache,
            Rosenbrock34ConstantCache,
            Rosenbrock4ConstantCache})
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator,
        cache::Union{Rosenbrock33Cache,
            Rosenbrock34Cache,
            Rosenbrock4Cache})
    integrator.kshortsize = 2
    @unpack fsalfirst, fsallast = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = fsalfirst
    integrator.k[2] = fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::Rosenbrock33ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tf, uf = cache
    @unpack a21, a31, a32, C21, C31, C32, b1, b2, b3, btilde1, btilde2, btilde3, gamma, c2, c3, d1, d2, d3 = cache.tab

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
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return nothing
end

@muladd function perform_step!(integrator, cache::Rosenbrock33Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack du, du1, du2, fsalfirst, fsallast, k1, k2, k3, dT, J, W, uf, tf, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter! = cache
    @unpack a21, a31, a32, C21, C31, C32, b1, b2, b3, btilde1, btilde2, btilde3, gamma, c2, c3, d1, d2, d3 = cache.tab

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

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    end

    vecu = _vec(linres.u)
    veck1 = _vec(k1)

    @.. broadcast=false veck1=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a21 * k1
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(du, u, p, t + c2 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd2 * dT + dtC21 * k1
    else
        @.. broadcast=false du1=dtC21 * k1
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd2 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck2 = _vec(k2)

    @.. broadcast=false veck2=-vecu

    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a31 * k1 + a32 * k2
    stage_limiter!(u, integrator, p, t + c3 * dt)
    f(du, u, p, t + c3 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd3 * dT + dtC31 * k1 + dtC32 * k2
    else
        @.. broadcast=false du1=dtC31 * k1 + dtC32 * k2
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd3 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck3 = _vec(k3)

    @.. broadcast=false veck3=-vecu

    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + b1 * k1 + b2 * k2 + b3 * k3

    step_limiter!(u, integrator, p, t + dt)

    f(fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        @.. broadcast=false utilde=btilde1 * k1 + btilde2 * k2 + btilde3 * k3
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    cache.linsolve = linres.cache
end

################################################################################

@muladd function perform_step!(integrator, cache::Rosenbrock34ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tf, uf = cache
    @unpack a21, a31, a32, a41, a42, a43, C21, C31, C32, C41, C42, C43, b1, b2, b3, b4, btilde1, btilde2, btilde3, btilde4, gamma, c2, c3, d1, d2, d3, d4 = cache.tab

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
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return nothing
end

@muladd function perform_step!(integrator, cache::Rosenbrock34Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack du, du1, du2, fsalfirst, fsallast, k1, k2, k3, k4, dT, J, W, uf, tf, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter! = cache
    @unpack a21, a31, a32, a41, a42, a43, C21, C31, C32, C41, C42, C43, b1, b2, b3, b4, btilde1, btilde2, btilde3, btilde4, gamma, c2, c3, d1, d2, d3, d4 = cache.tab

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

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    end

    vecu = _vec(linres.u)
    veck1 = _vec(k1)

    @.. broadcast=false veck1=-vecu
    integrator.stats.nsolve += 1

    #=
    a21 == 0 and c2 == 0
    so du = integrator.fsalfirst!
    @.. broadcast=false u = uprev + a21*k1

    f(du, u, p, t+c2*dt)
    =#

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=fsalfirst + dtd2 * dT + dtC21 * k1
    else
        @.. broadcast=false du1=dtC21 * k1
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=fsalfirst + dtd2 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck2 = _vec(k2)
    @.. broadcast=false veck2=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a31 * k1 + a32 * k2
    stage_limiter!(u, integrator, p, t + c3 * dt)
    f(du, u, p, t + c3 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd3 * dT + dtC31 * k1 + dtC32 * k2
    else
        @.. broadcast=false du1=dtC31 * k1 + dtC32 * k2
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd3 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck3 = _vec(k3)
    @.. broadcast=false veck3=-vecu
    integrator.stats.nsolve += 1
    @.. broadcast=false u=uprev + a41 * k1 + a42 * k2 + a43 * k3
    stage_limiter!(u, integrator, p, t + dt)
    f(du, u, p, t + dt) #-- c4 = 1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd4 * dT + dtC41 * k1 + dtC42 * k2 +
                                         dtC43 * k3
    else
        @.. broadcast=false du1=dtC41 * k1 + dtC42 * k2 + dtC43 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd4 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck4 = _vec(k4)
    @.. broadcast=false veck4=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4

    step_limiter!(u, integrator, p, t + dt)

    f(fsallast, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        @.. broadcast=false utilde=btilde1 * k1 + btilde2 * k2 + btilde3 * k3 + btilde4 * k4
        calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end
    cache.linsolve = linres.cache
end

################################################################################

#### ROS2 type method

@ROS2(:init)
@ROS2(:performstep)

################################################################################

#### ROS23 type method

@ROS23(:init)
@ROS23(:performstep)

################################################################################

#### ROS34PW type method

@ROS34PW(:init)
@ROS34PW(:performstep)

################################################################################

#### ROS4 type method

@Rosenbrock4(:performstep)

################################################################################

#### Rodas3P type method

function initialize!(integrator, cache::Union{Rodas23WConstantCache, Rodas3PConstantCache})
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    # Avoid undefined entries if k is an array of arrays
    integrator.k[1] = zero(integrator.u)
    integrator.k[2] = zero(integrator.u)
    integrator.k[3] = zero(integrator.u)
end

@muladd function perform_step!(
        integrator, cache::Union{Rodas23WConstantCache, Rodas3PConstantCache},
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tf, uf = cache
    @unpack a21, a41, a42, a43, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, gamma, c2, c3, d1, d2, d3 = cache.tab

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
        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35, h2_21, h2_22, h2_23, h2_24, h2_25 = cache.tab
        integrator.k[1] = h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 + h25 * k5
        integrator.k[2] = h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 + h35 * k5
        integrator.k[3] = h2_21 * k1 + h2_22 * k2 + h2_23 * k3 + h2_24 * k4 + h2_25 * k5
        if integrator.opts.adaptive
            if isa(linsolve_tmp, AbstractFloat)
                u_int, u_diff = calculate_interpoldiff(
                    uprev, du, u, integrator.k[1], integrator.k[2], integrator.k[3])
            else
                u_int = linsolve_tmp
                u_diff = linsolve_tmp .+ 0
                calculate_interpoldiff!(u_int, u_diff, uprev, du, u, integrator.k[1],
                    integrator.k[2], integrator.k[3])
            end
            atmp = calculate_residuals(u_diff, uprev, u_int, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
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
        atmp = calculate_residuals(u - du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = max(EEst, integrator.opts.internalnorm(atmp, t))
    end

    integrator.u = u
    return nothing
end

function initialize!(integrator, cache::Union{Rodas23WCache, Rodas3PCache})
    integrator.kshortsize = 3
    @unpack dense1, dense2, dense3 = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = dense1
    integrator.k[2] = dense2
    integrator.k[3] = dense3
end

@muladd function perform_step!(
        integrator, cache::Union{Rodas23WCache, Rodas3PCache}, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack du, du1, du2, dT, J, W, uf, tf, k1, k2, k3, k4, k5, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter! = cache
    @unpack a21, a41, a42, a43, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, gamma, c2, c3, d1, d2, d3 = cache.tab

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

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = cache.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = cache.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    end

    @.. broadcast=false $(_vec(k1))=-linres.u

    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a21 * k1
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(du, u, p, t + c2 * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd2 * dT + dtC21 * k1
    else
        @.. broadcast=false du1=dtC21 * k1
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd2 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k2))=-linres.u
    integrator.stats.nsolve += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=cache.fsalfirst + dtd3 * dT +
                                         (dtC31 * k1 + dtC32 * k2)
    else
        @.. broadcast=false du1=dtC31 * k1 + dtC32 * k2
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=cache.fsalfirst + dtd3 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k3))=-linres.u
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a41 * k1 + a42 * k2 + a43 * k3
    stage_limiter!(u, integrator, p, t + c2 * dt)
    f(du, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du +
                                         (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    else
        @.. broadcast=false du1=dtC41 * k1 + dtC42 * k2 + dtC43 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k4))=-linres.u
    integrator.stats.nsolve += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du +
                                         (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    else
        @.. broadcast=false du1=dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k5))=-linres.u
    integrator.stats.nsolve += 1

    du = u + k4 #-- p=2 solution
    u .+= k5

    step_limiter!(u, integrator, p, t + dt)

    EEst = 0.0
    if integrator.opts.calck
        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35, h2_21, h2_22, h2_23, h2_24, h2_25 = cache.tab
        @.. broadcast=false integrator.k[1]=h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 +
                                            h25 * k5
        @.. broadcast=false integrator.k[2]=h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 +
                                            h35 * k5
        @.. broadcast=false integrator.k[3]=h2_21 * k1 + h2_22 * k2 + h2_23 * k3 +
                                            h2_24 * k4 + h2_25 * k5
        if integrator.opts.adaptive
            calculate_interpoldiff!(
                du1, du2, uprev, du, u, integrator.k[1], integrator.k[2], integrator.k[3])
            calculate_residuals!(atmp, du2, uprev, du1, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
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
        calculate_residuals!(atmp, u - du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
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
end

#### Rodas4 type method

function initialize!(integrator, cache::RosenbrockCombinedConstantCache)
    integrator.kshortsize = size(cache.tab.H, 1)
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    # Avoid undefined entries if k is an array of arrays
    for i in 1:(integrator.kshortsize)
        integrator.k[i] = zero(integrator.u)
    end
end

@muladd function perform_step!(
        integrator, cache::RosenbrockCombinedConstantCache, repeat_step = false)
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
    ks = ntuple(Returns(k1), 10)
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
    #@show ks
    u = u .+ ks[num_stages]

    if integrator.opts.adaptive
        atmp = calculate_residuals(ks[num_stages], uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
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
            k2 = 0.5 * (uprev + u +
                  0.5 * (integrator.k[1] + 0.5 * (integrator.k[2] + 0.5 * integrator.k[3])))
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
end

@muladd function perform_step!(integrator, cache::RosenbrockCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; du, du1, du2, dT, J, W, uf, tf, ks, linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter!) = cache
    (; A, C, gamma, c, d, H) = cache.tab

    # Assignments
    sizeu = size(u)
    uidx = eachindex(integrator.uprev)
    mass_matrix = integrator.f.mass_matrix

    # Precalculations
    dtC = C .* inv(dt)
    dtd = dt .* d
    dtgamma = dt * gamma

    f(cache.fsalfirst, uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    calc_rosenbrock_differentiation!(integrator, cache, dtd[1], dtgamma, repeat_step)

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = cache.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = cache.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = dtgamma))
    end

    @.. $(_vec(ks[1])) = -linres.u
    integrator.stats.nsolve += 1

    for stage in 2:length(ks)
        u .= uprev
        for i in 1:(stage - 1)
            @.. u += A[stage, i] * ks[i]
        end

        stage_limiter!(u, integrator, p, t + c[stage] * dt)
        f(du, u, p, t + c[stage] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        du1 .= 0
        if mass_matrix === I
            for i in 1:(stage - 1)
                @.. du1 += dtC[stage, i] * ks[i]
            end
        else
            for i in 1:(stage - 1)
                @.. du1 += dtC[stage, i] * ks[i]
            end
            mul!(_vec(du2), mass_matrix, _vec(du1))
            du1 .= du2
        end
        @.. linsolve_tmp = du + dtd[stage] * dT + du1

        linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
        @.. $(_vec(ks[stage])) = -linres.u
        integrator.stats.nsolve += 1
    end
    du .= ks[end]
    u .+= ks[end]

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        if (integrator.alg isa Rodas5Pe)
            @.. du = 0.2606326497975715 * ks[1] - 0.005158627295444251 * ks[2] +
                     1.3038988631109731 * ks[3] + 1.235000722062074 * ks[4] +
                     -0.7931985603795049 * ks[5] - 1.005448461135913 * ks[6] -
                     0.18044626132120234 * ks[7] + 0.17051519239113755 * ks[8]
        end
        calculate_residuals!(atmp, ks[end], uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
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
            ks[2] = 0.5 * (uprev + u +
                     0.5 *
                     (integrator.k[1] + 0.5 * (integrator.k[2] + 0.5 * integrator.k[3])))
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

@RosenbrockW6S4OS(:init)
@RosenbrockW6S4OS(:performstep)
