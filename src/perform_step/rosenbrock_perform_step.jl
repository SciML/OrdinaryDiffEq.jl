function initialize!(integrator, cache::Union{Rosenbrock23Cache,
        Rosenbrock32Cache})
    integrator.kshortsize = 2
    @unpack k₁, k₂, fsalfirst, fsallast = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = fsallast
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = k₁
    integrator.k[2] = k₂
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
end

function initialize!(integrator,
        cache::Union{Rosenbrock23ConstantCache,
            Rosenbrock32ConstantCache})
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = zero(integrator.fsalfirst)
    integrator.k[2] = zero(integrator.fsalfirst)
end

@muladd function perform_step!(integrator, cache::Rosenbrock23Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p, opts = integrator
    @unpack k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W, tmp, uf, tf, linsolve_tmp, jac_config, atmp, weight = cache
    @unpack c₃₂, d = cache.tab

    # Assignments
    sizeu = size(u)
    mass_matrix = integrator.f.mass_matrix

    # Precalculations
    γ = dt * d
    dto2 = dt / 2
    dto6 = dt / 6

    if repeat_step
        f(integrator.fsalfirst, uprev, p, t)
        integrator.stats.nf += 1
    end

    calc_rosenbrock_differentiation!(integrator, cache, γ, γ, repeat_step, false)

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = γ))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = γ))
    end

    vecu = _vec(linres.u)
    veck₁ = _vec(k₁)

    @.. broadcast=false veck₁=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + dto2 * k₁
    f(f₁, u, p, t + dto2)
    integrator.stats.nf += 1

    if mass_matrix === I
        copyto!(tmp, k₁)
    else
        mul!(_vec(tmp), mass_matrix, _vec(k₁))
    end

    @.. broadcast=false linsolve_tmp=f₁ - tmp

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck2 = _vec(k₂)

    @.. broadcast=false veck2=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false k₂+=k₁
    @.. broadcast=false u=uprev + dt * k₂

    if integrator.opts.adaptive
        f(fsallast, u, p, t + dt)
        integrator.stats.nf += 1

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
        @.. broadcast=false veck3=-vecu

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
            @.. broadcast=false atmp=ifelse(cache.algebraic_vars, fsallast, false) /
                                     integrator.opts.abstol
            integrator.EEst += integrator.opts.internalnorm(atmp, t)
        end
    end
    cache.linsolve = linres.cache
end

@muladd function perform_step!(integrator, cache::Rosenbrock23Cache{<:Array},
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p, opts = integrator
    @unpack k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W, tmp, uf, tf, linsolve_tmp, jac_config, atmp, weight = cache
    @unpack c₃₂, d = cache.tab

    # Assignments
    sizeu = size(u)
    mass_matrix = integrator.f.mass_matrix

    # Precalculations
    γ = dt * d
    dto2 = dt / 2
    dto6 = dt / 6

    if repeat_step
        f(integrator.fsalfirst, uprev, p, t)
        integrator.stats.nf += 1
    end

    calc_rosenbrock_differentiation!(integrator, cache, γ, γ, repeat_step, false)

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)
    linsolve = cache.linsolve

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = γ))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = γ))
    end

    @inbounds @simd ivdep for i in eachindex(u)
        k₁[i] = -linres.u[i]
    end

    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + dto2 * k₁[i]
    end
    f(f₁, u, p, t + dto2)
    integrator.stats.nf += 1

    if mass_matrix === I
        copyto!(tmp, k₁)
    else
        mul!(_vec(tmp), mass_matrix, _vec(k₁))
    end

    @inbounds @simd ivdep for i in eachindex(u)
        linsolve_tmp[i] = f₁[i] - tmp[i]
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))

    @inbounds @simd ivdep for i in eachindex(u)
        k₂[i] = -linres.u[i]
    end

    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        k₂[i] += k₁[i]
    end
    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + dt * k₂[i]
    end

    if integrator.opts.adaptive
        f(fsallast, u, p, t + dt)
        integrator.stats.nf += 1

        if mass_matrix === I
            @inbounds @simd ivdep for i in eachindex(u)
                linsolve_tmp[i] = fsallast[i] - c₃₂ * (k₂[i] - f₁[i]) -
                                  2(k₁[i] - fsalfirst[i]) + dt * dT[i]
            end
        else
            @inbounds @simd ivdep for i in eachindex(u)
                du2[i] = c₃₂ * k₂[i] + 2k₁[i]
            end
            mul!(_vec(du1), mass_matrix, _vec(du2))

            @inbounds @simd ivdep for i in eachindex(u)
                linsolve_tmp[i] = fsallast[i] - du1[i] + c₃₂ * f₁[i] + 2fsalfirst[i] +
                                  dt * dT[i]
            end
        end

        linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))

        @inbounds @simd ivdep for i in eachindex(u)
            k₃[i] = -linres.u[i]
        end
        integrator.stats.nsolve += 1

        if mass_matrix === I
            @inbounds @simd ivdep for i in eachindex(u)
                tmp[i] = dto6 * (k₁[i] - 2 * k₂[i] + k₃[i])
            end
        else
            veck₁ = _vec(k₁)
            veck₂ = _vec(k₂)
            veck₃ = _vec(k₃)
            vectmp = _vec(tmp)
            @inbounds @simd ivdep for i in 1:length(vectmp)
                vectmp[i] = ifelse(cache.algebraic_vars[i], false,
                    dto6 * (veck₁[i] - 2 * veck₂[i] + veck₃[i]))
            end
        end
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            if integrator.opts.abstol isa AbstractVector
                @inbounds @simd ivdep for i in 1:length(vectmp)
                    vectmp[i] = ifelse(cache.algebraic_vars[i], fsallast[i], false) /
                                integrator.opts.abstol[i]
                end
            else
                @inbounds @simd ivdep for i in 1:length(vectmp)
                    vectmp[i] = ifelse(cache.algebraic_vars[i], fsallast[i], false) /
                                integrator.opts.abstol
                end
            end
            integrator.EEst += integrator.opts.internalnorm(vectmp, t)
        end
    end
    cache.linsolve = linres.cache
end

@muladd function perform_step!(integrator, cache::Rosenbrock32Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p, opts = integrator
    @unpack k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W, tmp, uf, tf, linsolve_tmp, jac_config, atmp, weight = cache
    @unpack c₃₂, d = cache.tab

    # Assignments
    sizeu = size(u)
    mass_matrix = integrator.f.mass_matrix

    # Precalculations
    γ = dt * d
    dto2 = dt / 2
    dto6 = dt / 6

    if repeat_step
        f(integrator.fsalfirst, uprev, p, t)
        integrator.stats.nf += 1
    end

    calc_rosenbrock_differentiation!(integrator, cache, γ, γ, repeat_step, false)

    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t)

    if repeat_step
        linres = dolinsolve(
            integrator, cache.linsolve; A = nothing, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = γ))
    else
        linres = dolinsolve(integrator, cache.linsolve; A = W, b = _vec(linsolve_tmp),
            du = integrator.fsalfirst, u = u, p = p, t = t, weight = weight,
            solverdata = (; gamma = γ))
    end

    vecu = _vec(linres.u)
    veck₁ = _vec(k₁)

    @.. broadcast=false veck₁=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + dto2 * k₁
    f(f₁, u, p, t + dto2)
    integrator.stats.nf += 1

    if mass_matrix === I
        tmp .= k₁
    else
        mul!(_vec(tmp), mass_matrix, _vec(k₁))
    end

    @.. broadcast=false linsolve_tmp=f₁ - tmp

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    vecu = _vec(linres.u)
    veck2 = _vec(k₂)

    @.. broadcast=false veck2=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false k₂+=k₁
    @.. broadcast=false tmp=uprev + dt * k₂
    f(fsallast, tmp, p, t + dt)
    integrator.stats.nf += 1

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

    @.. broadcast=false veck3=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + dto6 * (k₁ + 4k₂ + k₃)

    if integrator.opts.adaptive
        @.. broadcast=false tmp=dto6 * (k₁ - 2 * k₂ + k₃)
        calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            @.. broadcast=false atmp=ifelse(cache.algebraic_vars, fsallast, false) /
                                     integrator.opts.abstol
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
    γ = dt * d
    dto2 = dt / 2
    dto6 = dt / 6

    if repeat_step
        integrator.fsalfirst = f(uprev, p, t)
        integrator.stats.nf += 1
    end

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, γ, repeat_step)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    k₁ = _reshape(W \ -_vec((integrator.fsalfirst + γ * dT)), axes(uprev))
    integrator.stats.nsolve += 1
    f₁ = f(uprev + dto2 * k₁, p, t + dto2)
    integrator.stats.nf += 1

    if mass_matrix === I
        k₂ = _reshape(W \ -_vec(f₁ - k₁), axes(uprev)) + k₁
    else
        k₂ = _reshape(W \ -_vec(f₁ - mass_matrix * k₁), axes(uprev)) + k₁
    end
    integrator.stats.nsolve += 1
    u = uprev + dt * k₂

    if integrator.opts.adaptive
        integrator.fsallast = f(u, p, t + dt)
        integrator.stats.nf += 1

        if mass_matrix === I
            k₃ = _reshape(
                W \
                -_vec((integrator.fsallast - c₃₂ * (k₂ - f₁) -
                       2 * (k₁ - integrator.fsalfirst) + dt * dT)),
                axes(uprev))
        else
            linsolve_tmp = integrator.fsallast - mass_matrix * (c₃₂ * k₂ + 2 * k₁) +
                           c₃₂ * f₁ + 2 * integrator.fsalfirst + dt * dT
            k₃ = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
        end
        integrator.stats.nsolve += 1

        if u isa Number
            utilde = dto6 * f.mass_matrix[1, 1] * (k₁ - 2 * k₂ + k₃)
        else
            utilde = dto6 * f.mass_matrix * (k₁ - 2 * k₂ + k₃)
        end
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            atmp = @. ifelse(!integrator.differential_vars, integrator.fsallast, false) ./
                      integrator.opts.abstol
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
    γ = dt * d
    dto2 = dt / 2
    dto6 = dt / 6

    mass_matrix = integrator.f.mass_matrix

    if repeat_step
        integrator.fsalfirst = f(uprev, p, t)
        integrator.stats.nf += 1
    end

    # Time derivative
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, γ, repeat_step)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    k₁ = _reshape(W \ -_vec((integrator.fsalfirst + γ * dT)), axes(uprev))
    integrator.stats.nsolve += 1
    f₁ = f(uprev + dto2 * k₁, p, t + dto2)
    integrator.stats.nf += 1

    if mass_matrix === I
        k₂ = _reshape(W \ -_vec(f₁ - k₁), axes(uprev)) + k₁
    else
        linsolve_tmp = f₁ - mass_matrix * k₁
        k₂ = _reshape(W \ -_vec(linsolve_tmp), axes(uprev)) + k₁
    end

    integrator.stats.nsolve += 1
    tmp = uprev + dt * k₂
    integrator.fsallast = f(tmp, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        k₃ = _reshape(
            W \
            -_vec((integrator.fsallast - c₃₂ * (k₂ - f₁) -
                   2(k₁ - integrator.fsalfirst) + dt * dT)),
            axes(uprev))
    else
        linsolve_tmp = integrator.fsallast - mass_matrix * (c₃₂ * k₂ + 2k₁) + c₃₂ * f₁ +
                       2 * integrator.fsalfirst + dt * dT
        k₃ = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    end
    integrator.stats.nsolve += 1
    u = uprev + dto6 * (k₁ + 4k₂ + k₃)

    if integrator.opts.adaptive
        utilde = dto6 * (k₁ - 2k₂ + k₃)
        atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)

        if mass_matrix !== I
            atmp = @. ifelse(!integrator.differential_vars, integrator.fsallast, false) ./
                      integrator.opts.abstol
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
    integrator.stats.nf += 1

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
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = fsallast
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = fsalfirst
    integrator.k[2] = fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
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

    W = calc_W(integrator, cache, dtgamma, repeat_step, true)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    linsolve_tmp = integrator.fsalfirst + dtd1 * dT

    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a21 * k1
    du = f(u, p, t + c2 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd2 * dT + dtC21 * k1
    else
        linsolve_tmp = du + dtd2 * dT + mass_matrix * (dtC21 * k1)
    end

    k2 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a31 * k1 + a32 * k2
    du = f(u, p, t + c3 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd3 * dT + dtC31 * k1 + dtC32 * k2
    else
        linsolve_tmp = du + dtd3 * dT + mass_matrix * (dtC31 * k1 + dtC32 * k2)
    end

    k3 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + b1 * k1 + b2 * k2 + b3 * k3
    integrator.fsallast = f(u, p, t + dt)
    integrator.stats.nf += 1

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
    @unpack du, du1, du2, fsalfirst, fsallast, k1, k2, k3, dT, J, W, uf, tf, linsolve_tmp, jac_config, atmp, weight = cache
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

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

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
    f(du, u, p, t + c2 * dt)
    integrator.stats.nf += 1

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
    f(du, u, p, t + c3 * dt)
    integrator.stats.nf += 1

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
    f(fsallast, u, p, t + dt)
    integrator.stats.nf += 1

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

    W = calc_W(integrator, cache, dtgamma, repeat_step, true)
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
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd3 * dT + dtC31 * k1 + dtC32 * k2
    else
        linsolve_tmp = du + dtd3 * dT + mass_matrix * (dtC31 * k1 + dtC32 * k2)
    end

    k3 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a41 * k1 + a42 * k2 + a43 * k3
    du = f(u, p, t + dt) #-- c4 = 1
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd4 * dT + dtC41 * k1 + dtC42 * k2 + dtC43 * k3
    else
        linsolve_tmp = du + dtd4 * dT + mass_matrix * (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    end

    k4 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4
    integrator.fsallast = f(u, p, t + dt)
    integrator.stats.nf += 1

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
    @unpack du, du1, du2, fsalfirst, fsallast, k1, k2, k3, k4, dT, J, W, uf, tf, linsolve_tmp, jac_config, atmp, weight = cache
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

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

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
    f(du, u, p, t + c3 * dt)
    integrator.stats.nf += 1

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
    f(du, u, p, t + dt) #-- c4 = 1
    integrator.stats.nf += 1

    @.. broadcast=false u=uprev + a41 * k1 + a42 * k2 + a43 * k3
    f(du, u, p, t + dt) #-- c4 = 1
    integrator.stats.nf += 1

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
    f(fsallast, u, p, t + dt)
    integrator.stats.nf += 1

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

    W = calc_W(integrator, cache, dtgamma, repeat_step, true)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    du = f(uprev, p, t)
    integrator.stats.nf += 1
    k3 = copy(du)  #-- save for stage 3

    linsolve_tmp = du + dtd1 * dT

    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a21 * k1
    du = f(u, p, t + c2 * dt)
    integrator.stats.nf += 1

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
    integrator.stats.nf += 1

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
            if  isa(linsolve_tmp,AbstractFloat)
                u_int, u_diff = calculate_interpoldiff(uprev, du, u, integrator.k[1], integrator.k[2], integrator.k[3])
            else
                u_int = linsolve_tmp 
                u_diff = linsolve_tmp .+ 0
                calculate_interpoldiff!(u_int, u_diff, uprev, du, u, integrator.k[1], integrator.k[2], integrator.k[3])
            end
            atmp = calculate_residuals(u_diff, uprev, u_int, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            EEst = max(EEst,integrator.opts.internalnorm(atmp, t))  #-- role of t unclear
        end
    end

    if (integrator.alg isa Rodas23W)
        k1 = u .+ 0
        u = du .+ 0
        du = k1 .+ 0
        if integrator.opts.calck
            integrator.k[1] = integrator.k[3] .+ 0
            integrator.k[2] = 0*integrator.k[2]
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
    @unpack du, du1, du2, dT, J, W, uf, tf, k1, k2, k3, k4, k5, linsolve_tmp, jac_config, atmp, weight = cache
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
    integrator.stats.nf += 1

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

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
    f(du, u, p, t + c2 * dt)
    integrator.stats.nf += 1

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
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

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

@muladd function perform_step!(
        integrator, cache::Union{Rodas23WCache{<:Array}, Rodas3PCache{<:Array}},
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack du, du1, du2, dT, J, W, uf, tf, k1, k2, k3, k4, k5, linsolve_tmp, jac_config, atmp, weight = cache
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
    integrator.stats.nf += 1

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

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

    @inbounds @simd ivdep for i in eachindex(u)
        k1[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + a21 * k1[i]
    end
    f(du, u, p, t + c2 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd2 * dT[i] + dtC21 * k1[i]
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC21 * k1[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd2 * dT[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k2[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = cache.fsalfirst[i] + dtd3 * dT[i] +
                              (dtC31 * k1[i] + dtC32 * k2[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC31 * k1[i] + dtC32 * k2[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = cache.fsalfirst[i] + dtd3 * dT[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k3[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + a41 * k1[i] + a42 * k2[i] + a43 * k3[i]
    end
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] +
                              (dtC41 * k1[i] + dtC42 * k2[i] + dtC43 * k3[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC41 * k1[i] + dtC42 * k2[i] + dtC43 * k3[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k4[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + (dtC52 * k2[i] + dtC54 * k4[i] + dtC51 * k1[i] +
                               dtC53 * k3[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC52 * k2[i] + dtC54 * k4[i] + dtC51 * k1[i] + dtC53 * k3[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k5[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)  #-- sol p=2
        du[i] = u[i] + k4[i]
    end
    @inbounds @simd ivdep for i in eachindex(u)
        u[i] += k5[i]
    end

    EEst = 0.0
    if integrator.opts.calck
        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35, h2_21, h2_22, h2_23, h2_24, h2_25 = cache.tab
        @inbounds @simd ivdep for i in eachindex(u)
            integrator.k[1][i] = h21 * k1[i] + h22 * k2[i] + h23 * k3[i] + h24 * k4[i] +
                                 h25 * k5[i]
            integrator.k[2][i] = h31 * k1[i] + h32 * k2[i] + h33 * k3[i] + h34 * k4[i] +
                                 h35 * k5[i]
            integrator.k[3][i] = h2_21 * k1[i] + h2_22 * k2[i] + h2_23 * k3[i] +
                                 h2_24 * k4[i] + h2_25 * k5[i]
        end
        if integrator.opts.adaptive
            calculate_interpoldiff!(
                du1, du2, uprev, du, u, integrator.k[1], integrator.k[2], integrator.k[3])
            calculate_residuals!(atmp, du2, uprev, du1, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            EEst = max(EEst, integrator.opts.internalnorm(atmp, t))  #-- role of t unclear
        end
    end

    if (integrator.alg isa Rodas23W)
        @inbounds @simd ivdep for i in eachindex(u)
            tt = u[i]
            u[i] = du[i]
            du[i] = tt
            if integrator.opts.calck
                integrator.k[1][i] = integrator.k[3][i]
                integrator.k[2][i] = 0.0
            end
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
    dis = a2^2 - 3*a1*a3
    u_int = up3
    u_diff = 0.0
    if dis > 0.0 #-- Min/Max occurs
        tau1 = (-a2 - sqrt(dis))/(3*a3)
        tau2 = (-a2 + sqrt(dis))/(3*a3)
        if tau1 > tau2 
            tau1,tau2 = tau2,tau1 
        end
        for tau in (tau1,tau2)
            if (tau > 0.0) && (tau < 1.0)
                y_tau = (1 - tau)*uprev + 
                        tau*(up3 + (1 - tau)*(c_koeff + tau*d_koeff))
                dy_tau = ((a3*tau + a2)*tau + a1)*tau
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

function initialize!(integrator, cache::Rodas4ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    # Avoid undefined entries if k is an array of arrays
    integrator.k[1] = zero(integrator.u)
    integrator.k[2] = zero(integrator.u)
end

@muladd function perform_step!(integrator, cache::Rodas4ConstantCache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tf, uf = cache
    @unpack a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, C61, C62, C63, C64, C65, gamma, c2, c3, c4, d1, d2, d3, d4 = cache.tab

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
    dtC61 = C61 / dt
    dtC62 = C62 / dt
    dtC63 = C63 / dt
    dtC64 = C64 / dt
    dtC65 = C65 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtd4 = dt * d4
    dtgamma = dt * gamma

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    tf.u = uprev
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step, true)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    du = f(uprev, p, t)
    integrator.stats.nf += 1

    linsolve_tmp = du + dtd1 * dT

    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a21 * k1
    du = f(u, p, t + c2 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd2 * dT + dtC21 * k1
    else
        linsolve_tmp = du + dtd2 * dT + mass_matrix * (dtC21 * k1)
    end

    k2 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a31 * k1 + a32 * k2
    du = f(u, p, t + c3 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd3 * dT + (dtC31 * k1 + dtC32 * k2)
    else
        linsolve_tmp = du + dtd3 * dT + mass_matrix * (dtC31 * k1 + dtC32 * k2)
    end

    k3 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a41 * k1 + a42 * k2 + a43 * k3
    du = f(u, p, t + c4 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd4 * dT + (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    else
        linsolve_tmp = du + dtd4 * dT + mass_matrix * (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    end

    k4 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4
    du = f(u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    else
        linsolve_tmp = du +
                       mass_matrix * (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    end

    k5 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = u + k5
    du = f(u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + (dtC61 * k1 + dtC62 * k2 + dtC65 * k5 + dtC64 * k4 + dtC63 * k3)
    else
        linsolve_tmp = du +
                       mass_matrix *
                       (dtC61 * k1 + dtC62 * k2 + dtC65 * k5 + dtC64 * k4 + dtC63 * k3)
    end

    k6 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = u + k6

    if integrator.opts.adaptive
        atmp = calculate_residuals(k6, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35 = cache.tab
        integrator.k[1] = h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 + h25 * k5
        integrator.k[2] = h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 + h35 * k5
    end
    integrator.u = u
    return nothing
end

function initialize!(integrator, cache::Rodas4Cache)
    integrator.kshortsize = 2
    @unpack dense1, dense2 = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = dense1
    integrator.k[2] = dense2
end

@muladd function perform_step!(integrator, cache::Rodas4Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack du, du1, du2, dT, J, W, uf, tf, k1, k2, k3, k4, k5, k6, linsolve_tmp, jac_config, atmp, weight = cache
    @unpack a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, C61, C62, C63, C64, C65, gamma, c2, c3, c4, d1, d2, d3, d4 = cache.tab

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
    dtC61 = C61 / dt
    dtC62 = C62 / dt
    dtC63 = C63 / dt
    dtC64 = C64 / dt
    dtC65 = C65 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtd4 = dt * d4
    dtgamma = dt * gamma

    f(cache.fsalfirst, uprev, p, t) # used in calc_rosenbrock_differentiation!
    integrator.stats.nf += 1

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

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
    f(du, u, p, t + c2 * dt)
    integrator.stats.nf += 1

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

    @.. broadcast=false u=uprev + a31 * k1 + a32 * k2
    f(du, u, p, t + c3 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd3 * dT + (dtC31 * k1 + dtC32 * k2)
    else
        @.. broadcast=false du1=dtC31 * k1 + dtC32 * k2
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd3 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k3))=-linres.u
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a41 * k1 + a42 * k2 + a43 * k3
    f(du, u, p, t + c4 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd4 * dT +
                                         (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    else
        @.. broadcast=false du1=dtC41 * k1 + dtC42 * k2 + dtC43 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd4 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k4))=-linres.u
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

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

    u .+= k5
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + (dtC61 * k1 + dtC62 * k2 + dtC65 * k5 +
                                          dtC64 * k4 + dtC63 * k3)
    else
        @.. broadcast=false du1=dtC61 * k1 + dtC62 * k2 + dtC65 * k5 + dtC64 * k4 +
                                dtC63 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. broadcast=false $(_vec(k6))=-linres.u
    integrator.stats.nsolve += 1

    u .+= k6

    if integrator.opts.adaptive
        calculate_residuals!(atmp, k6, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35 = cache.tab
        @.. broadcast=false integrator.k[1]=h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 +
                                            h25 * k5
        @.. broadcast=false integrator.k[2]=h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 +
                                            h35 * k5
    end
    cache.linsolve = linres.cache
end

@muladd function perform_step!(integrator, cache::Rodas4Cache{<:Array}, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack du, du1, du2, dT, J, W, uf, tf, k1, k2, k3, k4, k5, k6, linsolve_tmp, jac_config, atmp, weight = cache
    @unpack a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, C61, C62, C63, C64, C65, gamma, c2, c3, c4, d1, d2, d3, d4 = cache.tab

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
    dtC61 = C61 / dt
    dtC62 = C62 / dt
    dtC63 = C63 / dt
    dtC64 = C64 / dt
    dtC65 = C65 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtd4 = dt * d4
    dtgamma = dt * gamma

    f(cache.fsalfirst, uprev, p, t) # used in calc_rosenbrock_differentiation!
    integrator.stats.nf += 1

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

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

    @inbounds @simd ivdep for i in eachindex(u)
        k1[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + a21 * k1[i]
    end
    f(du, u, p, t + c2 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd2 * dT[i] + dtC21 * k1[i]
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC21 * k1[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd2 * dT[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k2[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + a31 * k1[i] + a32 * k2[i]
    end
    f(du, u, p, t + c3 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd3 * dT[i] + (dtC31 * k1[i] + dtC32 * k2[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC31 * k1[i] + dtC32 * k2[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd3 * dT[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k3[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + a41 * k1[i] + a42 * k2[i] + a43 * k3[i]
    end
    f(du, u, p, t + c4 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd4 * dT[i] +
                              (dtC41 * k1[i] + dtC42 * k2[i] + dtC43 * k3[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC41 * k1[i] + dtC42 * k2[i] + dtC43 * k3[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd4 * dT[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k4[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]
    end
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + (dtC52 * k2[i] + dtC54 * k4[i] + dtC51 * k1[i] +
                               dtC53 * k3[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC52 * k2[i] + dtC54 * k4[i] + dtC51 * k1[i] + dtC53 * k3[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k5[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] += k5[i]
    end
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + (dtC61 * k1[i] + dtC62 * k2[i] + dtC65 * k5[i] +
                               dtC64 * k4[i] + dtC63 * k3[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC61 * k1[i] + dtC62 * k2[i] + dtC65 * k5[i] + dtC64 * k4[i] +
                     dtC63 * k3[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k6[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] += k6[i]
    end

    if integrator.opts.adaptive
        calculate_residuals!(atmp, k6, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        @unpack h21, h22, h23, h24, h25, h31, h32, h33, h34, h35 = cache.tab
        @inbounds @simd ivdep for i in eachindex(u)
            integrator.k[1][i] = h21 * k1[i] + h22 * k2[i] + h23 * k3[i] + h24 * k4[i] +
                                 h25 * k5[i]
            integrator.k[2][i] = h31 * k1[i] + h32 * k2[i] + h33 * k3[i] + h34 * k4[i] +
                                 h35 * k5[i]
        end
    end
    cache.linsolve = linres.cache
end

###############################################################################

### Rodas5 Method

function initialize!(integrator, cache::Rosenbrock5ConstantCache)
    integrator.kshortsize = 3
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    # Avoid undefined entries if k is an array of arrays
    integrator.k[1] = zero(integrator.u)
    integrator.k[2] = zero(integrator.u)
    integrator.k[3] = zero(integrator.u)
end

@muladd function perform_step!(integrator, cache::Rosenbrock5ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tf, uf = cache
    @unpack a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, C61, C62, C63, C64, C65, C71, C72, C73, C74, C75, C76, C81, C82, C83, C84, C85, C86, C87, gamma, d1, d2, d3, d4, d5, c2, c3, c4, c5 = cache.tab

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
    dtC61 = C61 / dt
    dtC62 = C62 / dt
    dtC63 = C63 / dt
    dtC64 = C64 / dt
    dtC65 = C65 / dt
    dtC71 = C71 / dt
    dtC72 = C72 / dt
    dtC73 = C73 / dt
    dtC74 = C74 / dt
    dtC75 = C75 / dt
    dtC76 = C76 / dt
    dtC81 = C81 / dt
    dtC82 = C82 / dt
    dtC83 = C83 / dt
    dtC84 = C84 / dt
    dtC85 = C85 / dt
    dtC86 = C86 / dt
    dtC87 = C87 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtd4 = dt * d4
    dtd5 = dt * d5
    dtgamma = dt * gamma

    mass_matrix = integrator.f.mass_matrix

    # Time derivative
    dT = calc_tderivative(integrator, cache)

    W = calc_W(integrator, cache, dtgamma, repeat_step, true)
    if !issuccess_W(W)
        integrator.EEst = 2
        return nothing
    end

    du1 = f(uprev, p, t)
    integrator.stats.nf += 1

    linsolve_tmp = du1 + dtd1 * dT

    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a21 * k1
    du = f(u, p, t + c2 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd2 * dT + dtC21 * k1
    else
        linsolve_tmp = du + dtd2 * dT + mass_matrix * (dtC21 * k1)
    end

    k2 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a31 * k1 + a32 * k2
    du = f(u, p, t + c3 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd3 * dT + (dtC31 * k1 + dtC32 * k2)
    else
        linsolve_tmp = du + dtd3 * dT + mass_matrix * (dtC31 * k1 + dtC32 * k2)
    end

    k3 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a41 * k1 + a42 * k2 + a43 * k3
    du = f(u, p, t + c4 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd4 * dT + (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    else
        linsolve_tmp = du + dtd4 * dT + mass_matrix * (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    end

    k4 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4
    du = f(u, p, t + c5 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + dtd5 * dT + (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    else
        linsolve_tmp = du + dtd5 * dT +
                       mass_matrix * (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    end

    k5 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = uprev + a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5
    du = f(u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du + (dtC61 * k1 + dtC62 * k2 + dtC63 * k3 + dtC64 * k4 + dtC65 * k5)
    else
        linsolve_tmp = du +
                       mass_matrix *
                       (dtC61 * k1 + dtC62 * k2 + dtC63 * k3 + dtC64 * k4 + dtC65 * k5)
    end

    k6 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = u + k6
    du = f(u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du +
                       (dtC71 * k1 + dtC72 * k2 + dtC73 * k3 + dtC74 * k4 + dtC75 * k5 +
                        dtC76 * k6)
    else
        linsolve_tmp = du +
                       mass_matrix *
                       (dtC71 * k1 + dtC72 * k2 + dtC73 * k3 + dtC74 * k4 + dtC75 * k5 +
                        dtC76 * k6)
    end

    k7 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = u + k7
    du = f(u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        linsolve_tmp = du +
                       (dtC81 * k1 + dtC82 * k2 + dtC83 * k3 + dtC84 * k4 + dtC85 * k5 +
                        dtC86 * k6 + dtC87 * k7)
    else
        linsolve_tmp = du +
                       mass_matrix *
                       (dtC81 * k1 + dtC82 * k2 + dtC83 * k3 + dtC84 * k4 + dtC85 * k5 +
                        dtC86 * k6 + dtC87 * k7)
    end

    k8 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    integrator.stats.nsolve += 1
    u = u + k8
    linsolve_tmp = k8

    if integrator.opts.adaptive
	    if (integrator.alg isa Rodas5Pe)
            linsolve_tmp = 0.2606326497975715*k1 - 0.005158627295444251*k2 + 1.3038988631109731*k3 + 1.235000722062074*k4 +
               - 0.7931985603795049*k5 - 1.005448461135913*k6 - 0.18044626132120234*k7 + 0.17051519239113755*k8
	    end
        atmp = calculate_residuals(linsolve_tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        @unpack h21, h22, h23, h24, h25, h26, h27, h28, h31, h32, h33, h34, h35, h36, h37, h38, h41, h42, h43, h44, h45, h46, h47, h48 = cache.tab
        integrator.k[1] = h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 + h25 * k5 + h26 * k6 +
                          h27 * k7 + h28 * k8
        integrator.k[2] = h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 + h35 * k5 + h36 * k6 +
                          h37 * k7 + h38 * k8
        integrator.k[3] = h41 * k1 + h42 * k2 + h43 * k3 + h44 * k4 + h45 * k5 + h46 * k6 +
                          h47 * k7 + h48 * k8
        if (integrator.alg isa Rodas5Pr) && integrator.opts.adaptive && (integrator.EEst < 1.0)
            k2 = 0.5*(uprev + u + 0.5 * (integrator.k[1] + 0.5 * (integrator.k[2] + 0.5 * integrator.k[3])))
            du1 = ( 0.25*(integrator.k[2] + integrator.k[3]) - uprev + u) / dt
            du = f(k2, p, t + dt/2)
            integrator.stats.nf += 1
            if mass_matrix === I
                du2 = du1 - du
            else
                du2 = mass_matrix*du1 - du
            end
	        EEst = norm(du2) / (integrator.opts.abstol + integrator.opts.reltol*norm(k2)) 
            integrator.EEst = max(EEst,integrator.EEst)
        end
    end

    integrator.u = u
    return nothing
end

function initialize!(integrator, cache::Rosenbrock5Cache)
    integrator.kshortsize = 3
    @unpack dense1, dense2, dense3 = cache
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = dense1
    integrator.k[2] = dense2
    integrator.k[3] = dense3
end

@muladd function perform_step!(integrator, cache::Rosenbrock5Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack du, du1, du2, k1, k2, k3, k4, k5, k6, k7, k8, dT, J, W, uf, tf, linsolve_tmp, jac_config, atmp, weight = cache
    @unpack a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, C61, C62, C63, C64, C65, C71, C72, C73, C74, C75, C76, C81, C82, C83, C84, C85, C86, C87, gamma, d1, d2, d3, d4, d5, c2, c3, c4, c5 = cache.tab

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
    dtC61 = C61 / dt
    dtC62 = C62 / dt
    dtC63 = C63 / dt
    dtC64 = C64 / dt
    dtC65 = C65 / dt
    dtC71 = C71 / dt
    dtC72 = C72 / dt
    dtC73 = C73 / dt
    dtC74 = C74 / dt
    dtC75 = C75 / dt
    dtC76 = C76 / dt
    dtC81 = C81 / dt
    dtC82 = C82 / dt
    dtC83 = C83 / dt
    dtC84 = C84 / dt
    dtC85 = C85 / dt
    dtC86 = C86 / dt
    dtC87 = C87 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtd4 = dt * d4
    dtd5 = dt * d5
    dtgamma = dt * gamma

    f(cache.fsalfirst, uprev, p, t) # used in calc_rosenbrock_differentiation!
    integrator.stats.nf += 1

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

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

    vecu = _vec(linres.u)
    veck1 = _vec(k1)

    @.. broadcast=false veck1=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a21 * k1
    f(du, u, p, t + c2 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd2 * dT + dtC21 * k1
    else
        @.. broadcast=false du1=dtC21 * k1
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd2 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck2 = _vec(k2)
    @.. broadcast=false veck2=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a31 * k1 + a32 * k2
    f(du, u, p, t + c3 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd3 * dT + (dtC31 * k1 + dtC32 * k2)
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
    f(du, u, p, t + c4 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd4 * dT +
                                         (dtC41 * k1 + dtC42 * k2 + dtC43 * k3)
    else
        @.. broadcast=false du1=dtC41 * k1 + dtC42 * k2 + dtC43 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd4 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck4 = _vec(k4)
    @.. broadcast=false veck4=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4
    f(du, u, p, t + c5 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + dtd5 * dT +
                                         (dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3)
    else
        @.. broadcast=false du1=dtC52 * k2 + dtC54 * k4 + dtC51 * k1 + dtC53 * k3
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + dtd5 * dT + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck5 = _vec(k5)
    @.. broadcast=false veck5=-vecu
    integrator.stats.nsolve += 1

    @.. broadcast=false u=uprev + a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + (dtC61 * k1 + dtC62 * k2 + dtC63 * k3 +
                                          dtC64 * k4 + dtC65 * k5)
    else
        @.. broadcast=false du1=dtC61 * k1 + dtC62 * k2 + dtC63 * k3 + dtC64 * k4 +
                                dtC65 * k5
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck6 = _vec(k6)
    @.. broadcast=false veck6=-vecu
    integrator.stats.nsolve += 1

    u .+= k6
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + (dtC71 * k1 + dtC72 * k2 + dtC73 * k3 +
                                          dtC74 * k4 + dtC75 * k5 + dtC76 * k6)
    else
        @.. broadcast=false du1=dtC71 * k1 + dtC72 * k2 + dtC73 * k3 + dtC74 * k4 +
                                dtC75 * k5 + dtC76 * k6
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck7 = _vec(k7)
    @.. broadcast=false veck7=-vecu
    integrator.stats.nsolve += 1

    u .+= k7
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @.. broadcast=false linsolve_tmp=du + (dtC81 * k1 + dtC82 * k2 + dtC83 * k3 +
                                          dtC84 * k4 + dtC85 * k5 + dtC86 * k6 + dtC87 * k7)
    else
        @.. broadcast=false du1=dtC81 * k1 + dtC82 * k2 + dtC83 * k3 + dtC84 * k4 +
                                dtC85 * k5 + dtC86 * k6 + dtC87 * k7
        mul!(_vec(du2), mass_matrix, _vec(du1))
        @.. broadcast=false linsolve_tmp=du + du2
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    veck8 = _vec(k8)
    @.. broadcast=false veck8=-vecu
    integrator.stats.nsolve += 1

    du .= k8
    u .+= k8

    if integrator.opts.adaptive
	    if (integrator.alg isa Rodas5Pe)
            du = 0.2606326497975715*k1 - 0.005158627295444251*k2 + 1.3038988631109731*k3 + 1.235000722062074*k4 +
               - 0.7931985603795049*k5 - 1.005448461135913*k6 - 0.18044626132120234*k7 + 0.17051519239113755*k8
	    end
        calculate_residuals!(atmp, du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        @unpack h21, h22, h23, h24, h25, h26, h27, h28, h31, h32, h33, h34, h35, h36, h37, h38, h41, h42, h43, h44, h45, h46, h47, h48 = cache.tab
        @.. broadcast=false integrator.k[1]=h21 * k1 + h22 * k2 + h23 * k3 + h24 * k4 +
                                            h25 * k5 + h26 * k6 + h27 * k7 + h28 * k8
        @.. broadcast=false integrator.k[2]=h31 * k1 + h32 * k2 + h33 * k3 + h34 * k4 +
                                            h35 * k5 + h36 * k6 + h37 * k7 + h38 * k8
        @.. broadcast=false integrator.k[3]=h41 * k1 + h42 * k2 + h43 * k3 + h44 * k4 +
                                            h45 * k5 + h46 * k6 + h47 * k7 + h48 * k8
        if (integrator.alg isa Rodas5Pr) && integrator.opts.adaptive && (integrator.EEst < 1.0)
            k2 = 0.5*(uprev + u + 0.5 * (integrator.k[1] + 0.5 * (integrator.k[2] + 0.5 * integrator.k[3])))
            du1 = ( 0.25*(integrator.k[2] + integrator.k[3]) - uprev + u) / dt
            f(du, k2, p, t + dt/2)
            integrator.stats.nf += 1
            if mass_matrix === I
                du2 = du1 - du
            else
                mul!(_vec(du2), mass_matrix, _vec(du1))
                du2 = du2 - du
            end
	        EEst = norm(du2) / (integrator.opts.abstol + integrator.opts.reltol*norm(k2)) 
            integrator.EEst = max(EEst,integrator.EEst)
        end
    end
    cache.linsolve = linres.cache
end

@muladd function perform_step!(integrator, cache::Rosenbrock5Cache{<:Array},
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack du, du1, du2, k1, k2, k3, k4, k5, k6, k7, k8, dT, J, W, uf, tf, linsolve_tmp, jac_config, atmp, weight = cache
    @unpack a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, C21, C31, C32, C41, C42, C43, C51, C52, C53, C54, C61, C62, C63, C64, C65, C71, C72, C73, C74, C75, C76, C81, C82, C83, C84, C85, C86, C87, gamma, d1, d2, d3, d4, d5, c2, c3, c4, c5 = cache.tab

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
    dtC61 = C61 / dt
    dtC62 = C62 / dt
    dtC63 = C63 / dt
    dtC64 = C64 / dt
    dtC65 = C65 / dt
    dtC71 = C71 / dt
    dtC72 = C72 / dt
    dtC73 = C73 / dt
    dtC74 = C74 / dt
    dtC75 = C75 / dt
    dtC76 = C76 / dt
    dtC81 = C81 / dt
    dtC82 = C82 / dt
    dtC83 = C83 / dt
    dtC84 = C84 / dt
    dtC85 = C85 / dt
    dtC86 = C86 / dt
    dtC87 = C87 / dt

    dtd1 = dt * d1
    dtd2 = dt * d2
    dtd3 = dt * d3
    dtd4 = dt * d4
    dtd5 = dt * d5
    dtgamma = dt * gamma

    f(cache.fsalfirst, uprev, p, t) # used in calc_rosenbrock_differentiation!
    integrator.stats.nf += 1

    calc_rosenbrock_differentiation!(integrator, cache, dtd1, dtgamma, repeat_step, true)

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

    vecu = _vec(linres.u)
    veck1 = _vec(k1)

    @inbounds @simd ivdep for i in eachindex(u)
        veck1[i] = -vecu[i]
    end

    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + a21 * k1[i]
    end
    f(du, u, p, t + c2 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd2 * dT[i] + dtC21 * k1[i]
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC21 * k1[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd2 * dT[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k2[i] = -linres.u[i]
    end

    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + a31 * k1[i] + a32 * k2[i]
    end
    f(du, u, p, t + c3 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd3 * dT[i] + (dtC31 * k1[i] + dtC32 * k2[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC31 * k1[i] + dtC32 * k2[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd3 * dT[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k3[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + a41 * k1[i] + a42 * k2[i] + a43 * k3[i]
    end
    f(du, u, p, t + c4 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd4 * dT[i] +
                              (dtC41 * k1[i] + dtC42 * k2[i] + dtC43 * k3[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC41 * k1[i] + dtC42 * k2[i] + dtC43 * k3[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd4 * dT[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k4[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]
    end
    f(du, u, p, t + c5 * dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd5 * dT[i] +
                              (dtC52 * k2[i] + dtC54 * k4[i] + dtC51 * k1[i] +
                               dtC53 * k3[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC52 * k2[i] + dtC54 * k4[i] + dtC51 * k1[i] + dtC53 * k3[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + dtd5 * dT[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k5[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] = uprev[i] + a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] +
               a65 * k5[i]
    end
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + (dtC61 * k1[i] + dtC62 * k2[i] + dtC63 * k3[i] +
                               dtC64 * k4[i] + dtC65 * k5[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC61 * k1[i] + dtC62 * k2[i] + dtC63 * k3[i] + dtC64 * k4[i] +
                     dtC65 * k5[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k6[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] += k6[i]
    end
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + (dtC71 * k1[i] + dtC72 * k2[i] + dtC73 * k3[i] +
                               dtC74 * k4[i] + dtC75 * k5[i] + dtC76 * k6[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC71 * k1[i] + dtC72 * k2[i] + dtC73 * k3[i] + dtC74 * k4[i] +
                     dtC75 * k5[i] + dtC76 * k6[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k7[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] += k7[i]
    end
    f(du, u, p, t + dt)
    integrator.stats.nf += 1

    if mass_matrix === I
        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + (dtC81 * k1[i] + dtC82 * k2[i] + dtC83 * k3[i] +
                               dtC84 * k4[i] + dtC85 * k5[i] + dtC86 * k6[i] +
                               dtC87 * k7[i])
        end
    else
        @inbounds @simd ivdep for i in eachindex(u)
            du1[i] = dtC81 * k1[i] + dtC82 * k2[i] + dtC83 * k3[i] + dtC84 * k4[i] +
                     dtC85 * k5[i] + dtC86 * k6[i] + dtC87 * k7[i]
        end
        mul!(_vec(du2), mass_matrix, _vec(du1))

        @inbounds @simd ivdep for i in eachindex(u)
            linsolve_tmp[i] = du[i] + du2[i]
        end
    end

    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @inbounds @simd ivdep for i in eachindex(u)
        k8[i] = -linres.u[i]
    end
    integrator.stats.nsolve += 1

    @inbounds @simd ivdep for i in eachindex(u)
        u[i] += k8[i]
        du[i] = k8[i]
    end

    if integrator.opts.adaptive
	    if (integrator.alg isa Rodas5Pe)
	    	@inbounds @simd ivdep for i in eachindex(u)
               	du[i] = 0.2606326497975715*k1[i] - 0.005158627295444251*k2[i] + 1.3038988631109731*k3[i] + 1.235000722062074*k4[i] +
                   	- 0.7931985603795049*k5[i] - 1.005448461135913*k6[i] - 0.18044626132120234*k7[i] + 0.17051519239113755*k8[i]
	        end					 
        end
        calculate_residuals!(atmp, du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.opts.calck
        @unpack h21, h22, h23, h24, h25, h26, h27, h28, h31, h32, h33, h34, h35, h36, h37, h38, h41, h42, h43, h44, h45, h46, h47, h48 = cache.tab
        @inbounds @simd ivdep for i in eachindex(u)
            integrator.k[1][i] = h21 * k1[i] + h22 * k2[i] + h23 * k3[i] + h24 * k4[i] +
                                 h25 * k5[i] + h26 * k6[i] + h27 * k7[i] + h28 * k8[i]
            integrator.k[2][i] = h31 * k1[i] + h32 * k2[i] + h33 * k3[i] + h34 * k4[i] +
                                 h35 * k5[i] + h36 * k6[i] + h37 * k7[i] + h38 * k8[i]
            integrator.k[3][i] = h41 * k1[i] + h42 * k2[i] + h43 * k3[i] + h44 * k4[i] +
                                 h45 * k5[i] + h46 * k6[i] + h47 * k7[i] + h48 * k8[i]
	        if (integrator.alg isa Rodas5Pr)
                k2[i] = 0.5*(uprev[i] + u[i] + 0.5 * (integrator.k[1][i] + 0.5 * (integrator.k[2][i] + 0.5 * integrator.k[3][i])))
                du1[i] = ( 0.25*(integrator.k[2][i] + integrator.k[3][i]) - uprev[i] + u[i]) / dt
	        end
        end
        if integrator.opts.adaptive && (integrator.EEst < 1.0) && (integrator.alg isa Rodas5Pr) 
            f(du, k2, p, t + dt/2)
            integrator.stats.nf += 1
            if mass_matrix === I
                @inbounds @simd ivdep for i in eachindex(u)
                    du2[i] = du1[i] - du[i]
                end
            else
                mul!(_vec(du2), mass_matrix, _vec(du1))
                @inbounds @simd ivdep for i in eachindex(u)
                    du2[i] = du2[i] - du[i]
                end
            end
	        EEst = norm(du2) / (integrator.opts.abstol + integrator.opts.reltol*norm(k2)) 
            integrator.EEst = max(EEst,integrator.EEst)
        end
    end
    cache.linsolve = linres.cache
end

@RosenbrockW6S4OS(:init)
@RosenbrockW6S4OS(:performstep)
