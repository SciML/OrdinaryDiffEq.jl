################################################################################
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
            atmp = ifelse(integrator.differential_vars, false, integrator.fsallast) .* invatol
            integrator.EEst += integrator.opts.internalnorm(atmp, t)
        end
    end

    integrator.k[1] = k₁
    integrator.k[2] = k₂
    integrator.u = u
    return nothing
end

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

#### Arbirary order method
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
    ks = ntuple(Returns(k1), 19)
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
      u = @.. u + ks[16]
    else
      du = ks[end]
      u = @.. u + ks[end]
    end

    if integrator.opts.adaptive
        atmp = calculate_residuals(du, uprev, u, integrator.opts.abstol,
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
    utilde = du

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
    num_stages = length(ks)
    for stage in 2:num_stages
        u .= uprev
        for i in 1:(stage - 1)
            @.. u += A[stage, i] * ks[i]
        end

        stage_limiter!(u, integrator, p, t + c[stage] * dt)
        f(du, u, p, t + c[stage] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

        du1 .= 0
        for i in 1:(stage - 1)
            @.. du1 += dtC[stage, i] * ks[i]
        end
        if mass_matrix !== I
            mul!(_vec(du2), mass_matrix, _vec(du1))
            du1 .= du2
        end
        @.. linsolve_tmp = du + dtd[stage] * dT + du1

        linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
        @.. $(_vec(ks[stage])) = -linres.u
        integrator.stats.nsolve += 1
    end
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
        @.. utilde = 0 * u
        for i in 1:num_stages
            @.. utilde += btilde[i] * ks[i]
        end
        if (integrator.alg isa Rodas5Pe)
            @.. du = 0.2606326497975715 * ks[1] - 0.005158627295444251 * ks[2] +
                     1.3038988631109731 * ks[3] + 1.235000722062074 * ks[4] +
                     -0.7931985603795049 * ks[5] - 1.005448461135913 * ks[6] -
                     0.18044626132120234 * ks[7] + 0.17051519239113755 * ks[8]
        end
        calculate_residuals!(atmp, du, uprev, u, integrator.opts.abstol,
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
