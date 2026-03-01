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
    (;
        k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W, tmp, uf, tf,
        linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter!,
    ) = cache
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

    new_W = calc_rosenbrock_differentiation!(integrator, cache, dtγ, dtγ, repeat_step)

    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )

    if repeat_step || !new_W
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
    (;
        k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W, tmp, uf, tf,
        linsolve_tmp, jac_config, atmp, weight, stage_limiter!, step_limiter!,
    ) = cache
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

    new_W = calc_rosenbrock_differentiation!(integrator, cache, dtγ, dtγ, repeat_step)

    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )

    if repeat_step || !new_W
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

    # Time derivative and W (with Jacobian reuse for W-methods)
    dT, W = calc_rosenbrock_differentiation(integrator, cache, dtγ, repeat_step)
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

    # Time derivative and W (with Jacobian reuse for W-methods)
    dT, W = calc_rosenbrock_differentiation(integrator, cache, dtγ, repeat_step)
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

#### Rodas4 type method — Helper functions for type-dispatched update logic
# These helpers avoid value-based branching on tableau_b/tableau_btilde,
# which Julia 1.11 cannot narrow through inference in large functions.

# IIP: Update u and du after stage loop, dispatched on tableau type
@inline function _rosenbrock_update_result!(
        u, du, uprev, ks, A, tab::RodasTableau, alg
    )
    if alg isa Rodas6P
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
end

@inline function _rosenbrock_update_result!(
        u, du, uprev, ks, A, tab::RosenbrockUnifiedTableau, alg
    )
    u .= uprev
    for i in eachindex(ks)
        if !iszero(tab.b[i])
            @.. u += tab.b[i] * ks[i]
        end
    end
    if tab.btilde !== nothing
        du .= 0
        for i in eachindex(ks)
            if !iszero(tab.btilde[i])
                @.. du += tab.btilde[i] * ks[i]
            end
        end
    end
end

# OOP: Update u and du after stage loop, dispatched on tableau type
@inline function _rosenbrock_update_result(
        u, du, uprev, ks, A, tab::RodasTableau, alg, num_stages
    )
    if alg isa Rodas6P
        du_new = ks[16]
        u_new = uprev
        for i in 1:15
            u_new = @.. u_new + A[16, i] * ks[i]
        end
        u_new = u_new .+ ks[16]
        return u_new, du_new
    else
        du_new = ks[num_stages]
        u_new = u .+ ks[num_stages]
        return u_new, du_new
    end
end

@inline function _rosenbrock_update_result(
        u, du, uprev, ks, A, tab::RosenbrockUnifiedTableau, alg, num_stages
    )
    u_new = copy(uprev)
    for i in 1:num_stages
        if !iszero(tab.b[i])
            u_new = @.. u_new + tab.b[i] * ks[i]
        end
    end
    du_new = du
    if tab.btilde !== nothing
        du_new = zero(uprev)
        for i in 1:num_stages
            if !iszero(tab.btilde[i])
                du_new = @.. du_new + tab.btilde[i] * ks[i]
            end
        end
    end
    return u_new, du_new
end

# Whether error estimate should be computed (Rodas encoding always computes error)
@inline _rosenbrock_should_est_error(::RodasTableau) = true
@inline _rosenbrock_should_est_error(tab::RosenbrockUnifiedTableau) = tab.btilde !== nothing

#### Rodas4 type method

function initialize!(integrator, cache::RosenbrockCombinedConstantCache)
    H_rows = size(cache.tab.H, 1)
    integrator.kshortsize = H_rows > 0 ? H_rows : 2
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

    # Time derivative and W (with Jacobian reuse for W-methods)
    dT, W = calc_rosenbrock_differentiation(integrator, cache, dtgamma, repeat_step)
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
    u, du = _rosenbrock_update_result(u, du, uprev, ks, A, cache.tab, integrator.alg, num_stages)

    if integrator.opts.adaptive
        if (integrator.alg isa Rodas5Pe)
            du = 0.2606326497975715 * ks[1] - 0.005158627295444251 * ks[2] +
                1.3038988631109731 * ks[3] + 1.235000722062074 * ks[4] +
                -0.7931985603795049 * ks[5] - 1.005448461135913 * ks[6] -
                0.18044626132120234 * ks[7] + 0.17051519239113755 * ks[8]
        end
        if _rosenbrock_should_est_error(cache.tab)
            atmp = calculate_residuals(
                du, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        end
    end

    if integrator.opts.calck
        if size(H, 1) > 0
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
                        0.5 *
                        (integrator.k[1] + 0.5 * (integrator.k[2] + 0.5 * integrator.k[3]))
                )
                du1 = (0.25 * (integrator.k[2] + integrator.k[3]) - uprev + u) / dt
                du = f(k2, p, t + dt / 2)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
                if mass_matrix === I
                    du2 = du1 - du
                else
                    du2 = mass_matrix * du1 - du
                end
                EEst = norm(du2) /
                    norm(integrator.opts.abstol .+ integrator.opts.reltol .* k2)
                integrator.EEst = max(EEst, integrator.EEst)
            end
        else
            # No H matrix: use basic Hermite interpolation (kshortsize=2)
            integrator.k[1] = f(uprev, p, t)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            integrator.k[2] = f(u, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
    end

    integrator.u = u
    return nothing
end

function initialize!(integrator, cache::RosenbrockCache)
    integrator.kshortsize = length(cache.dense)
    resize!(integrator.k, integrator.kshortsize)
    for i in 1:(integrator.kshortsize)
        integrator.k[i] = cache.dense[i]
    end
    return
end

@muladd function perform_step!(integrator, cache::RosenbrockCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (;
        du, du1, du2, dT, dtC, dtd, J, W, uf, tf, ks, linsolve_tmp,
        jac_config, atmp, weight, stage_limiter!, step_limiter!,
    ) = cache
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

    new_W = calc_rosenbrock_differentiation!(integrator, cache, dtd[1], dtgamma, repeat_step)

    calculate_residuals!(
        weight, fill!(weight, one(eltype(u))), uprev, uprev,
        integrator.opts.abstol, integrator.opts.reltol,
        integrator.opts.internalnorm, t
    )

    if repeat_step || !new_W
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
    _rosenbrock_update_result!(u, du, uprev, ks, A, cache.tab, integrator.alg)

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        if (integrator.alg isa Rodas5Pe)
            @.. du = 0.2606326497975715 * ks[1] - 0.005158627295444251 * ks[2] +
                1.3038988631109731 * ks[3] + 1.235000722062074 * ks[4] +
                -0.7931985603795049 * ks[5] - 1.005448461135913 * ks[6] -
                0.18044626132120234 * ks[7] + 0.17051519239113755 * ks[8]
        end
        if _rosenbrock_should_est_error(cache.tab)
            calculate_residuals!(
                atmp, du, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
            integrator.EEst = integrator.opts.internalnorm(atmp, t)
        end
    end

    if integrator.opts.calck
        if size(H, 1) > 0
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
                        (
                        integrator.k[1] +
                            0.5 * (integrator.k[2] + 0.5 * integrator.k[3])
                    )
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
        else
            # No H matrix: use basic Hermite interpolation (kshortsize=2)
            f(integrator.k[1], uprev, p, t)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            f(integrator.k[2], u, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
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
        _ = calc_rosenbrock_differentiation!(integrator, cache, dt * d[1], dtgamma, repeat_step)

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
