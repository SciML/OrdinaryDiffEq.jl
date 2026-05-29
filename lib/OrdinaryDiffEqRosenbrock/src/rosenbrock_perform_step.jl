#### Rodas4 type method — unified perform_step for all Rosenbrock/Rodas methods

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

# Rosenbrock23/32 (Shampine order-2/3 pair): stage 3 uses the Shampine RHS and the
# error estimate handles the DAE algebraic row.
@muladd function perform_step!(
        integrator,
        cache::RosenbrockCombinedConstantCache{TF, UF, Tab},
        repeat_step = false
    ) where {TF, UF, Tab <: ShampineRosenbrockTableau}
    (; t, dt, uprev, u, f, p) = integrator
    (; A, C, gamma, c, d, H, b, btilde) = cache.tab

    dtC = C ./ dt
    dtd = dt .* d
    dtgamma = dt * gamma
    mass_matrix = integrator.f.mass_matrix

    dT, W = calc_rosenbrock_differentiation(integrator, cache, dtgamma, repeat_step)
    if !issuccess_W(W)
        OrdinaryDiffEqCore.set_EEst!(integrator, 2)
        return nothing
    end

    # Stage 1
    fsalfirst = f(uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    k1 = _reshape(W \ -_vec(@.. fsalfirst + dtd[1] * dT), axes(uprev))
    integrator.stats.nsolve += 1

    # Stage 2 (standard Rosenbrock pattern)
    u = @.. uprev + A[2, 1] * k1
    f2 = f(u, p, t + c[2] * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    coupling2 = mass_matrix === I ? (@.. dtC[2, 1] * k1) : mass_matrix * (@.. dtC[2, 1] * k1)
    k2 = _reshape(W \ -_vec(@.. f2 + dtd[2] * dT + coupling2), axes(uprev))
    integrator.stats.nsolve += 1

    # Stage 3 (Shampine RHS: f₃ + c₃₂·f₂ + 2·f₁ + dt·dT + M·coupling)
    u = @.. uprev + A[3, 1] * k1 + A[3, 2] * k2
    f3 = f(u, p, t + c[3] * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    c32 = -C[3, 2] * gamma
    two = -C[3, 1] * gamma
    raw_coupling = @.. dtC[3, 1] * k1 + dtC[3, 2] * (k1 + k2)
    coupling3 = mass_matrix === I ? raw_coupling : mass_matrix * raw_coupling
    k3 = _reshape(
        W \ -_vec(@.. f3 + c32 * f2 + two * fsalfirst + dt * dT + coupling3), axes(uprev)
    )
    integrator.stats.nsolve += 1

    ks = (k1, k2, k3)
    u = @.. uprev + b[1] * k1 + b[2] * k2 + b[3] * k3

    if integrator.opts.adaptive
        utilde = @.. btilde[1] * k1 + btilde[2] * k2 + btilde[3] * k3
        if mass_matrix !== I && integrator.differential_vars isa AbstractArray &&
                iszero(b[end])
            # Rosenbrock23: b[3]=0 so k3 only appears in btilde. Zero algebraic
            # components to keep the poorly-conditioned algebraic row from inflating EEst.
            utilde = @.. ifelse(integrator.differential_vars, utilde, false)
        end
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        EEst = integrator.opts.internalnorm(atmp, t)
        if mass_matrix !== I && integrator.differential_vars isa AbstractArray
            # Algebraic residual: f(u_stage3)[alg] scaled by inv(abstol) measures drift
            invatol = inv(integrator.opts.abstol)
            atmp = @.. ifelse(integrator.differential_vars, false, f3) * invatol
            EEst += integrator.opts.internalnorm(atmp, t)
        end
        OrdinaryDiffEqCore.set_EEst!(integrator, EEst)
    end

    if integrator.opts.calck
        for j in eachindex(integrator.k)
            integrator.k[j] = zero(integrator.k[1])
            for i in 1:3
                integrator.k[j] = @.. integrator.k[j] + H[j, i] * ks[i]
            end
        end
    end

    integrator.u = u
    return nothing
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

    dT, W = calc_rosenbrock_differentiation(integrator, cache, dtgamma, repeat_step)
    if !issuccess_W(W)
        OrdinaryDiffEqCore.set_EEst!(integrator, 2)
        return nothing
    end

    num_stages = size(A, 1)
    du = f(uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    fsalfirst_cache = du
    linsolve_tmp = @.. du + dtd[1] * dT
    k1 = _reshape(W \ -_vec(linsolve_tmp), axes(uprev))
    # Val(20): fixed size for type stability; must exceed max num_stages across all methods
    ks = ntuple(Returns(k1), Val(20))

    for stage in 2:num_stages
        u = uprev
        for i in 1:(stage - 1)
            u = @.. u + A[stage, i] * ks[i]
        end

        # Some Rodas4 methods have c[stage]==c[stage-1] and identical A rows → reuse du
        if stage > 2 && c[stage] == c[stage - 1] &&
                all(j -> A[stage, j] == A[stage - 1, j], 1:(stage - 1))
        else
            du = f(u, p, t + c[stage] * dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end

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
    tab = cache.tab
    u = copy(uprev)
    for i in 1:num_stages
        if !iszero(tab.b[i])
            u = @.. u + tab.b[i] * ks[i]
        end
    end

    if integrator.opts.adaptive && tab.btilde !== nothing
        # Error estimate using explicit btilde weights
        du = zero(uprev)
        for i in 1:num_stages
            if !iszero(tab.btilde[i])
                du = @.. du + tab.btilde[i] * ks[i]
            end
        end
        atmp = calculate_residuals(
            du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        EEst = integrator.opts.internalnorm(atmp, t)
        OrdinaryDiffEqCore.set_EEst!(integrator, EEst)
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
                    (OrdinaryDiffEqCore.get_EEst(integrator) < 1.0)
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
                OrdinaryDiffEqCore.set_EEst!(integrator, max(EEst, OrdinaryDiffEqCore.get_EEst(integrator)))
            end
        else
            f1 = f(u, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            integrator.k[1] = fsalfirst_cache
            integrator.k[2] = f1
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

# Rosenbrock23/32 (Shampine order-2/3 pair), in-place. See the constant-cache
# method above for the formula; dispatch is by tableau type.
@muladd function perform_step!(
        integrator,
        cache::RosenbrockCache{
            uType, rateType, tabType, uNoUnitsType, JType, WType, TabType,
        },
        repeat_step = false
    ) where {
        uType, rateType, tabType, uNoUnitsType, JType, WType,
        TabType <: ShampineRosenbrockTableau,
    }
    (; t, dt, uprev, u, f, p) = integrator
    (;
        du, du1, du2, dT, dtC, dtd, W, ks, linsolve_tmp,
        atmp, weight, stage_limiter!, step_limiter!,
    ) = cache
    (; A, C, gamma, c, d, H, b, btilde) = cache.tab

    mass_matrix = integrator.f.mass_matrix

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

    # Stage 1 (linsolve_tmp already holds fsalfirst + dtd[1]·dT from calc_tderivative!)
    linres = dolinsolve(
        integrator, cache.linsolve;
        A = (repeat_step || !new_W) ? nothing : W, b = _vec(linsolve_tmp)
    )
    @.. $(_vec(ks[1])) = -linres.u
    integrator.stats.nsolve += 1

    # Stage 2 (standard Rosenbrock pattern)
    u .= uprev
    @.. u += A[2, 1] * ks[1]
    stage_limiter!(u, integrator, p, t + c[2] * dt)
    f(du, u, p, t + c[2] * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    if mass_matrix === I
        @.. du1 = dtC[2, 1] * ks[1]
    else
        @.. du2 = dtC[2, 1] * ks[1]
        mul!(_vec(du1), mass_matrix, _vec(du2))
    end
    @.. linsolve_tmp = du + dtd[2] * dT + du1
    copyto!(cache.fsallast, du)  # stash f₂ for the stage-3 Shampine RHS
    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. $(_vec(ks[2])) = -linres.u
    integrator.stats.nsolve += 1

    # Stage 3 (Shampine RHS: f₃ + c₃₂·f₂ + 2·f₁ + dt·dT + M·coupling)
    u .= uprev
    @.. u += A[3, 1] * ks[1] + A[3, 2] * ks[2]
    stage_limiter!(u, integrator, p, t + c[3] * dt)
    f(du, u, p, t + c[3] * dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    c32 = -C[3, 2] * gamma
    two = -C[3, 1] * gamma
    if mass_matrix === I
        @.. du1 = dtC[3, 1] * ks[1] + dtC[3, 2] * (ks[1] + ks[2])
    else
        @.. du2 = dtC[3, 1] * ks[1] + dtC[3, 2] * (ks[1] + ks[2])
        mul!(_vec(du1), mass_matrix, _vec(du2))
    end
    @.. linsolve_tmp = du + c32 * cache.fsallast + two * cache.fsalfirst + dt * dT + du1
    linres = dolinsolve(integrator, linres.cache; b = _vec(linsolve_tmp))
    @.. $(_vec(ks[3])) = -linres.u
    integrator.stats.nsolve += 1

    # du still holds f(u_stage3) for the DAE algebraic residual below
    u .= uprev
    @.. u += b[1] * ks[1] + b[2] * ks[2] + b[3] * ks[3]
    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        @.. du1 = btilde[1] * ks[1] + btilde[2] * ks[2] + btilde[3] * ks[3]
        if mass_matrix !== I && integrator.differential_vars isa AbstractArray &&
                iszero(b[end])
            # Rosenbrock23: b[3]=0 so ks[3] only appears in btilde. Zero algebraic
            # components to keep the poorly-conditioned algebraic row from inflating EEst.
            dv = integrator.differential_vars
            @.. du1 = ifelse(dv, du1, false)
        end
        calculate_residuals!(
            atmp, du1, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        EEst = integrator.opts.internalnorm(atmp, t)
        if mass_matrix !== I && integrator.differential_vars isa AbstractArray
            # Algebraic residual: f(u_stage3)[alg] scaled by inv(abstol) measures drift
            dv = integrator.differential_vars
            invatol = inv(integrator.opts.abstol)
            @.. atmp = ifelse(dv, false, du) * invatol
            EEst += integrator.opts.internalnorm(atmp, t)
        end
        OrdinaryDiffEqCore.set_EEst!(integrator, EEst)
    end

    if integrator.opts.calck
        for j in eachindex(integrator.k)
            integrator.k[j] .= 0
            for i in 1:3
                @.. integrator.k[j] += H[j, i] * ks[i]
            end
        end
    end

    cache.linsolve = linres.cache
    return nothing
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

    linres = dolinsolve(
        integrator, cache.linsolve;
        A = (repeat_step || !new_W) ? nothing : W, b = _vec(linsolve_tmp)
    )

    @.. $(_vec(ks[1])) = -linres.u
    integrator.stats.nsolve += 1

    for stage in 2:length(ks)
        u .= uprev
        for i in 1:(stage - 1)
            @.. u += A[stage, i] * ks[i]
        end

        stage_limiter!(u, integrator, p, t + c[stage] * dt)
        # Some Rodas4 methods have c[stage]==c[stage-1] and identical A rows → reuse du
        if stage > 2 && c[stage] == c[stage - 1] &&
                all(j -> A[stage, j] == A[stage - 1, j], 1:(stage - 1))
        else
            f(du, u, p, t + c[stage] * dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end

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

    tab = cache.tab
    u .= uprev
    for i in eachindex(ks)
        if !iszero(tab.b[i])
            @.. u += tab.b[i] * ks[i]
        end
    end

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive && tab.btilde !== nothing
        # Error estimate using explicit btilde weights
        du .= 0
        for i in eachindex(ks)
            if !iszero(tab.btilde[i])
                @.. du += tab.btilde[i] * ks[i]
            end
        end
        calculate_residuals!(
            atmp, du, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        EEst = integrator.opts.internalnorm(atmp, t)
        OrdinaryDiffEqCore.set_EEst!(integrator, EEst)
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
                    (OrdinaryDiffEqCore.get_EEst(integrator) < 1.0)
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
                OrdinaryDiffEqCore.set_EEst!(integrator, max(EEst, OrdinaryDiffEqCore.get_EEst(integrator)))
            end
        else
            f(du, u, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            @.. integrator.k[1] = cache.fsalfirst
            @.. integrator.k[2] = du
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
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
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
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
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
