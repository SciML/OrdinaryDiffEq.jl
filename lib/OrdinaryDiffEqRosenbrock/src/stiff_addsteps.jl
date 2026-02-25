function _ode_addsteps!(
        k, t, uprev, u, dt, f, p,
        cache::Union{
            Rosenbrock23ConstantCache,
            Rosenbrock32ConstantCache,
        },
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    if length(k) < 2 || always_calc_begin
        (; tf, uf, d) = cache
        dtγ = dt * d
        neginvdtγ = -inv(dtγ)
        dto2 = dt / 2
        tf.u = uprev

        autodiff_alg = cache.autodiff

        if autodiff_alg isa AutoFiniteDiff
            autodiff_alg = SciMLBase.@set autodiff_alg.dir = sign(dt)
        end

        dT = DI.derivative(tf, autodiff_alg, t)

        mass_matrix = f.mass_matrix
        if uprev isa Number
            J = DI.derivative(uf, autodiff_alg, uprev)
            W = neginvdtγ .+ J
        else
            J = DI.jacobian(uf, autodiff_alg, uprev)
            if mass_matrix isa UniformScaling
                W = neginvdtγ * mass_matrix + J
            else
                W = @.. neginvdtγ * mass_matrix .+ J
            end
        end
        f₀ = f(uprev, p, t)
        k₁ = _reshape(W \ _vec((f₀ + dtγ * dT)), axes(uprev)) * neginvdtγ
        tmp = @.. uprev + dto2 * k₁
        f₁ = f(tmp, p, t + dto2)
        if mass_matrix === I
            k₂ = _reshape(W \ _vec(f₁ - k₁), axes(uprev))
        else
            k₂ = _reshape(W \ _vec(f₁ - mass_matrix * k₁), axes(uprev))
        end
        k₂ = @.. k₂ * neginvdtγ + k₁
        copyat_or_push!(k, 1, k₁)
        copyat_or_push!(k, 2, k₂)
    end
    return nothing
end

function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::RosenbrockCombinedConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    if length(k) < size(cache.tab.H, 1) || always_calc_begin
        (; tf, uf) = cache
        (; A, C, gamma, c, d, H) = cache.tab

        # Precalculations
        dtC = C ./ dt
        dtd = dt .* d
        dtgamma = dt * gamma
        mass_matrix = f.mass_matrix

        # Time derivative
        tf.u = uprev

        autodiff_alg = cache.autodiff

        if autodiff_alg isa AutoFiniteDiff
            autodiff_alg = SciMLBase.@set autodiff_alg.dir = sign(dt)
        end

        dT = DI.derivative(tf, autodiff_alg, t)

        # Jacobian
        uf.t = t
        if uprev isa AbstractArray
            J = DI.jacobian(uf, autodiff_alg, uprev)
            W = mass_matrix / dtgamma - J
        else
            J = DI.derivative(uf, autodiff_alg, uprev)
            W = 1 / dtgamma - J
        end

        num_stages = size(A, 1)
        du = f(u, p, t)
        linsolve_tmp = @.. du + dtd[1] * dT
        k1 = _reshape(W \ _vec(linsolve_tmp), axes(uprev))
        # constant number for type stability make sure this is greater than num_stages
        ks = ntuple(Returns(k1), Val(20))
        # Last stage affect's ks for Rodas5,5P,6P
        for stage in 2:num_stages
            u = uprev
            for i in 1:(stage - 1)
                u = @.. u + A[stage, i] * ks[i]
            end

            du = f(u, p, t + c[stage] * dt)

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
            ks = Base.setindex(ks, _reshape(W \ _vec(linsolve_tmp), axes(uprev)), stage)
        end

        for j in 1:size(H, 1)
            kj = zero(ks[1])
            # Last stage affect's ks for Rodas5,5P,6P
            for i in 1:num_stages
                kj = @.. kj + H[j, i] * ks[i]
            end
            copyat_or_push!(k, j, kj)
        end
    end
    return nothing
end

function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::RosenbrockCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    if length(k) < 2 || always_calc_begin
        (; du, du1, du2, tmp, ks, dT, J, W, uf, tf, linsolve_tmp, jac_config, fsalfirst, weight) = cache
        (; A, C, gamma, c, d, H) = cache.tab

        # Assignments
        sizeu = size(u)
        uidx = eachindex(uprev)
        mass_matrix = f.mass_matrix
        tmp = ks[end]

        # Precalculations
        dtC = C ./ dt
        dtd = dt .* d
        dtgamma = dt * gamma

        @.. linsolve_tmp = @muladd fsalfirst + dtgamma * dT

        # Jacobian does not need to be re-evaluated after an event since it's unchanged
        jacobian2W!(W, mass_matrix, dtgamma, J)

        linsolve = cache.linsolve

        linres = dolinsolve(
            cache, linsolve; A = W, b = _vec(linsolve_tmp), reltol = cache.reltol
        )
        @.. $(_vec(ks[1])) = -linres.u
        # Last stage affect's ks for Rodas5,5P,6P
        for stage in 2:length(ks)
            tmp .= uprev
            for i in 1:(stage - 1)
                @.. tmp += A[stage, i] * _vec(ks[i])
            end
            f(du, tmp, p, t + c[stage] * dt)

            if mass_matrix === I
                @.. linsolve_tmp = du + dtd[stage] * dT
                for i in 1:(stage - 1)
                    @.. linsolve_tmp += dtC[stage, i] * _vec(ks[i])
                end
            else
                du1 .= du
                for i in 1:(stage - 1)
                    @.. du1 += dtC[stage, i] * _vec(ks[i])
                end
                mul!(_vec(du2), mass_matrix, _vec(du1))
                @.. linsolve_tmp = du + dtd[stage] * dT + du2
            end

            linres = dolinsolve(
                cache, linres.cache; b = _vec(linsolve_tmp), reltol = cache.reltol
            )
            @.. $(_vec(ks[stage])) = -linres.u
        end

        for j in 1:size(H, 1)
            copyat_or_push!(k, j, zero(du))
            # Last stage affect's ks for Rodas5,5P,6P
            for i in 1:length(ks)
                @.. k[j] += H[j, i] * _vec(ks[i])
            end
        end
    end
    return nothing
end

# Tsit5DA: explicit RK stages for addsteps (pure ODE only, ignores algebraic coupling)
function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::HybridExplicitImplicitConstantCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    if length(k) < size(cache.tab.H, 1) || always_calc_begin
        (; tf, uf) = cache
        (; A, C, gamma, b, c, d, H) = cache.tab

        num_stages = size(A, 1)

        # Compute explicit RK stages
        du = f(uprev, p, t)
        ks1 = dt .* du
        ks = ntuple(Returns(ks1), Val(12))
        for stage in 2:num_stages
            u_stage = uprev
            for i in 1:(stage - 1)
                u_stage = @.. u_stage + A[stage, i] * ks[i]
            end
            du = f(u_stage, p, t + c[stage] * dt)
            ks = Base.setindex(ks, dt .* du, stage)
        end

        # Dense output
        for j in 1:size(H, 1)
            kj = zero(ks[1])
            for i in 1:num_stages
                kj = @.. kj + H[j, i] * ks[i]
            end
            copyat_or_push!(k, j, kj)
        end
    end
    return nothing
end

function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::HybridExplicitImplicitCache,
        always_calc_begin = false, allow_calc_end = true,
        force_calc_end = false
    )
    if length(k) < size(cache.tab.H, 1) || always_calc_begin
        (; du, du1, du2, tmp, ks, dT, J, W, uf, tf, linsolve_tmp, jac_config, fsalfirst, weight) = cache
        (; A, C, gamma, b, c, d, H) = cache.tab

        num_stages = size(A, 1)
        mass_matrix = f.mass_matrix

        # Compute explicit RK stages
        f(fsalfirst, uprev, p, t)
        @.. ks[1] = dt * fsalfirst

        for stage in 2:num_stages
            tmp .= uprev
            for i in 1:(stage - 1)
                @.. tmp += A[stage, i] * ks[i]
            end
            f(du, tmp, p, t + c[stage] * dt)
            @.. ks[stage] = dt * du
        end

        # Dense output
        for j in 1:size(H, 1)
            copyat_or_push!(k, j, zero(du))
            for i in 1:num_stages
                @.. k[j] += H[j, i] * _vec(ks[i])
            end
        end
    end
    return nothing
end
