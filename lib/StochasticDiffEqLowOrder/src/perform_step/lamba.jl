@muladd function perform_step!(integrator, cache::LambaEMConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    du1 = integrator.f(uprev, p, t)
    K = @.. uprev + dt * du1

    if is_split_step(integrator.alg)
        L = integrator.f.g(uprev, p, t + dt)
    else
        L = integrator.f.g(K, p, t + dt)
    end

    mil_correction = zero(u)
    if !is_diagonal_noise(integrator.sol.prob) || W.dW isa Number
        noise = L * W.dW
    else
        noise = L .* W.dW
    end

    u = K + noise

    if integrator.opts.adaptive
        du2 = integrator.f(K, p, t + dt)
        Ed = dt * (du2 - du1) / 2

        utilde = K + noise * integrator.sqdt #L*integrator.sqdt
        ggprime = (integrator.f.g(utilde, p, t) .- L) ./ (integrator.sqdt)
        if !is_diagonal_noise(integrator.sol.prob) || W.dW isa Number
            En = ggprime * (W.dW .^ 2 .- dt) ./ 2
        else
            En = ggprime .* (W.dW .^ 2 .- dt) ./ 2
        end
        resids = calculate_residuals(
            Ed, En, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(resids, t))
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::LambaEMCache)
    (; du1, du2, K, tmp, L, gtmp, dW_cache) = cache
    (; t, dt, uprev, u, W, p, f) = integrator

    integrator.f(du1, uprev, p, t)
    @.. K = uprev + dt * du1

    if is_split_step(integrator.alg)
        integrator.f.g(L, K, p, t + dt)
    else
        integrator.f.g(L, uprev, p, t + dt)
    end

    if is_diagonal_noise(integrator.sol.prob)
        @.. tmp = L * W.dW
    else
        mul!(tmp, L, W.dW)
    end

    @.. u = K + tmp

    if integrator.opts.adaptive
        if !is_diagonal_noise(integrator.sol.prob)
            g_sized = norm(L, 2)
        else
            g_sized = L
        end

        if !is_diagonal_noise(integrator.sol.prob)
            @.. tmp = K + integrator.sqdt * g_sized
            integrator.f.g(gtmp, tmp, p, t)
            g_sized2 = norm(gtmp, 2)
            @.. dW_cache = W.dW .^ 2 - dt
            diff_tmp = integrator.opts.internalnorm(dW_cache, t)
            En = (g_sized2 - g_sized) / (2integrator.sqdt) * diff_tmp
            @.. tmp = En
        else
            @.. tmp = K + integrator.sqdt * L
            integrator.f.g(gtmp, tmp, p, t)
            @.. tmp = (gtmp - L) / (2integrator.sqdt) * (W.dW .^ 2 - dt)
        end

        # Ed
        integrator.f(du2, K, p, t + dt)
        @.. tmp += integrator.opts.internalnorm(
            integrator.opts.delta * dt * (du2 - du1) /
                2, t
        )

        calculate_residuals!(
            tmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(tmp, t))
    end
end

@muladd function perform_step!(integrator, cache::LambaEulerHeunConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    du1 = integrator.f(uprev, p, t)
    K = uprev + dt * du1
    L = integrator.f.g(uprev, p, t)

    if is_diagonal_noise(integrator.sol.prob)
        noise = L .* W.dW
    else
        noise = L * W.dW
    end
    tmp = K .+ noise
    gtmp2 = (1 / 2) .* (L .+ integrator.f.g(tmp, p, t + dt))
    if is_diagonal_noise(integrator.sol.prob)
        noise2 = gtmp2 .* W.dW
    else
        noise2 = gtmp2 * W.dW
    end

    u = uprev .+ (dt / 2) .* (du1 .+ integrator.f(tmp, p, t + dt)) .+ noise2

    if integrator.opts.adaptive
        du2 = integrator.f(K, p, t + dt)
        Ed = dt * (du2 - du1) / 2

        utilde = uprev + L * integrator.sqdt
        ggprime = (integrator.f.g(utilde, p, t) .- L) ./ (integrator.sqdt)
        En = ggprime .* (W.dW .^ 2) ./ 2

        resids = calculate_residuals(
            Ed, En, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(resids, t))
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::LambaEulerHeunCache)
    (; du1, du2, K, tmp, L, gtmp, dW_cache) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    integrator.f(du1, uprev, p, t)
    integrator.f.g(L, uprev, p, t)
    @.. K = uprev + dt * du1

    if is_diagonal_noise(integrator.sol.prob)
        @.. tmp = L * W.dW
    else
        mul!(tmp, L, W.dW)
    end

    @.. tmp = K + tmp

    integrator.f(du2, tmp, p, t + dt)
    integrator.f.g(gtmp, tmp, p, t + dt)

    if is_diagonal_noise(integrator.sol.prob)
        @.. tmp = (1 / 2) * W.dW * (L + gtmp)
    else
        @.. gtmp = (1 / 2) * (L + gtmp)
        mul!(tmp, gtmp, W.dW)
    end

    dto2 = dt * (1 / 2)
    @.. u = uprev + dto2 * (du1 + du2) + tmp

    if integrator.opts.adaptive
        if !is_diagonal_noise(integrator.sol.prob)
            g_sized = norm(L, 2)
        else
            g_sized = L
        end

        if !is_diagonal_noise(integrator.sol.prob)
            @.. tmp = uprev + integrator.sqdt * g_sized
            integrator.f.g(gtmp, tmp, p, t)
            g_sized2 = norm(gtmp, 2)
            @.. dW_cache = W.dW .^ 2
            diff_tmp = integrator.opts.internalnorm(dW_cache, t)
            En = (g_sized2 - g_sized) / (2integrator.sqdt) * diff_tmp
            @.. tmp = En
        else
            @.. tmp = uprev + integrator.sqdt * L
            integrator.f.g(gtmp, tmp, p, t)
            @.. tmp = (gtmp - L) / (2integrator.sqdt) * (W.dW .^ 2)
        end

        # Ed
        integrator.f(du2, K, p, t + dt)
        @.. tmp += integrator.opts.internalnorm(
            integrator.opts.delta * dt * (du2 - du1) /
                2, t
        )

        calculate_residuals!(
            tmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(tmp, t))
    end
end
