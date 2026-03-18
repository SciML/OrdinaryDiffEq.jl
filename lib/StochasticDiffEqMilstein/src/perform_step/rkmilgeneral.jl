@muladd function perform_step!(integrator, cache::RKMilGeneralConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    Jalg = cache.Jalg
    dW = W.dW

    J = get_iterated_I(
        dt, dW, W.dZ, Jalg, integrator.alg.p, integrator.alg.c, alg_order(integrator.alg)
    )

    if SciMLBase.alg_interpretation(integrator.alg) == SciMLBase.AlgorithmInterpretation.Ito
        if dW isa Number || is_diagonal_noise(integrator.sol.prob)
            J = J .- 1 // 2 .* abs(dt)
        else
            J -= 1 // 2 .* UniformScaling(abs(dt))
        end
    end

    du1 = integrator.f(uprev, p, t)
    L = integrator.f.g(uprev, p, t)
    mil_correction = zero(u)
    ggprime_norm = zero(eltype(u)) #0
    # sqdt = integrator.sqdt

    if dW isa Number || is_diagonal_noise(integrator.sol.prob)
        K = @.. uprev + dt * du1
        utilde = (
            SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Ito ? K : uprev
        ) + L * integrator.sqdt
        ggprime = (integrator.f.g(utilde, p, t) .- L) ./ (integrator.sqdt)
        mil_correction = ggprime .* J
        u = K + L .* dW + mil_correction
    else
        for i in 1:length(dW)
            K = uprev + dt * du1 + integrator.sqdt * @view(L[:, i])
            gtmp = integrator.f.g(K, p, t)
            ggprime = @.. (gtmp - L) / integrator.sqdt
            ggprime_norm = zero(eltype(u))
            if integrator.opts.adaptive
                ggprime_norm += integrator.opts.internalnorm(ggprime, t)
            end
            mil_correction += ggprime * @view(J[:, i])
        end
        if integrator.opts.adaptive
            K = @.. uprev + dt * du1
            u = K + L * dW + mil_correction
        else
            u = uprev + dt * du1 + L * dW + mil_correction
        end
    end

    if integrator.opts.adaptive
        du2 = integrator.f(K, p, t + dt)
        if dW isa Number || is_diagonal_noise(integrator.sol.prob)
            tmp = dt * (du2 - du1) / 2
            En = W.dW .^ 3 .* ((du2 - L) / (integrator.sqdt)) .^ 2 / 6
        else
            En = integrator.opts.internalnorm(W.dW, t)^3 * ggprime_norm^2 / 6
            tmp = integrator.opts.internalnorm((@.. dt * (du2 - du1) / 2), t)
        end
        tmp = calculate_residuals(
            tmp, En, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(tmp, t)
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RKMilGeneralCache)
    (; du1, du2, K, tmp, ggprime, L, mil_correction) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    dW = W.dW
    sqdt = integrator.sqdt
    Jalg = cache.Jalg

    get_iterated_I!(
        dt, dW, W.dZ, Jalg, integrator.alg.p, integrator.alg.c, alg_order(integrator.alg)
    )
    J = Jalg.J

    integrator.f(du1, uprev, p, t)
    integrator.f.g(L, uprev, p, t)
    @.. mil_correction = zero(eltype(u))
    ggprime_norm = zero(eltype(ggprime))

    if SciMLBase.alg_interpretation(integrator.alg) == SciMLBase.AlgorithmInterpretation.Ito
        if dW isa Number || is_diagonal_noise(integrator.sol.prob)
            @.. J -= 1 // 2 * abs(dt)
        else
            J -= 1 // 2 .* UniformScaling(abs(dt))
        end
    end

    if dW isa Number || is_diagonal_noise(integrator.sol.prob)
        @.. K = uprev + dt * du1
        @.. du2 = zero(eltype(u))
        tmp .= (
            SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Ito ? K : uprev
        ) .+ integrator.sqdt .* L
        integrator.f.g(du2, tmp, p, t)
        @.. ggprime = (du2 - L) / sqdt
        ggprime_norm = integrator.opts.internalnorm(ggprime, t)
        @.. u = K + L * dW + ggprime * J
    else
        for i in 1:length(dW)
            @.. K = uprev + dt * du1 + sqdt * @view(L[:, i])
            integrator.f.g(ggprime, K, p, t)
            @.. ggprime = (ggprime - L) / sqdt
            if integrator.opts.adaptive
                ggprime_norm += integrator.opts.internalnorm(ggprime, t)
            end
            mul!(tmp, ggprime, @view(J[:, i]))
            @.. mil_correction += tmp
        end
        mul!(tmp, L, dW)
        if integrator.opts.adaptive
            @.. K = uprev + dt * du1
            @.. u = K + tmp + mil_correction
        else
            @.. u = uprev + dt * du1 + tmp + mil_correction
        end
    end

    if integrator.opts.adaptive
        En = integrator.opts.internalnorm(W.dW, t)^3 * ggprime_norm^2 / 6
        integrator.f(du2, K, p, t + dt)
        @.. tmp = integrator.opts.internalnorm(
            integrator.opts.delta * dt * (du2 - du1) /
                2, t
        ) + En

        calculate_residuals!(
            tmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(tmp, t)
    end
    return nothing
end
