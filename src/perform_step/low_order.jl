@muladd function perform_step!(integrator, cache::EMConstantCache)
    (; t, dt, uprev, u, W, P, c, p, f) = integrator

    K = uprev .+ dt .* integrator.f(uprev, p, t)

    if is_split_step(integrator.alg)
        u_choice = K
    else
        u_choice = uprev
    end

    if !is_diagonal_noise(integrator.sol.prob) || W.dW isa Number
        noise = integrator.f.g(u_choice, p, t) * W.dW
    else
        noise = integrator.f.g(u_choice, p, t) .* W.dW
    end

    if P !== nothing
        tmp = c(uprev, p, t, P.dW, nothing)
        u = K + noise + tmp
    else
        u = K + noise
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::EMCache)
    (; tmp, rtmp1, rtmp2) = cache
    (; t, dt, uprev, u, W, P, c, p) = integrator
    integrator.f(rtmp1, uprev, p, t)

    @.. u = uprev + dt * rtmp1

    if is_split_step(integrator.alg)
        u_choice = u
    else
        u_choice = uprev
    end

    integrator.f.g(rtmp2, u_choice, p, t)

    if P !== nothing
        c(tmp, uprev, p, t, P.dW, nothing)
    end

    if is_diagonal_noise(integrator.sol.prob)
        @.. rtmp2 *= W.dW
        if P !== nothing
            @.. u += rtmp2 + tmp
        else
            @.. u += rtmp2
        end
    else
        mul!(rtmp1, rtmp2, W.dW)
        if P !== nothing
            @.. u += rtmp1 + tmp
        else
            @.. u += rtmp1
        end
    end
end

@muladd function perform_step!(integrator, cache::EulerHeunConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    ftmp = integrator.f(uprev, p, t)
    gtmp = integrator.f.g(uprev, p, t)
    if is_diagonal_noise(integrator.sol.prob)
        noise = gtmp .* W.dW
    else
        noise = gtmp * W.dW
    end
    tmp = @.. uprev + dt * ftmp + noise
    gtmp2 = (gtmp .+ integrator.f.g(tmp, p, t + dt)) ./ 2
    if is_diagonal_noise(integrator.sol.prob)
        noise2 = gtmp2 .* W.dW
    else
        noise2 = gtmp2 * W.dW
    end
    u = uprev .+ (dt / 2) .* (ftmp .+ integrator.f(tmp, p, t + dt)) .+ noise2
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::EulerHeunCache)
    (; ftmp1, ftmp2, gtmp1, gtmp2, tmp, nrtmp) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    integrator.f(ftmp1, uprev, p, t)
    integrator.f.g(gtmp1, uprev, p, t)

    if is_diagonal_noise(integrator.sol.prob)
        @.. nrtmp = gtmp1 * W.dW
    else
        mul!(nrtmp, gtmp1, W.dW)
    end

    @.. tmp = uprev + dt * ftmp1 + nrtmp

    integrator.f(ftmp2, tmp, p, t + dt)
    integrator.f.g(gtmp2, tmp, p, t + dt)

    if is_diagonal_noise(integrator.sol.prob)
        @.. nrtmp = W.dW * (gtmp1 + gtmp2) / 2
    else
        # nrtmp already contains gtmp1 * W.dW from stage 1.
        # By linearity: 0.5*(gtmp1+gtmp2)*W.dW == 0.5*(gtmp2*W.dW) + 0.5*(gtmp1*W.dW).
        # Avoid forming (gtmp1 + gtmp2), which would allocate a temporary SparseMatrixCSC.
        # Use 5-arg mul! to accumulate directly into the cached vector (allocation-free).
        mul!(nrtmp, gtmp2, W.dW, convert(eltype(nrtmp), 0.5), convert(eltype(nrtmp), 0.5))
    end

    dto2 = dt / 2
    @.. u = uprev + dto2 * (ftmp1 + ftmp2) + nrtmp
end

@muladd function perform_step!(integrator, cache::RandomEMConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    u = uprev .+ dt .* integrator.f(uprev, p, t, W.curW)
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomEMCache)
    (; rtmp) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    integrator.f(rtmp, uprev, p, t, W.curW)
    @.. u = uprev + dt * rtmp
end

@muladd function perform_step!(integrator, cache::RandomTamedEMConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    ftmp = integrator.f(uprev, p, t, W.curW)
    u = uprev .+ dt .* ftmp ./ (1 .+ dt .* norm(ftmp))
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomTamedEMCache)
    (; rtmp) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    integrator.f(rtmp, uprev, p, t, W.curW)
    @.. u = uprev + dt * rtmp / (1 + dt * norm(rtmp))
end

@muladd function perform_step!(integrator, cache::RandomHeunConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    ftmp = integrator.f(uprev, p, t, W.curW)
    tmp = @.. uprev + dt * ftmp
    wtmp = @.. W.curW + W.dW
    u = uprev .+ (dt / 2) .* (ftmp .+ integrator.f(tmp, p, t + dt, wtmp))
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RandomHeunCache)
    (; tmp, rtmp1, rtmp2, wtmp) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    integrator.f(rtmp1, uprev, p, t, W.curW)
    @.. tmp = uprev + dt * rtmp1
    if W.dW isa Number
        wtmp = W.curW + W.dW
    else
        @.. wtmp = W.curW + W.dW
    end
    integrator.f(rtmp2, tmp, p, t + dt, wtmp)
    @.. u = uprev + (dt / 2) * (rtmp1 + rtmp2)
end

# weak approximation EM
@muladd function perform_step!(integrator, cache::SimplifiedEMConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator

    K = uprev .+ dt .* integrator.f(uprev, p, t)

    _dW = map(x -> calc_twopoint_random(integrator.sqdt, x), W.dW)

    if !is_diagonal_noise(integrator.sol.prob) || W.dW isa Number
        noise = integrator.f.g(uprev, p, t) * _dW
    else
        noise = integrator.f.g(uprev, p, t) .* _dW
    end

    u = K + noise
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SimplifiedEMCache)
    (; rtmp1, rtmp2, _dW) = cache
    (; t, dt, uprev, u, W, p, f) = integrator

    integrator.f(rtmp1, uprev, p, t)

    @.. u = uprev + dt * rtmp1

    integrator.f.g(rtmp2, uprev, p, t)

    if W.dW isa Union{SArray, Number}
        _dW = map(x -> calc_twopoint_random(integrator.sqdt, x), W.dW)
    else
        calc_twopoint_random!(_dW, integrator.sqdt, W.dW)
    end

    if is_diagonal_noise(integrator.sol.prob)
        @.. rtmp2 *= _dW
        @.. u += rtmp2
    else
        mul!(rtmp1, rtmp2, _dW)
        @.. u += rtmp1
    end
end

@muladd function perform_step!(integrator, cache::RKMilConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    du1 = integrator.f(uprev, p, t)
    K = @.. uprev + dt * du1
    L = integrator.f.g(uprev, p, t)
    mil_correction = zero(u)
    if SciMLBase.alg_interpretation(integrator.alg) == SciMLBase.AlgorithmInterpretation.Ito
        utilde = K + L * integrator.sqdt
        ggprime = (integrator.f.g(utilde, p, t) .- L) ./ (integrator.sqdt)
        mil_correction = ggprime .* (W.dW .^ 2 .- abs(dt)) ./ 2
    elseif SciMLBase.alg_interpretation(integrator.alg) ==
            SciMLBase.AlgorithmInterpretation.Stratonovich
        utilde = uprev + L * integrator.sqdt
        ggprime = (integrator.f.g(utilde, p, t) .- L) ./ (integrator.sqdt)
        mil_correction = ggprime .* (W.dW .^ 2) ./ 2
    else
        error("Algorithm interpretation invalid. Use either SciMLBase.AlgorithmInterpretation.Ito or SciMLBase.AlgorithmInterpretation.Stratonovich")
    end
    u = K + L .* W.dW + mil_correction

    if integrator.opts.adaptive
        du2 = integrator.f(K, p, t + dt)
        Ed = dt * (du2 - du1) / 2
        En = W.dW .^ 3 .* ((du2 - L) / (integrator.sqdt)) .^ 2 / 6

        resids = calculate_residuals(
            Ed, En, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.delta,
            integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(resids, t)
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RKMilCache)
    (; du1, du2, K, tmp, L) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    integrator.f(du1, uprev, p, t)
    integrator.f.g(L, uprev, p, t)

    @.. K = uprev + dt * du1
    @.. du2 = zero(eltype(u)) # This makes it safe to re-use the array
    if SciMLBase.alg_interpretation(integrator.alg) == SciMLBase.AlgorithmInterpretation.Ito
        @.. tmp = K + integrator.sqdt * L
        integrator.f.g(du2, tmp, p, t)
        @.. tmp = (du2 - L) / (2integrator.sqdt) * (W.dW .^ 2 - abs(dt))
    elseif SciMLBase.alg_interpretation(integrator.alg) ==
            SciMLBase.AlgorithmInterpretation.Stratonovich
        @.. tmp = uprev + integrator.sqdt * L
        integrator.f.g(du2, tmp, p, t)
        @.. tmp = (du2 - L) / (2integrator.sqdt) * (W.dW .^ 2)
    else
        error("Algorithm interpretation invalid. Use either SciMLBase.AlgorithmInterpretation.Ito or SciMLBase.AlgorithmInterpretation.Stratonovich")
    end
    @.. u = K + L * W.dW + tmp
    if integrator.opts.adaptive
        @.. tmp = integrator.opts.internalnorm(W.dW^3, t) *
            integrator.opts.internalnorm((du2 - L) / (integrator.sqdt), t)^2 / 6
        integrator.f(du2, K, p, t + dt)
        @.. tmp += integrator.opts.internalnorm(integrator.opts.delta * dt * (du2 - du1) / 2, t)

        calculate_residuals!(
            tmp, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(tmp, t)
    end
end

@muladd function perform_step!(integrator, cache::RKMilCommuteConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    dW = W.dW
    sqdt = integrator.sqdt
    Jalg = cache.Jalg

    ggprime_norm = 0.0

    J = get_iterated_I(dt, dW, W.dZ, Jalg)

    mil_correction = zero(u)
    if SciMLBase.alg_interpretation(integrator.alg) == SciMLBase.AlgorithmInterpretation.Ito
        if dW isa Number || is_diagonal_noise(integrator.sol.prob)
            J = J .- 1 // 2 .* abs(dt)
        else
            J -= 1 // 2 .* UniformScaling(abs(dt))
        end
    end

    du1 = integrator.f(uprev, p, t)
    L = integrator.f.g(uprev, p, t)

    K = uprev + dt * du1

    if is_diagonal_noise(integrator.sol.prob)
        tmp = (
            SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Ito ? K : uprev
        ) .+ integrator.sqdt .* L
        gtmp = integrator.f.g(tmp, p, t)
        Dgj = (gtmp - L) / sqdt
        ggprime_norm = integrator.opts.internalnorm(Dgj, t)
        u = @.. K + L * dW + Dgj * J
    else
        for j in 1:length(dW)
            if dW isa Number
                Kj = K + sqdt * L
            else
                Kj = K + sqdt * @view(L[:, j])
            end
            gtmp = integrator.f.g(Kj, p, t)
            Dgj = (gtmp - L) / sqdt
            if integrator.opts.adaptive
                ggprime_norm += integrator.opts.internalnorm(Dgj, t)
            end
            if dW isa Number
                tmp = Dgj * J
            else
                tmp = Dgj * @view(J[:, j])
            end
            mil_correction += tmp
        end
        tmp = L * dW
        u = uprev + dt * du1 + tmp + mil_correction
    end
    if integrator.opts.adaptive
        En = integrator.opts.internalnorm(dW, t)^3 * ggprime_norm^2 / 6
        du2 = integrator.f(K, p, t + dt)
        tmp = integrator.opts.internalnorm(integrator.opts.delta * dt * (du2 - du1) / 2, t) + En

        tmp = calculate_residuals(
            tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(tmp, t)
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::RKMilCommuteCache)
    (; du1, du2, K, gtmp, L) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    (; mil_correction, Kj, Dgj, tmp) = cache
    dW = W.dW
    sqdt = integrator.sqdt

    ggprime_norm = 0.0

    Jalg = cache.Jalg

    get_iterated_I!(dt, dW, W.dZ, Jalg)
    J = Jalg.J

    @.. mil_correction = zero(u)
    if SciMLBase.alg_interpretation(integrator.alg) == SciMLBase.AlgorithmInterpretation.Ito
        if dW isa Number || is_diagonal_noise(integrator.sol.prob)
            @.. J -= 1 // 2 * abs(dt)
        else
            J -= 1 // 2 .* UniformScaling(abs(dt))
        end
    end

    integrator.f(du1, uprev, p, t)
    integrator.f.g(L, uprev, p, t)

    @.. K = uprev + dt * du1

    if is_diagonal_noise(integrator.sol.prob)
        tmp .= (
            SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Ito ? K : uprev
        ) .+ integrator.sqdt .* L
        integrator.f.g(gtmp, tmp, p, t)
        @.. Dgj = (gtmp - L) / sqdt
        ggprime_norm = integrator.opts.internalnorm(Dgj, t)
        @.. u = K + L * dW + Dgj * J
    else
        for j in 1:length(dW)
            @.. Kj = K + sqdt * @view(L[:, j]) # This works too
            #Kj .= uprev .+ sqdt*L[:,j]
            integrator.f.g(gtmp, Kj, p, t)
            @.. Dgj = (gtmp - L) / sqdt
            if integrator.opts.adaptive
                ggprime_norm += integrator.opts.internalnorm(Dgj, t)
            end
            mul!(tmp, Dgj, @view(J[:, j]))
            mil_correction .+= tmp
        end
        mul!(tmp, L, dW)
        @.. u .= uprev + dt * du1 + tmp + mil_correction
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
end

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

    du₁ = integrator.f(uprev, p, t)
    L = integrator.f.g(uprev, p, t)
    mil_correction = zero(u)
    ggprime_norm = zero(eltype(u)) #0
    # sqdt = integrator.sqdt

    if dW isa Number || is_diagonal_noise(integrator.sol.prob)
        K = @.. uprev + dt * du₁
        utilde = (
            SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Ito ? K : uprev
        ) + L * integrator.sqdt
        ggprime = (integrator.f.g(utilde, p, t) .- L) ./ (integrator.sqdt)
        mil_correction = ggprime .* J
        u = K + L .* dW + mil_correction
    else
        for i in 1:length(dW)
            K = uprev + dt * du₁ + integrator.sqdt * @view(L[:, i])
            gtmp = integrator.f.g(K, p, t)
            ggprime = @.. (gtmp - L) / integrator.sqdt
            ggprime_norm = zero(eltype(u))
            if integrator.opts.adaptive
                ggprime_norm += integrator.opts.internalnorm(ggprime, t)
            end
            mil_correction += ggprime * @view(J[:, i])
        end
        if integrator.opts.adaptive
            K = @.. uprev + dt * du₁
            u = K + L * dW + mil_correction
        else
            u = uprev + dt * du₁ + L * dW + mil_correction
        end
    end

    if integrator.opts.adaptive
        du₂ = integrator.f(K, p, t + dt)
        if dW isa Number || is_diagonal_noise(integrator.sol.prob)
            tmp = dt * (du₂ - du₁) / 2
            En = W.dW .^ 3 .* ((du₂ - L) / (integrator.sqdt)) .^ 2 / 6
        else
            En = integrator.opts.internalnorm(W.dW, t)^3 * ggprime_norm^2 / 6
            tmp = integrator.opts.internalnorm((@.. dt * (du₂ - du₁) / 2), t)
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
    (; du₁, du₂, K, tmp, ggprime, L, mil_correction) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    dW = W.dW
    sqdt = integrator.sqdt
    Jalg = cache.Jalg

    get_iterated_I!(
        dt, dW, W.dZ, Jalg, integrator.alg.p, integrator.alg.c, alg_order(integrator.alg)
    )
    J = Jalg.J

    integrator.f(du₁, uprev, p, t)
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
        @.. K = uprev + dt * du₁
        @.. du₂ = zero(eltype(u))
        tmp .= (
            SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Ito ? K : uprev
        ) .+ integrator.sqdt .* L
        integrator.f.g(du₂, tmp, p, t)
        @.. ggprime = (du₂ - L) / sqdt
        ggprime_norm = integrator.opts.internalnorm(ggprime, t)
        @.. u = K + L * dW + ggprime * J
    else
        for i in 1:length(dW)
            @.. K = uprev + dt * du₁ + sqdt * @view(L[:, i])
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
            @.. K = uprev + dt * du₁
            @.. u = K + tmp + mil_correction
        else
            @.. u = uprev + dt * du₁ + tmp + mil_correction
        end
    end

    if integrator.opts.adaptive
        En = integrator.opts.internalnorm(W.dW, t)^3 * ggprime_norm^2 / 6
        integrator.f(du₂, K, p, t + dt)
        @.. tmp = integrator.opts.internalnorm(
            integrator.opts.delta * dt * (du₂ - du₁) /
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
