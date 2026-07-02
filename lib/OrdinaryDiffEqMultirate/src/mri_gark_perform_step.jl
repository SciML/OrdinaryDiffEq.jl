function initialize!(integrator, cache::MRIGARKCache)
    integrator.kshortsize = 2
    (; fsalfirst, k) = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    (; tmp) = cache.tmp_cache
    integrator.f.f1(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.f.f2(tmp, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    return integrator.fsalfirst .+= tmp
end

function initialize!(integrator, cache::MRIGARKConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f.f1(integrator.uprev, integrator.p, integrator.t) +
        integrator.f.f2(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::MRIGARKCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; v, f1eval, slow_f, Y, fS, tab) = cache
    (; tmp, atmp) = cache.tmp_cache
    (; Γ, γ, c) = tab
    alg = unwrap_alg(integrator, false)
    m = alg.m
    s = length(Γ)

    @.. broadcast = false Y[1] = uprev
    h_inner = dt / m

    for i in 1:s
        f.f2(fS[i], Y[i], p, t + c[i] * dt)
        integrator.stats.nf2 += 1

        @.. broadcast = false slow_f = γ[i, 1] * fS[1]
        for j in 2:i
            if !iszero(γ[i, j])
                @.. broadcast = false slow_f = slow_f + γ[i, j] * fS[j]
            end
        end

        if iszero(Γ[i])
            if i < s
                @.. broadcast = false Y[i + 1] = Y[i] + dt * slow_f
            else
                @.. broadcast = false u = Y[i] + dt * slow_f
            end
        else
            @.. broadcast = false v = Y[i]
            for k_step in 1:m
                t_a = t + (k_step - 1) * h_inner
                f.f1(f1eval, v, p, t_a)
                @.. broadcast = false tmp = v + (h_inner / 2) * (Γ[i] * f1eval + slow_f)
                f.f1(f1eval, tmp, p, t_a + h_inner / 2)
                @.. broadcast = false v = v + h_inner * (Γ[i] * f1eval + slow_f)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
            end
            if i < s
                @.. broadcast = false Y[i + 1] = v
            else
                @.. broadcast = false u = v
            end
        end
    end

    return if integrator.opts.adaptive
        @.. broadcast = false tmp = u - Y[s]
        calculate_residuals!(
            atmp, tmp, uprev, u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end

@muladd function perform_step!(
        integrator, cache::MRIGARKConstantCache, repeat_step = false
    )
    (; t, dt, uprev, f, p) = integrator
    (; Γ, γ, c) = cache.tab
    alg = unwrap_alg(integrator, false)
    m = alg.m
    s = length(Γ)

    h_inner = dt / m
    Y = Vector{typeof(uprev)}(undef, s)
    fS = Vector{typeof(f.f2(uprev, p, t))}(undef, s)
    Y[1] = uprev

    for i in 1:s
        fS[i] = f.f2(Y[i], p, t + c[i] * dt)
        integrator.stats.nf2 += 1

        slow_f = γ[i, 1] * fS[1]
        for j in 2:i
            if !iszero(γ[i, j])
                slow_f = @.. broadcast = false slow_f + γ[i, j] * fS[j]
            end
        end

        if iszero(Γ[i])
            result = @.. broadcast = false Y[i] + dt * slow_f
        else
            v = Y[i]
            for k_step in 1:m
                t_a = t + (k_step - 1) * h_inner
                k1 = f.f1(v, p, t_a)
                v_mid = @.. broadcast = false v + (h_inner / 2) * (Γ[i] * k1 + slow_f)
                k2 = f.f1(v_mid, p, t_a + h_inner / 2)
                v = @.. broadcast = false v + h_inner * (Γ[i] * k2 + slow_f)
                OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
            end
            result = v
        end

        if i < s
            Y[i + 1] = result
        else
            integrator.u = result
        end
    end

    if integrator.opts.adaptive
        utilde = @.. broadcast = false integrator.u - Y[s]
        atmp = calculate_residuals(
            utilde, uprev, integrator.u,
            integrator.opts.abstol, integrator.opts.reltol,
            integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
end
