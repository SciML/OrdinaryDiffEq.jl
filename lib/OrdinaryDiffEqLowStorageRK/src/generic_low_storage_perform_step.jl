@muladd function _perform_step_iip!(
        integrator, cache, tab::LowStorageRKTableau{TwoN}
    )
    (; t, dt, u, f, p) = integrator
    (; k, tmp, williamson_condition, stage_limiter!, step_limiter!, thread) = cache
    (; A2end, B1, B2end, c2end) = tab

    f(k, u, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false thread = thread tmp = dt * k
    @.. broadcast = false thread = thread u = u + B1 * tmp

    for i in eachindex(A2end)
        if williamson_condition
            f(ArrayFuse(tmp, u, (A2end[i], dt, B2end[i])), u, p, t + c2end[i] * dt)
        else
            @.. broadcast = false thread = thread tmp = A2end[i] * tmp
            stage_limiter!(u, integrator, p, t + c2end[i] * dt)
            f(k, u, p, t + c2end[i] * dt)
            @.. broadcast = false thread = thread tmp = tmp + dt * k
            @.. broadcast = false thread = thread u = u + B2end[i] * tmp
        end
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end
    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)
    return nothing
end

@muladd function _perform_step_oop!(
        integrator, tab::LowStorageRKTableau{TwoN}
    )
    (; t, dt, u, f, p) = integrator
    (; A2end, B1, B2end, c2end) = tab

    tmp = dt * integrator.fsalfirst
    u = u + B1 * tmp

    for i in eachindex(A2end)
        k = f(u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tmp = A2end[i] * tmp + dt * k
        u = u + B2end[i] * tmp
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.fsalfirst = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.u = u
    return nothing
end

@muladd function _perform_step_iip!(
        integrator, cache, tab::LowStorageRKTableau{TwoC}
    )
    (; t, dt, u, f, p) = integrator
    (; k, fsalfirst, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; A2end, B1, B2end, c2end) = tab

    @.. broadcast = false thread = thread k = integrator.fsalfirst
    @.. broadcast = false thread = thread u = u + B1 * dt * k

    for i in eachindex(A2end)
        @.. broadcast = false thread = thread tmp = u + A2end[i] * dt * k
        f(k, tmp, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false thread = thread u = u + B2end[i] * dt * k
    end
    step_limiter!(u, integrator, p, t + dt)
    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return nothing
end

@muladd function _perform_step_oop!(
        integrator, tab::LowStorageRKTableau{TwoC}
    )
    (; t, dt, u, f, p) = integrator
    (; A2end, B1, B2end, c2end) = tab

    k = integrator.fsalfirst = integrator.f(u, p, t)
    integrator.k[1] = integrator.fsalfirst
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    u = u + B1 * dt * k

    for i in eachindex(A2end)
        tmp = u + A2end[i] * dt * k
        k = integrator.f(tmp, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        u = u + B2end[i] * dt * k
    end

    integrator.u = u
    return nothing
end

@muladd function _perform_step_oop!(integrator, tab::LowStorageRK3SConstantCache)
    (; t, dt, uprev, u, f, p) = integrator
    (; γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end) = tab

    tmp = u
    u = tmp + β1 * dt * integrator.fsalfirst

    for i in eachindex(γ12end)
        k = f(u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tmp = tmp + δ2end[i] * u
        u = γ12end[i] * u + γ22end[i] * tmp + γ32end[i] * uprev + β2end[i] * dt * k
    end

    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.u = u
    return nothing
end

@muladd function _perform_step_iip!(integrator, cache, tab::LowStorageRK3SConstantCache)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, fsalfirst, tmp, stage_limiter!, step_limiter!, thread) = cache
    (; γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end) = tab

    @.. broadcast = false thread = thread tmp = u
    @.. broadcast = false thread = thread u = tmp + β1 * dt * integrator.fsalfirst

    for i in eachindex(γ12end)
        f(k, u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false thread = thread tmp = tmp + δ2end[i] * u
        @.. broadcast = false thread = thread u = γ12end[i] * u + γ22end[i] * tmp +
            γ32end[i] * uprev +
            β2end[i] * dt * k
    end

    step_limiter!(u, integrator, p, t + dt)
    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return nothing
end

@muladd function _perform_step_oop!(integrator, tab::LowStorageRK3SpConstantCache)
    (; t, dt, uprev, u, f, p) = integrator
    (; γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end) = tab

    integrator.fsalfirst = f(uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    tmp = uprev
    u = tmp + β1 * dt * integrator.fsalfirst
    utilde = u
    if integrator.opts.adaptive
        utilde = bhat1 * dt * integrator.fsalfirst
    end

    for i in eachindex(γ12end)
        k = f(u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tmp = tmp + δ2end[i] * u
        u = γ12end[i] * u + γ22end[i] * tmp + γ32end[i] * uprev + β2end[i] * dt * k
        if integrator.opts.adaptive
            utilde = utilde + bhat2end[i] * dt * k
        end
    end

    if integrator.opts.adaptive
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end

    integrator.u = u
    return nothing
end

@muladd function _perform_step_iip!(integrator, cache, tab::LowStorageRK3SpConstantCache)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, tmp, utilde, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end) = tab

    f(integrator.fsalfirst, uprev, p, t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false thread = thread tmp = uprev
    @.. broadcast = false thread = thread u = tmp + β1 * dt * integrator.fsalfirst
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = bhat1 * dt * integrator.fsalfirst
    end

    for i in eachindex(γ12end)
        stage_limiter!(u, integrator, p, t + c2end[i] * dt)
        f(k, u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false thread = thread tmp = tmp + δ2end[i] * u
        @.. broadcast = false thread = thread u = γ12end[i] * u + γ22end[i] * tmp +
            γ32end[i] * uprev + β2end[i] * dt * k
        if integrator.opts.adaptive
            @.. broadcast = false thread = thread utilde = utilde + bhat2end[i] * dt * k
        end
    end

    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    return nothing
end

@muladd function _perform_step_oop!(integrator, tab::LowStorageRK3SpFSALConstantCache)
    (; t, dt, uprev, u, f, p) = integrator
    (; γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end, bhatfsal) = tab

    tmp = uprev
    u = tmp + β1 * dt * integrator.fsalfirst
    utilde = u
    if integrator.opts.adaptive
        utilde = bhat1 * dt * integrator.fsalfirst
    end

    for i in eachindex(γ12end)
        k = f(u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tmp = tmp + δ2end[i] * u
        u = γ12end[i] * u + γ22end[i] * tmp + γ32end[i] * uprev + β2end[i] * dt * k
        if integrator.opts.adaptive
            utilde = utilde + bhat2end[i] * dt * k
        end
    end

    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        utilde = utilde + bhatfsal * dt * integrator.fsallast
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end

    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return nothing
end

@muladd function _perform_step_iip!(integrator, cache, tab::LowStorageRK3SpFSALConstantCache)
    (; t, dt, uprev, u, f, p) = integrator
    (; k, tmp, utilde, atmp, stage_limiter!, step_limiter!, thread) = cache
    (; γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end, bhatfsal) = tab

    @.. broadcast = false thread = thread tmp = uprev
    @.. broadcast = false thread = thread u = tmp + β1 * dt * integrator.fsalfirst
    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = bhat1 * dt * integrator.fsalfirst
    end

    for i in eachindex(γ12end)
        stage_limiter!(u, integrator, p, t + c2end[i] * dt)
        f(k, u, p, t + c2end[i] * dt)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false thread = thread tmp = tmp + δ2end[i] * u
        @.. broadcast = false thread = thread u = γ12end[i] * u + γ22end[i] * tmp +
            γ32end[i] * uprev + β2end[i] * dt * k
        if integrator.opts.adaptive
            @.. broadcast = false thread = thread utilde = utilde + bhat2end[i] * dt * k
        end
    end

    stage_limiter!(u, integrator, p, t + dt)
    step_limiter!(u, integrator, p, t + dt)

    f(k, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    if integrator.opts.adaptive
        @.. broadcast = false thread = thread utilde = utilde + bhatfsal * dt * k
        calculate_residuals!(
            atmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t,
            thread
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end
    return nothing
end
