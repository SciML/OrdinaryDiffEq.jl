################################################################################
# Implicit Weak Order 2 SRK Methods
# IRI1: Implicit Rößler 1 - drift-implicit version of RI1 with weak order 2.0

"""
IRI1 perform_step! implementation (constant cache, out-of-place)

This implements a drift-implicit weak order 2 stochastic Runge-Kutta method
based on the RI1 scheme with theta-method implicitization of the drift.
"""
@muladd function perform_step!(integrator, cache::IRI1ConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    (;
        nlsolver, a021, a031, a032, a121, a131, b021, b031, b121, b131,
        b221, b222, b223, b231, b232, b233, α1, α2, α3, c02, c03, c12, c13,
        beta11, beta12, beta13, beta22, beta23, beta31, beta32, beta33,
        beta42, beta43, NORMAL_ONESIX_QUANTILE,
    ) = cache

    alg = unwrap_alg(integrator, true)
    theta = alg.theta
    OrdinaryDiffEqNonlinearSolve.markfirststage!(nlsolver)
    repeat_step = false

    # Define three-point distributed random variables
    dW_scaled = W.dW / sqrt(dt)
    sq3dt = sqrt(3 * dt)
    _dW = map(x -> calc_threepoint_random(sq3dt, NORMAL_ONESIX_QUANTILE, x), dW_scaled)
    chi1 = map(x -> (x^2 - dt) / 2, _dW)  # diagonal of Ihat2

    if !(W.dW isa Number)
        m = length(W.dW)
        _dZ = map(x -> calc_twopoint_random(integrator.sqdt, x), W.dZ)
    else
        m = 1
    end

    # Stage 1: Compute k1 implicitly
    # For implicit treatment: k1 = f(Y1) where Y1 = uprev + theta*dt*k1
    # Using nlsolver to find k1

    # Initial guess for z (where z = theta*dt*k1)
    z = zero(uprev)
    nlsolver.z = z
    nlsolver.c = theta
    nlsolver.tmp = uprev
    z = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return nothing

    # k1 = z / (theta * dt) but we store the full drift contribution
    Y1 = uprev + z  # Y1 is the implicit stage value
    k1 = integrator.f(Y1, p, t)  # Evaluate drift at implicit stage

    # Compute g1 at uprev (diffusion remains explicit)
    g1 = integrator.f.g(uprev, p, t)

    # H^(0) Stage 2
    if !is_diagonal_noise(integrator.sol.prob) || W.dW isa Number
        H02 = uprev + a021 * k1 * dt + b021 * g1 * _dW
    else
        H02 = uprev + a021 * k1 * dt + b021 * g1 .* _dW
    end

    # Stage 2: Compute k2 implicitly at H02
    nlsolver.z = zero(uprev)
    nlsolver.tmp = H02
    nlsolver.c = theta
    z2 = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return nothing

    Y2 = H02 + z2
    k2 = integrator.f(Y2, p, t + c02 * dt)

    # H^(0) Stage 3
    H03 = uprev + a032 * k2 * dt + a031 * k1 * dt
    if !is_diagonal_noise(integrator.sol.prob) || W.dW isa Number
        H03 += b031 * g1 * _dW
    else
        H03 += b031 * g1 .* _dW
    end

    # Stage 3: Compute k3 implicitly at H03
    nlsolver.z = zero(uprev)
    nlsolver.tmp = H03
    nlsolver.c = theta
    z3 = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return nothing

    Y3 = H03 + z3
    k3 = integrator.f(Y3, p, t + c03 * dt)

    # H^(k) stages for diffusion (explicit as in original RI1)
    if W.dW isa Number
        H12 = uprev + a121 * k1 * dt + b121 * g1 * integrator.sqdt
        H13 = uprev + a131 * k1 * dt + b131 * g1 * integrator.sqdt
    else
        H12 = [uprev .+ a121 * k1 * dt for k in 1:m]
        H13 = [uprev .+ a131 * k1 * dt for k in 1:m]
        for k in 1:m
            if is_diagonal_noise(integrator.sol.prob)
                tmp = zero(integrator.u)
                tmp[k] = g1[k]
                H12[k] += b121 * tmp * integrator.sqdt
                H13[k] += b131 * tmp * integrator.sqdt
            else
                H12[k] += b121 * g1[:, k] * integrator.sqdt
                H13[k] += b131 * g1[:, k] * integrator.sqdt
            end
        end
    end

    # Compute g2 and g3 at H12 and H13
    if W.dW isa Number
        g2 = integrator.f.g(H12, p, t + c12 * dt)
        g3 = integrator.f.g(H13, p, t + c13 * dt)
    else
        g2 = [integrator.f.g(H12[k], p, t + c12 * dt) for k in 1:m]
        g3 = [integrator.f.g(H13[k], p, t + c13 * dt) for k in 1:m]
        H22 = [copy(uprev) for k in 1:m]
        H23 = [copy(uprev) for k in 1:m]
        for k in 1:m
            if is_diagonal_noise(integrator.sol.prob)
                tmp = zero(integrator.u)
                tmp[k] = b221 * g1[k] + b222 * g2[k][k] + b223 * g3[k][k]
                H22[k] += tmp * integrator.sqdt
                tmp[k] = b231 * g1[k] + b232 * g2[k][k] + b233 * g3[k][k]
                H23[k] += tmp * integrator.sqdt
            else
                H22[k] += (b221 * g1[:, k] + b222 * g2[k][:, k] + b223 * g3[k][:, k]) *
                    integrator.sqdt
                H23[k] += (b231 * g1[:, k] + b232 * g2[k][:, k] + b233 * g3[k][:, k]) *
                    integrator.sqdt
            end
        end
    end

    # Final update: combine drift (implicit) and diffusion (explicit)
    u = uprev + α1 * k1 * dt + α2 * k2 * dt + α3 * k3 * dt

    # Add noise contribution (same as explicit RI1)
    if W.dW isa Number
        u += g1 * (_dW * beta11) + g2 * (_dW * beta12 + chi1 * beta22 / integrator.sqdt) +
            g3 * (_dW * beta13 + chi1 * beta23 / integrator.sqdt)
    else
        if is_diagonal_noise(integrator.sol.prob)
            u += g1 .* _dW * beta11
            for k in 1:m
                u[k] += g2[k][k] * _dW[k] * beta12 + g3[k][k] * _dW[k] * beta13
                u[k] += g2[k][k] * chi1[k] * beta22 / integrator.sqdt +
                    g3[k][k] * chi1[k] * beta23 / integrator.sqdt
                u[k] += (m - 1) * beta31 * g1[k] * _dW[k]
                for l in 1:m
                    if l != k
                        ihat2 = Ihat2(cache, _dW, _dZ, integrator.sqdt, k, l)
                        tmpg = integrator.f.g(H22[l], p, t)
                        u[k] += tmpg[k] *
                            (_dW[k] * beta32 + ihat2 * beta42 / integrator.sqdt)
                        tmpg = integrator.f.g(H23[l], p, t)
                        u[k] += tmpg[k] *
                            (_dW[k] * beta33 + ihat2 * beta43 / integrator.sqdt)
                    end
                end
            end
        else
            for k in 1:m
                g1k = @view g1[:, k]
                g2k = @view g2[k][:, k]
                g3k = @view g3[k][:, k]
                @.. u = u + g1k * _dW[k] * (beta11 + (m - 1) * beta31)
                @.. u = u + g2k * _dW[k] * beta12 + g3k * _dW[k] * beta13
                @.. u = u + g2k * chi1[k] * beta22 / integrator.sqdt +
                    g3k * chi1[k] * beta23 / integrator.sqdt
                for l in 1:m
                    if l != k
                        ihat2 = Ihat2(cache, _dW, _dZ, integrator.sqdt, k, l)
                        tmpg = integrator.f.g(H22[l], p, t)
                        tmpgk = @view tmpg[:, k]
                        @.. u = u +
                            tmpgk * (_dW[k] * beta32 + ihat2 * beta42 / integrator.sqdt)
                        tmpg = integrator.f.g(H23[l], p, t)
                        tmpgk = @view tmpg[:, k]
                        @.. u = u +
                            tmpgk * (_dW[k] * beta33 + ihat2 * beta43 / integrator.sqdt)
                    end
                end
            end
        end
    end

    # Adaptive error estimation
    if integrator.opts.adaptive &&
            (W.dW isa Number || is_diagonal_noise(integrator.sol.prob) || m == 1)
        # Compare against lower order method (similar to DRI1)
        if c03 != 0.0
            rat = c02 / c03
            δ₁ = integrator.opts.delta * (rat - 1)
            δ₂ = -integrator.opts.delta * rat
            uhat = uprev +
                ((α1 + δ₁) * k1 + (α2 + integrator.opts.delta) * k2 + (α3 + δ₂) * k3) *
                dt
        else
            uhat = uprev + k1 * dt
        end
        if is_diagonal_noise(integrator.sol.prob)
            uhat += g1 .* _dW
        else
            uhat += g1 * _dW
        end

        resids = calculate_residuals(
            u - uhat, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(resids, t))
    end

    integrator.u = u
end

# Ihat2 for IRI1ConstantCache (same as DRI1)
function Ihat2(cache::IRI1ConstantCache, _dW, _dZ, sqdt, k, l)
    # Two-point approximation of iterated integral (same as DRI1)
    if k < l
        return (_dW[k] * _dW[l] - sqdt * _dZ[k]) / 2
    elseif l < k
        return (_dW[k] * _dW[l] + sqdt * _dZ[l]) / 2
    else
        # l == k
        return (_dW[k]^2 - sqdt^2) / 2
    end
end

"""
IRI1 perform_step! implementation (mutable cache, in-place)
"""
@muladd function perform_step!(integrator, cache::IRI1Cache)
    (; t, dt, uprev, u, W, p, f) = integrator
    (;
        _dW, _dZ, chi1, g1, g2, g3, k1, k2, k3, H02, H03, H12, H13, H22, H23,
        tmp, tmpg, resids, nlsolver, uhat,
        a021, a031, a032, a121, a131, b021, b031, b121, b131,
        b221, b222, b223, b231, b232, b233, α1, α2, α3, c02, c03, c12, c13,
        beta11, beta12, beta13, beta22, beta23, beta31, beta32, beta33,
        beta42, beta43, NORMAL_ONESIX_QUANTILE,
    ) = cache
    (; z) = nlsolver

    alg = unwrap_alg(integrator, true)
    theta = alg.theta
    OrdinaryDiffEqNonlinearSolve.markfirststage!(nlsolver)
    repeat_step = false

    m = is_diagonal_noise(integrator.sol.prob) ? length(W.dW) : size(g1, 2)
    sq3dt = sqrt(3 * dt)

    # Define three-point distributed random variables
    if W.dW isa Union{SArray, Number}
        _dW = map(
            x -> calc_threepoint_random(sq3dt, NORMAL_ONESIX_QUANTILE, x),
            W.dW / sqrt(dt)
        )
    else
        sqrtdt = sqrt(dt)
        @.. chi1 = W.dW / sqrtdt
        calc_threepoint_random!(_dW, sq3dt, NORMAL_ONESIX_QUANTILE, chi1)
        map!(x -> (x^2 - dt) / 2, chi1, _dW)

        if m > 1 && _dZ !== nothing
            calc_twopoint_random!(_dZ, integrator.sqdt, W.dZ)
        end
    end

    # Compute g1 (diffusion is explicit)
    integrator.f.g(g1, uprev, p, t)

    # Stage 1: Implicit solve for k1
    @.. z = zero(eltype(u))
    nlsolver.z = z
    nlsolver.c = theta
    @.. nlsolver.tmp = uprev

    z = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return

    @.. tmp = uprev + z  # Y1
    integrator.f(k1, tmp, p, t)

    # H^(0) Stage 2
    if is_diagonal_noise(integrator.sol.prob)
        @.. H02 = uprev + a021 * k1 * dt + b021 * g1 * _dW
    else
        mul!(tmp, g1, _dW)
        @.. H02 = uprev + a021 * k1 * dt + b021 * tmp
    end

    # Stage 2: Implicit solve for k2
    @.. z = zero(eltype(u))
    nlsolver.z = z
    nlsolver.c = theta
    @.. nlsolver.tmp = H02

    z = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return

    @.. tmp = H02 + z  # Y2
    integrator.f(k2, tmp, p, t + c02 * dt)

    # H^(0) Stage 3
    @.. H03 = uprev + a032 * k2 * dt + a031 * k1 * dt
    if is_diagonal_noise(integrator.sol.prob)
        @.. H03 += b031 * g1 * _dW
    else
        mul!(tmp, g1, _dW)
        @.. H03 += b031 * tmp
    end

    # Stage 3: Implicit solve for k3
    @.. z = zero(eltype(u))
    nlsolver.z = z
    nlsolver.c = theta
    @.. nlsolver.tmp = H03

    z = OrdinaryDiffEqNonlinearSolve.nlsolve!(nlsolver, integrator, cache, repeat_step)
    OrdinaryDiffEqNonlinearSolve.nlsolvefail(nlsolver) && return

    @.. tmp = H03 + z  # Y3
    integrator.f(k3, tmp, p, t + c03 * dt)

    # H^(k) stages for diffusion (explicit)
    for k in 1:m
        @.. H12[k] = uprev + a121 * k1 * dt
        @.. H13[k] = uprev + a131 * k1 * dt
        if is_diagonal_noise(integrator.sol.prob)
            H12[k][k] += b121 * g1[k] * integrator.sqdt
            H13[k][k] += b131 * g1[k] * integrator.sqdt
        else
            g1k = @view g1[:, k]
            @.. H12[k] += b121 * g1k * integrator.sqdt
            @.. H13[k] += b131 * g1k * integrator.sqdt
        end
    end

    # Compute g2 and g3
    for k in 1:m
        integrator.f.g(g2[k], H12[k], p, t + c12 * dt)
        integrator.f.g(g3[k], H13[k], p, t + c13 * dt)
    end

    # H22 and H23 for non-diagonal noise
    if m > 1
        for k in 1:m
            @.. H22[k] = uprev
            @.. H23[k] = uprev
            if is_diagonal_noise(integrator.sol.prob)
                H22[k][k] += (b221 * g1[k] + b222 * g2[k][k] + b223 * g3[k][k]) *
                    integrator.sqdt
                H23[k][k] += (b231 * g1[k] + b232 * g2[k][k] + b233 * g3[k][k]) *
                    integrator.sqdt
            else
                g1k = @view g1[:, k]
                g2kk = @view g2[k][:, k]
                g3kk = @view g3[k][:, k]
                @.. H22[k] += (b221 * g1k + b222 * g2kk + b223 * g3kk) * integrator.sqdt
                @.. H23[k] += (b231 * g1k + b232 * g2kk + b233 * g3kk) * integrator.sqdt
            end
        end
    end

    # Final update
    @.. u = uprev + α1 * k1 * dt + α2 * k2 * dt + α3 * k3 * dt

    # Add noise contribution
    if is_diagonal_noise(integrator.sol.prob)
        @.. u += g1 * _dW * beta11
        for k in 1:m
            u[k] += g2[k][k] * _dW[k] * beta12 + g3[k][k] * _dW[k] * beta13
            u[k] += g2[k][k] * chi1[k] * beta22 / integrator.sqdt +
                g3[k][k] * chi1[k] * beta23 / integrator.sqdt
            u[k] += (m - 1) * beta31 * g1[k] * _dW[k]
            for l in 1:m
                if l != k
                    ihat2 = Ihat2(cache, _dW, _dZ, integrator.sqdt, k, l)
                    integrator.f.g(tmpg, H22[l], p, t)
                    u[k] += tmpg[k] * (_dW[k] * beta32 + ihat2 * beta42 / integrator.sqdt)
                    integrator.f.g(tmpg, H23[l], p, t)
                    u[k] += tmpg[k] * (_dW[k] * beta33 + ihat2 * beta43 / integrator.sqdt)
                end
            end
        end
    else
        for k in 1:m
            g1k = @view g1[:, k]
            g2k = @view g2[k][:, k]
            g3k = @view g3[k][:, k]
            @.. u = u + g1k * _dW[k] * (beta11 + (m - 1) * beta31)
            @.. u = u + g2k * _dW[k] * beta12 + g3k * _dW[k] * beta13
            @.. u = u + g2k * chi1[k] * beta22 / integrator.sqdt +
                g3k * chi1[k] * beta23 / integrator.sqdt
            for l in 1:m
                if l != k
                    ihat2 = Ihat2(cache, _dW, _dZ, integrator.sqdt, k, l)
                    integrator.f.g(tmpg, H22[l], p, t)
                    tmpgk = @view tmpg[:, k]
                    @.. u = u + tmpgk * (_dW[k] * beta32 + ihat2 * beta42 / integrator.sqdt)
                    integrator.f.g(tmpg, H23[l], p, t)
                    tmpgk = @view tmpg[:, k]
                    @.. u = u + tmpgk * (_dW[k] * beta33 + ihat2 * beta43 / integrator.sqdt)
                end
            end
        end
    end

    # Adaptive error estimation
    if integrator.opts.adaptive &&
            (W.dW isa Number || is_diagonal_noise(integrator.sol.prob) || m == 1)
        if c03 != 0.0
            rat = c02 / c03
            δ₁ = integrator.opts.delta * (rat - 1)
            δ₂ = -integrator.opts.delta * rat
            @.. uhat = uprev +
                (
                (α1 + δ₁) * k1 + (α2 + integrator.opts.delta) * k2 +
                    (α3 + δ₂) * k3
            ) * dt
        else
            @.. uhat = uprev + k1 * dt
        end
        if is_diagonal_noise(integrator.sol.prob)
            @.. uhat += g1 * _dW
        else
            mul!(tmp, g1, _dW)
            @.. uhat += tmp
        end

        @.. tmp = u - uhat
        calculate_residuals!(
            resids, tmp, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(resids, t))
    end
end

# Ihat2 for IRI1Cache (same as DRI1)
function Ihat2(cache::IRI1Cache, _dW, _dZ, sqdt, k, l)
    # Two-point approximation of iterated integral (same as DRI1)
    if k < l
        return (_dW[k] * _dW[l] - sqdt * _dZ[k]) / 2
    elseif l < k
        return (_dW[k] * _dW[l] + sqdt * _dZ[l]) / 2
    else
        # l == k
        return (_dW[k]^2 - sqdt^2) / 2
    end
end
