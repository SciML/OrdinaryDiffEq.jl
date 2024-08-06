function do_newJ(integrator, alg, cache, repeat_step)::Bool # for FIRK
    integrator.iter <= 1 && return true
    repeat_step && return false
    first(islinearfunction(integrator)) && return false
    integrator.opts.adaptive || return true
    alg_can_repeat_jac(alg) || return true
    integrator.u_modified && return true
    # below is Newton specific logic, so we return non-Newton algs here
    alg isa NewtonAlgorithm || return true
    nlstatus = cache.status
    Int8(nlstatus) < 0 && return true
    # no reuse when the cutoff is 0
    fast_convergence_cutoff = alg.fast_convergence_cutoff
    iszero(fast_convergence_cutoff) && return true
    # reuse J when there is fast convergence
    fastconvergence = nlstatus === FastConvergence
    return !fastconvergence
end

function do_newW(integrator, nlsolver, new_jac, W_dt)::Bool # for FIRK
    nlsolver === nothing && return true
    new_jac && return true
    # reuse W when the change in stepsize is small enough
    dt = integrator.dt
    smallstepchange = abs((dt - W_dt) / W_dt) <= get_new_W_γdt_cutoff(nlsolver)
    return !smallstepchange
end

function initialize!(integrator, cache::RadauIIA3ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    nothing
end

function initialize!(integrator, cache::RadauIIA5ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.stats.nf += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    nothing
end

function initialize!(integrator, cache::RadauIIA9ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.stats.nf += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    nothing
end

function initialize!(integrator, cache::RadauIIA3Cache)
    integrator.kshortsize = 2
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    nothing
end

function initialize!(integrator, cache::RadauIIA5Cache)
    integrator.kshortsize = 2
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    if integrator.opts.adaptive
        @unpack abstol, reltol = integrator.opts
        if reltol isa Number
            cache.rtol = reltol^(2 / 3) / 10
            cache.atol = cache.rtol * (abstol / reltol)
        else
            @.. broadcast=false cache.rtol=reltol^(2 / 3) / 10
            @.. broadcast=false cache.atol=cache.rtol * (abstol / reltol)
        end
    end
    nothing
end

function initialize!(integrator, cache::RadauIIA9Cache)
    integrator.kshortsize = 2
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    if integrator.opts.adaptive
        @unpack abstol, reltol = integrator.opts
        if reltol isa Number
            cache.rtol = reltol^(5 / 8) / 10
            cache.atol = cache.rtol * (abstol / reltol)
        else
            @.. broadcast=false cache.rtol=reltol^(5 / 8) / 10
            @.. broadcast=false cache.atol=cache.rtol * (abstol / reltol)
        end
    end
    nothing
end

@muladd function perform_step!(integrator, cache::RadauIIA3ConstantCache)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack T11, T12, T21, T22, TI11, TI12, TI21, TI22 = cache.tab
    @unpack c1, c2, α, β, e1, e2 = cache.tab
    @unpack κ, cont1, cont2 = cache
    @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
    alg = unwrap_alg(integrator, true)
    @unpack maxiters = alg
    mass_matrix = integrator.f.mass_matrix

    # precalculations
    rtol = @. reltol^(2 / 3) / 10
    atol = @. rtol * (abstol / reltol)
    αdt, βdt = α / dt, β / dt
    J = calc_J(integrator, cache)

    cache.dtprev = one(cache.dtprev)
    z1 = w1 = map(zero, u)
    z2 = w2 = map(zero, u)
    cache.cont1 = map(zero, u)
    cache.cont2 = map(zero, u)

    if u isa Number
        LU1 = -(αdt + βdt * im) * mass_matrix + J
    else
        LU1 = lu(-(αdt + βdt * im) * mass_matrix + J)
    end
    integrator.stats.nw += 1

    # Newton iteration
    local ndw, ff1, ff2
    η = max(cache.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)
    fail_convergence = true
    iter = 0
    while iter < maxiters
        iter += 1
        integrator.stats.nnonliniter += 1
        # evaluate function
        ff1 = f(uprev + z1, p, t + c1 * dt)
        ff2 = f(uprev + z2, p, t + c2 * dt)
        integrator.stats.nf += 2

        fw1 = @. TI11 * ff1 + TI12 * ff2
        fw2 = @. TI21 * ff1 + TI22 * ff2

        if mass_matrix isa UniformScaling
            Mw1 = @. mass_matrix.λ * w1
            Mw2 = @. mass_matrix.λ * w2
        else
            Mw1 = mass_matrix * w1
            Mw2 = mass_matrix * w2
        end

        rhs1 = @. fw1 - αdt * Mw1 + βdt * Mw2
        rhs2 = @. fw2 - βdt * Mw1 - αdt * Mw2
        dw12 = LU1 \ (@. rhs1 + rhs2 * im)
        integrator.stats.nsolve += 1
        dw1 = real(dw12)
        dw2 = imag(dw12)

        # compute norm of residuals
        iter > 1 && (ndwprev = ndw)
        atmp1 = calculate_residuals(dw1, uprev, u, atol, rtol, internalnorm, t)
        atmp2 = calculate_residuals(dw2, uprev, u, atol, rtol, internalnorm, t)
        ndw = internalnorm(atmp1, t) + internalnorm(atmp2, t)

        # check divergence (not in initial step)
        if iter > 1
            θ = ndw / ndwprev
            (diverge = θ > 1) && (cache.status = Divergence)
            (veryslowconvergence = ndw * θ^(maxiters - iter) > κ * (1 - θ)) &&
                (cache.status = VerySlowConvergence)
            if diverge || veryslowconvergence
                break
            end
        end

        w1 = @. w1 - dw1
        w2 = @. w2 - dw2

        # transform `w` to `z`
        z1 = @. T11 * w1 + T12 * w2
        z2 = @. T21 * w1 + T22 * w2

        # check stopping criterion
        iter > 1 && (η = θ / (1 - θ))
        if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
            # Newton method converges
            cache.status = η < alg.fast_convergence_cutoff ? FastConvergence :
                           Convergence
            fail_convergence = false
            break
        end
    end

    cache.ηold = η
    cache.iter = iter

    u = @. uprev + z2

    if adaptive
        utilde = @. dt * (e1 * ff1 + e2 * ff2)
        atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm, t)
        integrator.EEst = internalnorm(atmp, t)
    end

    integrator.fsallast = f(u, p, t + dt)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return
end

@muladd function perform_step!(integrator, cache::RadauIIA3Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p, fsallast, fsalfirst = integrator
    @unpack T11, T12, T21, T22, TI11, TI12, TI21, TI22 = cache.tab
    @unpack c1, c2, α, β, e1, e2 = cache.tab
    @unpack κ, cont1, cont2 = cache
    @unpack z1, z2, w1, w2,
    dw12, cubuff,
    k, k2, fw1, fw2,
    J, W1,
    tmp, atmp, jac_config, rtol, atol, step_limiter! = cache
    @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
    alg = unwrap_alg(integrator, true)
    @unpack maxiters = alg
    mass_matrix = integrator.f.mass_matrix
    # precalculations
    αdt, βdt = α / dt, β / dt
    (new_jac = do_newJ(integrator, alg, cache, repeat_step)) &&
        (calc_J!(J, integrator, cache); cache.W_γdt = dt)
    if (new_W = do_newW(integrator, alg, new_jac, cache.W_γdt))
        @inbounds for II in CartesianIndices(J)
            W1[II] = -(αdt + βdt * im) * mass_matrix[Tuple(II)...] + J[II]
        end
        integrator.stats.nw += 1
    end

    #better initial guess
    uzero = zero(eltype(z1))
    @. z1 = uzero
    @. z2 = uzero
    @. w1 = uzero
    @. w2 = uzero
    @. cache.cont1 = uzero
    @. cache.cont2 = uzero

    # Newton iteration
    local ndw
    η = max(cache.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)
    fail_convergence = true
    iter = 0
    while iter < maxiters
        iter += 1
        integrator.stats.nnonliniter += 1
        # evaluate function
        @. tmp = uprev + z1
        f(fsallast, tmp, p, t + c1 * dt)
        @. tmp = uprev + z2
        f(k2, tmp, p, t + c2 * dt)
        integrator.stats.nf += 2

        @. fw1 = TI11 * fsallast + TI12 * k2
        @. fw2 = TI21 * fsallast + TI22 * k2

        if mass_matrix === I
            Mw1 = w1
            Mw2 = w2
        elseif mass_matrix isa UniformScaling
            mul!(z1, mass_matrix.λ, w1)
            mul!(z2, mass_matrix.λ, w2)
            Mw1 = z1
            Mw2 = z2
        else
            mul!(z1, mass_matrix, w1)
            mul!(z2, mass_matrix, w2)
            Mw1 = z1
            Mw2 = z2
        end

        @. cubuff = complex(fw1 - αdt * Mw1 + βdt * Mw2, fw2 - βdt * Mw1 - αdt * Mw2)
        needfactor = iter == 1

        linsolve = cache.linsolve

        if needfactor
            linres = dolinsolve(integrator, linsolve; A = W1, b = _vec(cubuff),
                linu = _vec(dw12))
        else
            linres = dolinsolve(integrator, linsolve; A = nothing, b = _vec(cubuff),
                linu = _vec(dw12))
        end

        cache.linsolve = linres.cache

        integrator.stats.nsolve += 1
        dw1 = real(dw12)
        dw2 = imag(dw12)

        # compute norm of residuals
        iter > 1 && (ndwprev = ndw)
        calculate_residuals!(atmp, dw1, uprev, u, atol, rtol, internalnorm, t)
        ndw1 = internalnorm(atmp, t)
        calculate_residuals!(atmp, dw2, uprev, u, atol, rtol, internalnorm, t)
        ndw2 = internalnorm(atmp, t)
        ndw = ndw1 + ndw2

        # check divergence (not in initial step)
        if iter > 1
            θ = ndw / ndwprev
            (diverge = θ > 2) && (cache.status = Divergence)
            if diverge
                break
            end
        end

        @. w1 = w1 - dw1
        @. w2 = w2 - dw2

        # transform `w` to `z`
        @. z1 = T11 * w1 + T12 * w2
        @. z2 = T21 * w1 + T22 * w2

        # check stopping criterion
        iter > 1 && (η = θ / (1 - θ))
        if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
            # Newton method converges
            cache.status = η < alg.fast_convergence_cutoff ? FastConvergence :
                           Convergence
            fail_convergence = false
            break
        end
    end
    if fail_convergence
        integrator.force_stepfail = true
        integrator.stats.nnonlinconvfail += 1
        return
    end
    cache.ηold = η
    cache.iter = iter

    @. u = uprev + z2
    step_limiter!(u, integrator, p, t + dt)

    if adaptive
        utilde = w2
        @. utilde = dt * (e1 * fsallast + e2 * k2)
        calculate_residuals!(atmp, utilde, uprev, u, atol, rtol, internalnorm, t)
        integrator.EEst = internalnorm(atmp, t)
    end
    f(fsallast, u, p, t + dt)
    integrator.stats.nf += 1
    return
end

@muladd function perform_step!(integrator, cache::RadauIIA5ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack T11, T12, T13, T21, T22, T23, T31, TI11, TI12, TI13, TI21, TI22, TI23, TI31, TI32, TI33 = cache.tab
    @unpack c1, c2, γ, α, β, e1, e2, e3 = cache.tab
    @unpack κ, cont1, cont2, cont3 = cache
    @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
    alg = unwrap_alg(integrator, true)
    @unpack maxiters = alg
    mass_matrix = integrator.f.mass_matrix

    # precalculations
    rtol = @.. broadcast=false reltol^(2 / 3)/10
    atol = @.. broadcast=false rtol*(abstol / reltol)
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    c1mc2 = c1 - c2
    γdt, αdt, βdt = γ / dt, α / dt, β / dt
    J = calc_J(integrator, cache)
    if u isa Number
        LU1 = -γdt * mass_matrix + J
        LU2 = -(αdt + βdt * im) * mass_matrix + J
    else
        LU1 = lu(-γdt * mass_matrix + J)
        LU2 = lu(-(αdt + βdt * im) * mass_matrix + J)
    end
    integrator.stats.nw += 1

    # TODO better initial guess
    if integrator.iter == 1 || integrator.u_modified || alg.extrapolant == :constant
        cache.dtprev = one(cache.dtprev)
        z1 = w1 = map(zero, u)
        z2 = w2 = map(zero, u)
        z3 = w3 = map(zero, u)
        cache.cont1 = map(zero, u)
        cache.cont2 = map(zero, u)
        cache.cont3 = map(zero, u)
    else
        c3′ = dt / cache.dtprev
        c1′ = c1 * c3′
        c2′ = c2 * c3′
        z1 = @.. broadcast=false c1′*(cont1 + (c1′ - c2m1) * (cont2 + (c1′ - c1m1) * cont3))
        z2 = @.. broadcast=false c2′*(cont1 + (c2′ - c2m1) * (cont2 + (c2′ - c1m1) * cont3))
        z3 = @.. broadcast=false c3′*(cont1 + (c3′ - c2m1) * (cont2 + (c3′ - c1m1) * cont3))
        w1 = @.. broadcast=false TI11*z1+TI12*z2+TI13*z3
        w2 = @.. broadcast=false TI21*z1+TI22*z2+TI23*z3
        w3 = @.. broadcast=false TI31*z1+TI32*z2+TI33*z3
    end

    # Newton iteration
    local ndw
    η = max(cache.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)
    fail_convergence = true
    iter = 0
    while iter < maxiters
        iter += 1
        integrator.stats.nnonliniter += 1

        # evaluate function
        ff1 = f(uprev + z1, p, t + c1 * dt)
        ff2 = f(uprev + z2, p, t + c2 * dt)
        ff3 = f(uprev + z3, p, t + dt) # c3 = 1
        integrator.stats.nf += 3

        fw1 = @.. broadcast=false TI11*ff1+TI12*ff2+TI13*ff3
        fw2 = @.. broadcast=false TI21*ff1+TI22*ff2+TI23*ff3
        fw3 = @.. broadcast=false TI31*ff1+TI32*ff2+TI33*ff3

        if mass_matrix isa UniformScaling # `UniformScaling` doesn't play nicely with broadcast
            Mw1 = @.. broadcast=false mass_matrix.λ*w1
            Mw2 = @.. broadcast=false mass_matrix.λ*w2
            Mw3 = @.. broadcast=false mass_matrix.λ*w3
        else
            Mw1 = mass_matrix * w1
            Mw2 = mass_matrix * w2
            Mw3 = mass_matrix * w3
        end

        rhs1 = @.. broadcast=false fw1-γdt * Mw1
        rhs2 = @.. broadcast=false fw2 - αdt * Mw2+βdt * Mw3
        rhs3 = @.. broadcast=false fw3 - βdt * Mw2-αdt * Mw3
        dw1 = LU1 \ rhs1
        dw23 = LU2 \ (@.. broadcast=false rhs2+rhs3 * im)
        integrator.stats.nsolve += 2
        dw2 = real(dw23)
        dw3 = imag(dw23)

        # compute norm of residuals
        iter > 1 && (ndwprev = ndw)
        atmp1 = calculate_residuals(dw1, uprev, u, atol, rtol, internalnorm, t)
        atmp2 = calculate_residuals(dw2, uprev, u, atol, rtol, internalnorm, t)
        atmp3 = calculate_residuals(dw3, uprev, u, atol, rtol, internalnorm, t)
        ndw = internalnorm(atmp1, t) + internalnorm(atmp2, t) + internalnorm(atmp3, t)
        # check divergence (not in initial step)
        if iter > 1
            θ = ndw / ndwprev
            (diverge = θ > 1) && (cache.status = Divergence)
            (veryslowconvergence = ndw * θ^(maxiters - iter) > κ * (1 - θ)) &&
                (cache.status = VerySlowConvergence)
            if diverge || veryslowconvergence
                break
            end
        end

        w1 = @.. broadcast=false w1-dw1
        w2 = @.. broadcast=false w2-dw2
        w3 = @.. broadcast=false w3-dw3

        # transform `w` to `z`
        z1 = @.. broadcast=false T11*w1+T12*w2+T13*w3
        z2 = @.. broadcast=false T21*w1+T22*w2+T23*w3
        z3 = @.. broadcast=false T31 * w1+w2           # T32 = 1, T33 = 0

        # check stopping criterion
        iter > 1 && (η = θ / (1 - θ))
        if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
            # Newton method converges
            cache.status = η < alg.fast_convergence_cutoff ? FastConvergence :
                           Convergence
            fail_convergence = false
            break
        end
    end
    if fail_convergence
        integrator.force_stepfail = true
        integrator.stats.nnonlinconvfail += 1
        return
    end
    cache.ηold = η
    cache.iter = iter

    u = @.. broadcast=false uprev+z3

    if adaptive
        e1dt, e2dt, e3dt = e1 / dt, e2 / dt, e3 / dt
        tmp = @.. broadcast=false e1dt*z1+e2dt*z2+e3dt*z3
        mass_matrix != I && (tmp = mass_matrix * tmp)
        utilde = @.. broadcast=false integrator.fsalfirst+tmp
        alg.smooth_est && (utilde = LU1 \ utilde; integrator.stats.nsolve += 1)
        # RadauIIA5 needs a transformed rtol and atol see
        # https://github.com/luchr/ODEInterface.jl/blob/0bd134a5a358c4bc13e0fb6a90e27e4ee79e0115/src/radau5.f#L399-L421
        atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm, t)
        integrator.EEst = internalnorm(atmp, t)

        if !(integrator.EEst < oneunit(integrator.EEst)) && integrator.iter == 1 ||
           integrator.u_modified
            f0 = f(uprev .+ utilde, p, t)
            integrator.stats.nf += 1
            utilde = @.. broadcast=false f0+tmp
            alg.smooth_est && (utilde = LU1 \ utilde; integrator.stats.nsolve += 1)
            atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm, t)
            integrator.EEst = internalnorm(atmp, t)
        end
    end

    if integrator.EEst <= oneunit(integrator.EEst)
        cache.dtprev = dt
        if alg.extrapolant != :constant
            cache.cont1 = @.. broadcast=false (z2 - z3)/c2m1
            tmp = @.. broadcast=false (z1 - z2)/c1mc2
            cache.cont2 = @.. broadcast=false (tmp - cache.cont1)/c1m1
            cache.cont3 = @.. broadcast=false cache.cont2-(tmp - z1 / c1) / c2
        end
    end

    integrator.fsallast = f(u, p, t + dt)
    integrator.stats.nf += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return
end

@muladd function perform_step!(integrator, cache::RadauIIA5Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p, fsallast, fsalfirst = integrator
    @unpack T11, T12, T13, T21, T22, T23, T31, TI11, TI12, TI13, TI21, TI22, TI23, TI31, TI32, TI33 = cache.tab
    @unpack c1, c2, γ, α, β, e1, e2, e3 = cache.tab
    @unpack κ, cont1, cont2, cont3 = cache
    @unpack z1, z2, z3, w1, w2, w3,
    dw1, ubuff, dw23, cubuff,
    k, k2, k3, fw1, fw2, fw3,
    J, W1, W2,
    tmp, atmp, jac_config, linsolve1, linsolve2, rtol, atol, step_limiter! = cache
    @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
    alg = unwrap_alg(integrator, true)
    @unpack maxiters = alg
    mass_matrix = integrator.f.mass_matrix

    # precalculations
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    c1mc2 = c1 - c2
    γdt, αdt, βdt = γ / dt, α / dt, β / dt
    (new_jac = do_newJ(integrator, alg, cache, repeat_step)) &&
        (calc_J!(J, integrator, cache); cache.W_γdt = dt)
    if (new_W = do_newW(integrator, alg, new_jac, cache.W_γdt))
        @inbounds for II in CartesianIndices(J)
            W1[II] = -γdt * mass_matrix[Tuple(II)...] + J[II]
            W2[II] = -(αdt + βdt * im) * mass_matrix[Tuple(II)...] + J[II]
        end
        integrator.stats.nw += 1
    end

    # TODO better initial guess
    if integrator.iter == 1 || integrator.u_modified || alg.extrapolant == :constant
        cache.dtprev = one(cache.dtprev)
        uzero = zero(eltype(u))
        @.. broadcast=false z1=uzero
        @.. broadcast=false z2=uzero
        @.. broadcast=false z3=uzero
        @.. broadcast=false w1=uzero
        @.. broadcast=false w2=uzero
        @.. broadcast=false w3=uzero
        @.. broadcast=false cache.cont1=uzero
        @.. broadcast=false cache.cont2=uzero
        @.. broadcast=false cache.cont3=uzero
    else
        c3′ = dt / cache.dtprev
        c1′ = c1 * c3′
        c2′ = c2 * c3′
        @.. broadcast=false z1=c1′ * (cont1 + (c1′ - c2m1) * (cont2 + (c1′ - c1m1) * cont3))
        @.. broadcast=false z2=c2′ * (cont1 + (c2′ - c2m1) * (cont2 + (c2′ - c1m1) * cont3))
        @.. broadcast=false z3=c3′ * (cont1 + (c3′ - c2m1) * (cont2 + (c3′ - c1m1) * cont3))
        @.. broadcast=false w1=TI11 * z1 + TI12 * z2 + TI13 * z3
        @.. broadcast=false w2=TI21 * z1 + TI22 * z2 + TI23 * z3
        @.. broadcast=false w3=TI31 * z1 + TI32 * z2 + TI33 * z3
    end

    # Newton iteration
    local ndw
    η = max(cache.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)
    fail_convergence = true
    iter = 0
    while iter < maxiters
        iter += 1
        integrator.stats.nnonliniter += 1

        # evaluate function
        @.. broadcast=false tmp=uprev + z1
        f(fsallast, tmp, p, t + c1 * dt)
        @.. broadcast=false tmp=uprev + z2
        f(k2, tmp, p, t + c2 * dt)
        @.. broadcast=false tmp=uprev + z3
        f(k3, tmp, p, t + dt) # c3 = 1
        integrator.stats.nf += 3

        @.. broadcast=false fw1=TI11 * fsallast + TI12 * k2 + TI13 * k3
        @.. broadcast=false fw2=TI21 * fsallast + TI22 * k2 + TI23 * k3
        @.. broadcast=false fw3=TI31 * fsallast + TI32 * k2 + TI33 * k3

        if mass_matrix === I
            Mw1 = w1
            Mw2 = w2
            Mw3 = w3
        elseif mass_matrix isa UniformScaling
            mul!(z1, mass_matrix.λ, w1)
            mul!(z2, mass_matrix.λ, w2)
            mul!(z3, mass_matrix.λ, w3)
            Mw1 = z1
            Mw2 = z2
            Mw3 = z3
        else
            mul!(z1, mass_matrix, w1)
            mul!(z2, mass_matrix, w2)
            mul!(z3, mass_matrix, w3)
            Mw1 = z1
            Mw2 = z2
            Mw3 = z3
        end

        @.. broadcast=false ubuff=fw1 - γdt * Mw1
        needfactor = iter == 1 && new_W

        linsolve1 = cache.linsolve1

        if needfactor
            linres1 = dolinsolve(integrator, linsolve1; A = W1, b = _vec(ubuff),
                linu = _vec(dw1))
        else
            linres1 = dolinsolve(integrator, linsolve1; A = nothing, b = _vec(ubuff),
                linu = _vec(dw1))
        end

        cache.linsolve1 = linres1.cache

        @.. broadcast=false cubuff=complex(fw2 - αdt * Mw2 + βdt * Mw3,
            fw3 - βdt * Mw2 - αdt * Mw3)

        linsolve2 = cache.linsolve2

        if needfactor
            linres2 = dolinsolve(integrator, linsolve2; A = W2, b = _vec(cubuff),
                linu = _vec(dw23))
        else
            linres2 = dolinsolve(integrator, linsolve2; A = nothing, b = _vec(cubuff),
                linu = _vec(dw23))
        end

        cache.linsolve2 = linres2.cache

        integrator.stats.nsolve += 2
        dw2 = z2
        dw3 = z3
        @.. broadcast=false dw2=real(dw23)
        @.. broadcast=false dw3=imag(dw23)

        # compute norm of residuals
        iter > 1 && (ndwprev = ndw)
        calculate_residuals!(atmp, dw1, uprev, u, atol, rtol, internalnorm, t)
        ndw1 = internalnorm(atmp, t)
        calculate_residuals!(atmp, dw2, uprev, u, atol, rtol, internalnorm, t)
        ndw2 = internalnorm(atmp, t)
        calculate_residuals!(atmp, dw3, uprev, u, atol, rtol, internalnorm, t)
        ndw3 = internalnorm(atmp, t)
        ndw = ndw1 + ndw2 + ndw3

        # check divergence (not in initial step)
        if iter > 1
            θ = ndw / ndwprev
            (diverge = θ > 1) && (cache.status = Divergence)
            (veryslowconvergence = ndw * θ^(maxiters - iter) > κ * (1 - θ)) &&
                (cache.status = VerySlowConvergence)
            if diverge || veryslowconvergence
                break
            end
        end

        @.. broadcast=false w1=w1 - dw1
        @.. broadcast=false w2=w2 - dw2
        @.. broadcast=false w3=w3 - dw3

        # transform `w` to `z`
        @.. broadcast=false z1=T11 * w1 + T12 * w2 + T13 * w3
        @.. broadcast=false z2=T21 * w1 + T22 * w2 + T23 * w3
        @.. broadcast=false z3=T31 * w1 + w2           # T32 = 1, T33 = 0

        # check stopping criterion
        iter > 1 && (η = θ / (1 - θ))
        if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
            # Newton method converges
            cache.status = η < alg.fast_convergence_cutoff ? FastConvergence :
                           Convergence
            fail_convergence = false
            break
        end
    end
    if fail_convergence
        integrator.force_stepfail = true
        integrator.stats.nnonlinconvfail += 1
        return
    end
    cache.ηold = η
    cache.iter = iter

    @.. broadcast=false u=uprev + z3
    step_limiter!(u, integrator, p, t + dt)

    if adaptive
        utilde = w2
        e1dt, e2dt, e3dt = e1 / dt, e2 / dt, e3 / dt
        @.. broadcast=false tmp=e1dt * z1 + e2dt * z2 + e3dt * z3
        mass_matrix != I && (mul!(w1, mass_matrix, tmp); copyto!(tmp, w1))
        @.. broadcast=false ubuff=integrator.fsalfirst + tmp

        if alg.smooth_est
            linres1 = dolinsolve(integrator, linres1.cache; b = _vec(ubuff),
                linu = _vec(utilde))
            cache.linsolve1 = linres1.cache
            integrator.stats.nsolve += 1
        end

        # RadauIIA5 needs a transformed rtol and atol see
        # https://github.com/luchr/ODEInterface.jl/blob/0bd134a5a358c4bc13e0fb6a90e27e4ee79e0115/src/radau5.f#L399-L421
        calculate_residuals!(atmp, utilde, uprev, u, atol, rtol, internalnorm, t)
        integrator.EEst = internalnorm(atmp, t)

        if !(integrator.EEst < oneunit(integrator.EEst)) && integrator.iter == 1 ||
           integrator.u_modified
            @.. broadcast=false utilde=uprev + utilde
            f(fsallast, utilde, p, t)
            integrator.stats.nf += 1
            @.. broadcast=false ubuff=fsallast + tmp

            if alg.smooth_est
                linres1 = dolinsolve(integrator, linres1.cache; b = _vec(ubuff),
                    linu = _vec(utilde))
                cache.linsolve1 = linres1.cache
                integrator.stats.nsolve += 1
            end

            calculate_residuals!(atmp, utilde, uprev, u, atol, rtol, internalnorm, t)
            integrator.EEst = internalnorm(atmp, t)
        end
    end

    if integrator.EEst <= oneunit(integrator.EEst)
        cache.dtprev = dt
        if alg.extrapolant != :constant
            @.. broadcast=false cache.cont1=(z2 - z3) / c2m1
            @.. broadcast=false tmp=(z1 - z2) / c1mc2
            @.. broadcast=false cache.cont2=(tmp - cache.cont1) / c1m1
            @.. broadcast=false cache.cont3=cache.cont2 - (tmp - z1 / c1) / c2
        end
    end
    f(fsallast, u, p, t + dt)
    integrator.stats.nf += 1
    return
end

@muladd function perform_step!(integrator, cache::RadauIIA9ConstantCache,
        repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack T11, T12, T13, T14, T15, T21, T22, T23, T24, T25, T31, T32, T33, T34, T35, T41, T42, T43, T44, T45, T51 = cache.tab #= T52 = 1, T53 = 0, T54 = 1, T55 = 0=#
    @unpack TI11, TI12, TI13, TI14, TI15, TI21, TI22, TI23, TI24, TI25, TI31, TI32, TI33, TI34, TI35, TI41, TI42, TI43, TI44, TI45, TI51, TI52, TI53, TI54, TI55 = cache.tab
    @unpack c1, c2, c3, c4, γ, α1, β1, α2, β2, e1, e2, e3, e4, e5 = cache.tab
    @unpack κ, cont1, cont2, cont3, cont4 = cache
    @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
    alg = unwrap_alg(integrator, true)
    @unpack maxiters = alg
    mass_matrix = integrator.f.mass_matrix

    # precalculations rtol pow is (num stages + 1)/(2*num stages)
    rtol = @.. broadcast=false reltol^(5 / 8)/10
    atol = @.. broadcast=false rtol*(abstol / reltol)
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    c3m1 = c3 - 1
    c4m1 = c4 - 1
    c1mc2 = c1 - c2
    c1mc3 = c1 - c3
    c1mc4 = c1 - c4
    c2mc3 = c2 - c3
    c2mc4 = c2 - c4
    c3mc4 = c3 - c4

    γdt, α1dt, β1dt, α2dt, β2dt = γ / dt, α1 / dt, β1 / dt, α2 / dt, β2 / dt
    J = calc_J(integrator, cache)
    if u isa Number
        LU1 = -γdt * mass_matrix + J
        LU2 = -(α1dt + β1dt * im) * mass_matrix + J
        LU3 = -(α2dt + β2dt * im) * mass_matrix + J
    else
        LU1 = lu(-γdt * mass_matrix + J)
        LU2 = lu(-(α1dt + β1dt * im) * mass_matrix + J)
        LU3 = lu(-(α2dt + β2dt * im) * mass_matrix + J)
    end
    integrator.stats.nw += 1

    # TODO better initial guess
    if integrator.iter == 1 || integrator.u_modified || alg.extrapolant == :constant
        cache.dtprev = one(cache.dtprev)
        z1 = w1 = map(zero, u)
        z2 = w2 = map(zero, u)
        z3 = w3 = map(zero, u)
        z4 = w4 = map(zero, u)
        z5 = w5 = map(zero, u)
        cache.cont1 = map(zero, u)
        cache.cont2 = map(zero, u)
        cache.cont3 = map(zero, u)
        cache.cont4 = map(zero, u)
    else
        c5′ = dt / cache.dtprev
        c1′ = c1 * c5′
        c2′ = c2 * c5′
        c3′ = c3 * c5′
        c4′ = c4 * c5′
        z1 = @.. broadcast=false c1′*(cont1 +
                                      (c1′ - c3m1) * (cont2 +
                                       (c1′ - c2m1) * (cont3 + (c1′ - c1m1) * cont4)))
        z2 = @.. broadcast=false c2′*(cont1 +
                                      (c2′ - c3m1) * (cont2 +
                                       (c2′ - c2m1) * (cont3 + (c2′ - c1m1) * cont4)))
        z3 = @.. broadcast=false c3′*(cont1 +
                                      (c3′ - c3m1) * (cont2 +
                                       (c3′ - c2m1) * (cont3 + (c3′ - c1m1) * cont4)))
        z4 = @.. broadcast=false c4′*(cont1 +
                                      (c4′ - c3m1) * (cont2 +
                                       (c4′ - c2m1) * (cont3 + (c4′ - c1m1) * cont4)))
        z5 = @.. broadcast=false c5′*(cont1 +
                                      (c5′ - c3m1) * (cont2 +
                                       (c5′ - c2m1) * (cont3 + (c5′ - c1m1) * cont4)))
        w1 = @.. broadcast=false TI11*z1+TI12*z2+TI13*z3+TI14*z4+TI15*z5
        w2 = @.. broadcast=false TI21*z1+TI22*z2+TI23*z3+TI24*z4+TI25*z5
        w3 = @.. broadcast=false TI31*z1+TI32*z2+TI33*z3+TI34*z4+TI35*z5
        w4 = @.. broadcast=false TI41*z1+TI42*z2+TI43*z3+TI44*z4+TI45*z5
        w5 = @.. broadcast=false TI51*z1+TI52*z2+TI53*z3+TI54*z4+TI55*z5
    end

    # Newton iteration
    local ndw
    η = max(cache.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)
    fail_convergence = true
    iter = 0
    while iter < maxiters
        iter += 1
        integrator.stats.nnonliniter += 1

        # evaluate function
        ff1 = f(uprev + z1, p, t + c1 * dt)
        ff2 = f(uprev + z2, p, t + c2 * dt)
        ff3 = f(uprev + z3, p, t + c3 * dt)
        ff4 = f(uprev + z4, p, t + c4 * dt)
        ff5 = f(uprev + z5, p, t + dt) # c5 = 1
        integrator.stats.nf += 5

        fw1 = @.. broadcast=false TI11*ff1+TI12*ff2+TI13*ff3+TI14*ff4+TI15*ff5
        fw2 = @.. broadcast=false TI21*ff1+TI22*ff2+TI23*ff3+TI24*ff4+TI25*ff5
        fw3 = @.. broadcast=false TI31*ff1+TI32*ff2+TI33*ff3+TI34*ff4+TI35*ff5
        fw4 = @.. broadcast=false TI41*ff1+TI42*ff2+TI43*ff3+TI44*ff4+TI45*ff5
        fw5 = @.. broadcast=false TI51*ff1+TI52*ff2+TI53*ff3+TI54*ff4+TI55*ff5

        if mass_matrix isa UniformScaling # `UniformScaling` doesn't play nicely with broadcast
            Mw1 = @.. broadcast=false mass_matrix.λ*w1
            Mw2 = @.. broadcast=false mass_matrix.λ*w2
            Mw3 = @.. broadcast=false mass_matrix.λ*w3
            Mw4 = @.. broadcast=false mass_matrix.λ*w4
            Mw5 = @.. broadcast=false mass_matrix.λ*w5
        else
            Mw1 = mass_matrix * w1
            Mw2 = mass_matrix * w2
            Mw3 = mass_matrix * w3
            Mw4 = mass_matrix * w4
            Mw5 = mass_matrix * w5
        end

        rhs1 = @.. broadcast=false fw1-γdt * Mw1
        rhs2 = @.. broadcast=false fw2 - α1dt * Mw2+β1dt * Mw3
        rhs3 = @.. broadcast=false fw3 - β1dt * Mw2-α1dt * Mw3
        rhs4 = @.. broadcast=false fw4 - α2dt * Mw4+β2dt * Mw5
        rhs5 = @.. broadcast=false fw5 - β2dt * Mw4-α2dt * Mw5
        dw1 = LU1 \ rhs1
        dw23 = LU2 \ (@.. broadcast=false rhs2+rhs3 * im)
        dw45 = LU3 \ (@.. broadcast=false rhs4+rhs5 * im)
        integrator.stats.nsolve += 3
        dw2 = real(dw23)
        dw3 = imag(dw23)
        dw4 = real(dw45)
        dw5 = imag(dw45)

        # compute norm of residuals
        iter > 1 && (ndwprev = ndw)
        atmp1 = calculate_residuals(dw1, uprev, u, atol, rtol, internalnorm, t)
        atmp2 = calculate_residuals(dw2, uprev, u, atol, rtol, internalnorm, t)
        atmp3 = calculate_residuals(dw3, uprev, u, atol, rtol, internalnorm, t)
        atmp4 = calculate_residuals(dw4, uprev, u, atol, rtol, internalnorm, t)
        atmp5 = calculate_residuals(dw5, uprev, u, atol, rtol, internalnorm, t)
        ndw = internalnorm(atmp1, t) + internalnorm(atmp2, t) + internalnorm(atmp3, t) +
              internalnorm(atmp4, t) + internalnorm(atmp5, t)

        # check divergence (not in initial step)

        if iter > 1
            θ = ndw / ndwprev
            (diverge = θ > 1) && (cache.status = Divergence)
            (veryslowconvergence = ndw * θ^(maxiters - iter) > κ * (1 - θ)) &&
                (cache.status = VerySlowConvergence)
            if diverge || veryslowconvergence
                break
            end
        end

        w1 = @.. broadcast=false w1-dw1
        w2 = @.. broadcast=false w2-dw2
        w3 = @.. broadcast=false w3-dw3
        w4 = @.. broadcast=false w4-dw4
        w5 = @.. broadcast=false w5-dw5

        # transform `w` to `z`
        z1 = @.. broadcast=false T11*w1+T12*w2+T13*w3+T14*w4+T15*w5
        z2 = @.. broadcast=false T21*w1+T22*w2+T23*w3+T24*w4+T25*w5
        z3 = @.. broadcast=false T31*w1+T32*w2+T33*w3+T34*w4+T35*w5
        z4 = @.. broadcast=false T41*w1+T42*w2+T43*w3+T44*w4+T45*w5
        z5 = @.. broadcast=false T51*w1+w2+w4 #= T52=1, T53=0, T54=1, T55=0 =#

        # check stopping criterion
        iter > 1 && (η = θ / (1 - θ))
        if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
            # Newton method converges
            cache.status = η < alg.fast_convergence_cutoff ? FastConvergence :
                           Convergence
            fail_convergence = false
            break
        end
    end

    if fail_convergence
        integrator.force_stepfail = true
        integrator.stats.nnonlinconvfail += 1
        return
    end
    cache.ηold = η
    cache.iter = iter

    u = @.. broadcast=false uprev+z5

    if adaptive
        e1dt, e2dt, e3dt, e4dt, e5dt = e1 / dt, e2 / dt, e3 / dt, e4 / dt, e5 / dt
        tmp = @.. broadcast=false e1dt*z1+e2dt*z2+e3dt*z3+e4dt*z4+e5dt*z5
        mass_matrix != I && (tmp = mass_matrix * tmp)
        utilde = @.. broadcast=false integrator.fsalfirst+tmp
        alg.smooth_est && (utilde = LU1 \ utilde; integrator.stats.nsolve += 1)
        atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm, t)
        integrator.EEst = internalnorm(atmp, t)

        if !(integrator.EEst < oneunit(integrator.EEst)) && integrator.iter == 1 ||
           integrator.u_modified
            f0 = f(uprev .+ utilde, p, t)
            integrator.stats.nf += 1
            utilde = @.. broadcast=false f0+tmp
            alg.smooth_est && (utilde = LU1 \ utilde; integrator.stats.nsolve += 1)
            atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm, t)
            integrator.EEst = internalnorm(atmp, t)
        end
    end

    if integrator.EEst <= oneunit(integrator.EEst)
        cache.dtprev = dt
        if alg.extrapolant != :constant
            cache.cont1 = @.. (z4 - z5) / c4m1 # first derivative on [c4, 1]
            tmp1 = @.. (z3 - z4) / c3mc4 # first derivative on [c3, c4]
            cache.cont2 = @.. (tmp1 - cache.cont1) / c3m1 # second derivative on [c3, 1]
            tmp2 = @.. (z2 - z3) / c2mc3 # first derivative on [c2, c3]
            tmp3 = @.. (tmp2 - tmp1) / c2mc4 # second derivative on [c2, c4]
            cache.cont3 = @.. (tmp3 - cache.cont2) / c2m1 # third derivative on [c2, 1]
            tmp4 = @.. (z1 - z2) / c1mc2 # first derivative on [c1, c2]
            tmp5 = @.. (tmp4 - tmp2) / c1mc3 # second derivative on [c1, c3]
            tmp6 = @.. (tmp5 - tmp3) / c1mc4 # third derivative on [c1, c4]
            cache.cont4 = @.. (tmp6 - cache.cont3) / c1m1 #fourth derivative on [c1, 1]
        end
    end

    integrator.fsallast = f(u, p, t + dt)
    integrator.stats.nf += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return
end

@muladd function perform_step!(integrator, cache::RadauIIA9Cache, repeat_step = false)
    @unpack t, dt, uprev, u, f, p, fsallast, fsalfirst = integrator
    @unpack T11, T12, T13, T14, T15, T21, T22, T23, T24, T25, T31, T32, T33, T34, T35, T41, T42, T43, T44, T45, T51 = cache.tab #= T52 = 1, T53 = 0, T54 = 1, T55 = 0=#
    @unpack TI11, TI12, TI13, TI14, TI15, TI21, TI22, TI23, TI24, TI25, TI31, TI32, TI33, TI34, TI35, TI41, TI42, TI43, TI44, TI45, TI51, TI52, TI53, TI54, TI55 = cache.tab
    @unpack c1, c2, c3, c4, γ, α1, β1, α2, β2, e1, e2, e3, e4, e5 = cache.tab
    @unpack κ, cont1, cont2, cont3, cont4 = cache
    @unpack z1, z2, z3, z4, z5, w1, w2, w3, w4, w5 = cache
    @unpack dw1, ubuff, dw23, dw45, cubuff1, cubuff2 = cache
    @unpack k, k2, k3, k4, k5, fw1, fw2, fw3, fw4, fw5 = cache
    @unpack J, W1, W2, W3 = cache
    @unpack tmp, tmp2, tmp3, tmp4, tmp5, tmp6, atmp, jac_config, linsolve1, linsolve2, rtol, atol, step_limiter! = cache
    @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
    alg = unwrap_alg(integrator, true)
    @unpack maxiters = alg
    mass_matrix = integrator.f.mass_matrix

    # precalculations
    c1m1 = c1 - 1
    c2m1 = c2 - 1
    c3m1 = c3 - 1
    c4m1 = c4 - 1
    c1mc2 = c1 - c2
    c1mc3 = c1 - c3
    c1mc4 = c1 - c4
    c2mc3 = c2 - c3
    c2mc4 = c2 - c4
    c3mc4 = c3 - c4

    γdt, α1dt, β1dt, α2dt, β2dt = γ / dt, α1 / dt, β1 / dt, α2 / dt, β2 / dt
    (new_jac = do_newJ(integrator, alg, cache, repeat_step)) &&
        (calc_J!(J, integrator, cache); cache.W_γdt = dt)
    if (new_W = do_newW(integrator, alg, new_jac, cache.W_γdt))
        @inbounds for II in CartesianIndices(J)
            W1[II] = -γdt * mass_matrix[Tuple(II)...] + J[II]
            W2[II] = -(α1dt + β1dt * im) * mass_matrix[Tuple(II)...] + J[II]
            W3[II] = -(α2dt + β2dt * im) * mass_matrix[Tuple(II)...] + J[II]
        end
        integrator.stats.nw += 1
    end

    # TODO better initial guess
    if integrator.iter == 1 || integrator.u_modified || alg.extrapolant == :constant
        cache.dtprev = one(cache.dtprev)
        uzero = zero(eltype(u))
        @.. broadcast=false z1=uzero
        @.. broadcast=false z2=uzero
        @.. broadcast=false z3=uzero
        @.. broadcast=false z4=uzero
        @.. broadcast=false z5=uzero
        @.. broadcast=false w1=uzero
        @.. broadcast=false w2=uzero
        @.. broadcast=false w3=uzero
        @.. broadcast=false w4=uzero
        @.. broadcast=false w5=uzero
        @.. broadcast=false cache.cont1=uzero
        @.. broadcast=false cache.cont2=uzero
        @.. broadcast=false cache.cont3=uzero
        @.. broadcast=false cache.cont4=uzero
    else
        c5′ = dt / cache.dtprev
        c1′ = c1 * c5′
        c2′ = c2 * c5′
        c3′ = c3 * c5′
        c4′ = c4 * c5′
        z1 = @.. broadcast=false c1′*(cont1 +
                                      (c1′ - c3m1) * (cont2 +
                                       (c1′ - c2m1) * (cont3 + (c1′ - c1m1) * cont4)))
        z2 = @.. broadcast=false c2′*(cont1 +
                                      (c2′ - c3m1) * (cont2 +
                                       (c2′ - c2m1) * (cont3 + (c2′ - c1m1) * cont4)))
        z3 = @.. broadcast=false c3′*(cont1 +
                                      (c3′ - c3m1) * (cont2 +
                                       (c3′ - c2m1) * (cont3 + (c3′ - c1m1) * cont4)))
        z4 = @.. broadcast=false c4′*(cont1 +
                                      (c4′ - c3m1) * (cont2 +
                                       (c4′ - c2m1) * (cont3 + (c4′ - c1m1) * cont4)))
        z5 = @.. broadcast=false c5′*(cont1 +
                                      (c5′ - c3m1) * (cont2 +
                                       (c5′ - c2m1) * (cont3 + (c5′ - c1m1) * cont4)))
        w1 = @.. broadcast=false TI11*z1+TI12*z2+TI13*z3+TI14*z4+TI15*z5
        w2 = @.. broadcast=false TI21*z1+TI22*z2+TI23*z3+TI24*z4+TI25*z5
        w3 = @.. broadcast=false TI31*z1+TI32*z2+TI33*z3+TI34*z4+TI35*z5
        w4 = @.. broadcast=false TI41*z1+TI42*z2+TI43*z3+TI44*z4+TI45*z5
        w5 = @.. broadcast=false TI51*z1+TI52*z2+TI53*z3+TI54*z4+TI55*z5
    end

    # Newton iteration
    local ndw
    η = max(cache.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)
    fail_convergence = true
    iter = 0
    while iter < maxiters
        iter += 1
        integrator.stats.nnonliniter += 1

        # evaluate function
        @.. broadcast=false tmp=uprev + z1
        f(fsallast, tmp, p, t + c1 * dt)
        @.. broadcast=false tmp=uprev + z2
        f(k2, tmp, p, t + c2 * dt)
        @.. broadcast=false tmp=uprev + z3
        f(k3, tmp, p, t + c3 * dt)
        @.. broadcast=false tmp=uprev + z4
        f(k4, tmp, p, t + c4 * dt)
        @.. broadcast=false tmp=uprev + z5
        f(k5, tmp, p, t + dt) # c5 = 1
        integrator.stats.nf += 5

        @.. broadcast=false fw1=TI11 * fsallast + TI12 * k2 + TI13 * k3 + TI14 * k4 +
                                TI15 * k5
        @.. broadcast=false fw2=TI21 * fsallast + TI22 * k2 + TI23 * k3 + TI24 * k4 +
                                TI25 * k5
        @.. broadcast=false fw3=TI31 * fsallast + TI32 * k2 + TI33 * k3 + TI34 * k4 +
                                TI35 * k5
        @.. broadcast=false fw4=TI41 * fsallast + TI42 * k2 + TI43 * k3 + TI44 * k4 +
                                TI45 * k5
        @.. broadcast=false fw5=TI51 * fsallast + TI52 * k2 + TI53 * k3 + TI54 * k4 +
                                TI55 * k5

        if mass_matrix === I
            Mw1 = w1
            Mw2 = w2
            Mw3 = w3
            Mw4 = w4
            Mw5 = w5
        elseif mass_matrix isa UniformScaling
            mul!(z1, mass_matrix.λ, w1)
            mul!(z2, mass_matrix.λ, w2)
            mul!(z3, mass_matrix.λ, w3)
            mul!(z4, mass_matrix.λ, w4)
            mul!(z5, mass_matrix.λ, w5)
            Mw1 = z1
            Mw2 = z2
            Mw3 = z3
            Mw4 = z4
            Mw5 = z5
        else
            mul!(z1, mass_matrix, w1)
            mul!(z2, mass_matrix, w2)
            mul!(z3, mass_matrix, w3)
            mul!(z4, mass_matrix, w4)
            mul!(z5, mass_matrix, w5)
            Mw1 = z1
            Mw2 = z2
            Mw3 = z3
            Mw4 = z4
            Mw5 = z5
        end

        @.. broadcast=false ubuff=fw1 - γdt * Mw1
        needfactor = iter == 1 && new_W

        linsolve1 = cache.linsolve1

        if needfactor
            linres1 = dolinsolve(integrator, linsolve1; A = W1, b = _vec(ubuff),
                linu = _vec(dw1))
        else
            linres1 = dolinsolve(integrator, linsolve1; A = nothing, b = _vec(ubuff),
                linu = _vec(dw1))
        end

        cache.linsolve1 = linres1.cache

        @.. broadcast=false cubuff1=complex(
            fw2 - α1dt * Mw2 + β1dt * Mw3, fw3 - β1dt * Mw2 - α1dt * Mw3)

        linsolve2 = cache.linsolve2

        if needfactor
            linres2 = dolinsolve(integrator, linsolve2; A = W2, b = _vec(cubuff1),
                linu = _vec(dw23))
        else
            linres2 = dolinsolve(integrator, linsolve2; A = nothing, b = _vec(cubuff1),
                linu = _vec(dw23))
        end

        cache.linsolve2 = linres2.cache

        @.. broadcast=false cubuff2=complex(
            fw4 - α2dt * Mw4 + β2dt * Mw5, fw5 - β2dt * Mw4 - α2dt * Mw5)

        linsolve3 = cache.linsolve3

        if needfactor
            linres3 = dolinsolve(integrator, linsolve3; A = W3, b = _vec(cubuff2),
                linu = _vec(dw45))
        else
            linres3 = dolinsolve(integrator, linsolve3; A = nothing, b = _vec(cubuff2),
                linu = _vec(dw45))
        end

        cache.linsolve3 = linres3.cache

        integrator.stats.nsolve += 3
        dw2 = z2
        dw3 = z3
        @.. broadcast=false dw2=real(dw23)
        @.. broadcast=false dw3=imag(dw23)
        dw4 = z4
        dw5 = z5
        @.. broadcast=false dw4=real(dw45)
        @.. broadcast=false dw5=imag(dw45)

        # compute norm of residuals
        iter > 1 && (ndwprev = ndw)
        calculate_residuals!(atmp, dw1, uprev, u, atol, rtol, internalnorm, t)
        ndw1 = internalnorm(atmp, t)
        calculate_residuals!(atmp, dw2, uprev, u, atol, rtol, internalnorm, t)
        ndw2 = internalnorm(atmp, t)
        calculate_residuals!(atmp, dw3, uprev, u, atol, rtol, internalnorm, t)
        ndw3 = internalnorm(atmp, t)
        calculate_residuals!(atmp, dw4, uprev, u, atol, rtol, internalnorm, t)
        ndw4 = internalnorm(atmp, t)
        calculate_residuals!(atmp, dw5, uprev, u, atol, rtol, internalnorm, t)
        ndw5 = internalnorm(atmp, t)
        ndw = ndw1 + ndw2 + ndw3 + ndw4 + ndw5

        # check divergence (not in initial step)

        if iter > 1
            θ = ndw / ndwprev
            (diverge = θ > 1) && (cache.status = Divergence)
            (veryslowconvergence = ndw * θ^(maxiters - iter) > κ * (1 - θ)) &&
                (cache.status = VerySlowConvergence)
            if diverge || veryslowconvergence
                break
            end
        end

        @.. broadcast=false w1=w1 - dw1
        @.. broadcast=false w2=w2 - dw2
        @.. broadcast=false w3=w3 - dw3
        @.. broadcast=false w4=w4 - dw4
        @.. broadcast=false w5=w5 - dw5

        # transform `w` to `z`
        @.. broadcast=false z1=T11 * w1 + T12 * w2 + T13 * w3 + T14 * w4 + T15 * w5
        @.. broadcast=false z2=T21 * w1 + T22 * w2 + T23 * w3 + T24 * w4 + T25 * w5
        @.. broadcast=false z3=T31 * w1 + T32 * w2 + T33 * w3 + T34 * w4 + T35 * w5
        @.. broadcast=false z4=T41 * w1 + T42 * w2 + T43 * w3 + T44 * w4 + T45 * w5
        @.. broadcast=false z5=T51 * w1 + w2 + w4 #= T52=1, T53=0, T54=1, T55=0 =#

        # check stopping criterion

        iter > 1 && (η = θ / (1 - θ))
        if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
            # Newton method converges
            cache.status = η < alg.fast_convergence_cutoff ? FastConvergence :
                           Convergence
            fail_convergence = false
            break
        end
    end
    if fail_convergence
        integrator.force_stepfail = true
        integrator.stats.nnonlinconvfail += 1
        return
    end

    cache.ηold = η
    cache.iter = iter

    @.. broadcast=false u=uprev + z5

    step_limiter!(u, integrator, p, t + dt)

    if adaptive
        utilde = w2
        e1dt, e2dt, e3dt, e4dt, e5dt = e1 / dt, e2 / dt, e3 / dt, e4 / dt, e5 / dt
        @.. broadcast=false tmp=e1dt * z1 + e2dt * z2 + e3dt * z3 + e4dt * z4 + e5dt * z5
        mass_matrix != I && (mul!(w1, mass_matrix, tmp); copyto!(tmp, w1))
        @.. broadcast=false ubuff=integrator.fsalfirst + tmp

        if alg.smooth_est
            linres1 = dolinsolve(integrator, linres1.cache; b = _vec(ubuff),
                linu = _vec(utilde))
            cache.linsolve1 = linres1.cache
            integrator.stats.nsolve += 1
        end

        # RadauIIA5 needs a transformed rtol and atol see
        # https://github.com/luchr/ODEInterface.jl/blob/0bd134a5a358c4bc13e0fb6a90e27e4ee79e0115/src/radau5.f#L399-L421
        calculate_residuals!(atmp, utilde, uprev, u, atol, rtol, internalnorm, t)
        integrator.EEst = internalnorm(atmp, t)

        if !(integrator.EEst < oneunit(integrator.EEst)) && integrator.iter == 1 ||
           integrator.u_modified
            @.. broadcast=false utilde=uprev + utilde
            f(fsallast, utilde, p, t)
            integrator.stats.nf += 1
            @.. broadcast=false ubuff=fsallast + tmp

            if alg.smooth_est
                linres1 = dolinsolve(integrator, linres1.cache; b = _vec(ubuff),
                    linu = _vec(utilde))
                cache.linsolve1 = linres1.cache
                integrator.stats.nsolve += 1
            end

            calculate_residuals!(atmp, utilde, uprev, u, atol, rtol, internalnorm, t)
            integrator.EEst = internalnorm(atmp, t)
        end
    end

    if integrator.EEst <= oneunit(integrator.EEst)
        cache.dtprev = dt
        if alg.extrapolant != :constant
            @.. cache.cont1 = (z4 - z5) / c4m1
            @.. tmp = (z3 - z4) / c3mc4
            @.. cache.cont2 = (tmp - cache.cont1) / c3m1
            @.. tmp2 = (z2 - z3) / c2mc3
            @.. tmp3 = (tmp2 - tmp) / c2mc4
            @.. cache.cont3 = (tmp3 - cache.cont2) / c2m1
            @.. tmp4 = (z1 - z2) / c1mc2
            @.. tmp5 = (tmp4 - tmp2) / c1mc3
            @.. tmp6 = (tmp5 - tmp3) / c1mc4
            @.. cache.cont4 = (tmp6 - cache.cont3) / c1m1
        end
    end

    f(fsallast, u, p, t + dt)
    integrator.stats.nf += 1
    return
end

@muladd function perform_step!(integrator, cache::adaptiveRadauConstantCache,
    repeat_step = false, s::Int64)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack T, TI, γ, α, β, c, e= cache.tab 
    @unpack κ, cont = cache
    @unpack internalnorm, abstol, reltol, adaptive = integrator.opts
    alg = unwrap_alg(integrator, true)
    @unpack maxiters = alg
    mass_matrix = integrator.f.mass_matrix

    # precalculations rtol pow is (num stages + 1)/(2*num stages)
    rtol = @.. broadcast=false reltol^(5 / 8)/10
    atol = @.. broadcast=false rtol*(abstol / reltol)

    γdt, αdt, βdt = γ / dt, α ./ dt, β ./ dt

    J = calc_J(integrator, cache)
    LU = Vector{Any}(undef, (s + 1) / 2)
    if u isa Number
        LU[1] = -γdt * mass_matrix + J
        for i in 2 : (s + 1) / 2
            LU[i] = -(α[i - 1]dt + β[i - 1]dt * im) * mass_matrix + J
        end
    else
        LU[1] = lu(-γdt * mass_matrix + J)
        for i in 2 : (s + 1) / 2
            LU[i] = lu(-(α[i - 1]dt + β[i - 1]dt * im) * mass_matrix + J)
        end
    end
    integrator.stats.nw += 1

    if integrator.iter == 1 || integrator.u_modified || alg.extrapolant == :constant
        cache.dtprev = one(cache.dtprev)
        for i in 1:s
            z[i] = w[i] = map(zero, u)
        end
        for i in 1:s-1
            cont[i] = map(zero, u)
        end
    else 
        c' = Vector{eltype(u)}(undef, s) #time stepping
        c'[s] = dt / cache.dtprev
        for i in 1 : s-1
            c'[i] = c[i] * c'[s]
        end
        for i in 1 : s # collocation polynomial
            z[i] = @.. cont[s-1] * (c'[i] - c[1] + 1) + cont[s-1]
            j = s - 2
            while j > 0
                z[i] = @.. z[i] * (c'[i] - c[s-j] + 1) + cont[j]
            end
            z[i] = @.. z[i] * c'[i]
        end
        w = @.. TI * z
    end

    # Newton iteration
    local ndw
    η = max(cache.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)
    fail_convergence = true
    iter = 0
    while iter < maxiters
        iter += 1
        integrator.stats.nnonliniter += 1

        # evaluate function
        ff = Vector{eltype(u)}(undef, s)
        for i in 1:s
            ff[i] = f(uprev + z[i], p, t + c[i] * dt)
        end    
        integrator.stats.nf += 5

        fw = @.. TI * ff
        Mw = Vector{eltype(u)}(undef, s)
        if mass_matrix isa UniformScaling # `UniformScaling` doesn't play nicely with broadcast
            for i in 1:s
                Mw[i] = @.. mass_matrix.λ * w[i] #scaling by eigenvalues
            end
        else
            Mw = mass_matrix * w #standard multiplication
        end

        rhs = Vector{eltype(u)}(undef, s)
        rhs[1] = @.. fw[1]-γdt * Mw[1]
        i = 2
        while i <= s #block by block multiplication
            rhs[i] = @.. fw[i] - α[i/2]dt * Mw[i] + β[i/2]dt * Mw[i + 1]
            rhs[i + 1] = @.. fw[i + 1] - β[i/2]dt * Mw[i] - α[i/2]dt * Mw[i + 1]
            i += 2
        end 

        dw = Vector{eltype(u)}(undef, s)
        dw[1] = LU1 \ rhs[1]
        for i in 2 : (s + 1) / 2 
            tmp = LU[i] \ (@.. rhs[2 * i - 2] + rhs[2 * i - 1] * im)
            dw[2 * i - 2] = real(tmp)
            dw[2 * i - 1] = imag(tmp)
        end
        integrator.stats.nlsolve += (s + 1) / 2

        # compute norm of residuals
        iter > 1 && (ndwprev = ndw)
        atmp = Vector{eltype(u)}(undef, s)
        for i in 1:s
            atmp[i] = calculate_residuals(dw[i], uprev, u, atol, rtol, internalnorm, t)
        end

        ndw = 0
        for i in 1:s
            ndw = ndw + internalnorm(atmp[i], t)
        end
        # check divergence (not in initial step)

        if iter > 1
            θ = ndw / ndwprev
            (diverge = θ > 1) && (cache.status = Divergence)
            (veryslowconvergence = ndw * θ^(maxiters - iter) > κ * (1 - θ)) &&
                (cache.status = VerySlowConvergence)
            if diverge || veryslowconvergence
                break
            end
        end

        for i in 1 : s
            w[i] = @.. w[i] - dw[i]
        end
        # transform `w` to `z`
        z = @.. T * w
        
        # check stopping criterion
        iter > 1 && (η = θ / (1 - θ))
        if η * ndw < κ && (iter > 1 || iszero(ndw) || !iszero(integrator.success_iter))
            # Newton method converges
            cache.status = η < alg.fast_convergence_cutoff ? FastConvergence :
                        Convergence
            fail_convergence = false
            break
        end
    end

    if fail_convergence
        integrator.force_stepfail = true
        integrator.stats.nnonlinconvfail += 1
        return
    end
    cache.ηold = η
    cache.iter = iter

    u = @.. uprev + z[s]

    if adaptive
        edt = e ./ dt
        tmp = @.. dot(edt, z)
        mass_matrix != I && (tmp = mass_matrix * tmp)
        utilde = @.. broadcast=false integrator.fsalfirst+tmp
        alg.smooth_est && (utilde = LU[1] \ utilde; integrator.stats.nsolve += 1)
        atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm, t)
        integrator.EEst = internalnorm(atmp, t)

        if !(integrator.EEst < oneunit(integrator.EEst)) && integrator.iter == 1 ||
        integrator.u_modified
            f0 = f(uprev .+ utilde, p, t)
            integrator.stats.nf += 1
            utilde = @.. broadcast=false f0+tmp
            alg.smooth_est && (utilde = LU[1] \ utilde; integrator.stats.nsolve += 1)
            atmp = calculate_residuals(utilde, uprev, u, atol, rtol, internalnorm, t)
            integrator.EEst = internalnorm(atmp, t)
        end
    end

    if integrator.EEst <= oneunit(integrator.EEst)
        cache.dtprev = dt
        if alg.extrapolant != :constant
            derivatives = Matrix{eltype(u)}(undef, s-1, s-1)
            for i in 1 : (s - 1)
                for j in i : (s-1)
                    if i == 1
                        derivatives[i, j] = @.. (z[i] - z[i + 1]) / (c[i] - c[i + 1]) #first derivatives
                    else
                        derivatives[i, j] = @.. (derivatives[i - 1, j - 1] - derivatives[i - 1, j]) / (c[j - i + 1] - c[j + 1]) #all others
                    end
                end
            end
            for i in 1 : (s-1)
                cache.cont[i] = derivatives[i, i]
            end
        end
    end

    integrator.fsallast = f(u, p, t + dt)
    integrator.stats.nf += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
    return
end
