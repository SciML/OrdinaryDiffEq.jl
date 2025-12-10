function initialize!(integrator, cache::AN5ConstantCache)
    integrator.kshortsize = 7
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::AN5ConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, differential_vars) = integrator
    (; z, l, m, c_LTE, dts, tsit5tab) = cache
    # handle callbacks, rewind back to order one.
    if integrator.u_modified
        cache.order = 1
    end
    # Nordsieck form needs to build the history vector
    if cache.order == 1
        # Start the Nordsieck vector in one shot!
        perform_step!(integrator, tsit5tab, repeat_step)
        cache.order = 4
        z[1] = integrator.uprev
        z[2] = integrator.k[1] * dt
        z[3] = ode_interpolant(t, dt, nothing, nothing, integrator.k, tsit5tab, nothing,
            Val{2}, differential_vars) * dt^2 / 2
        z[4] = ode_interpolant(t, dt, nothing, nothing, integrator.k, tsit5tab, nothing,
            Val{3}, differential_vars) * dt^3 / 6
        z[5] = ode_interpolant(t, dt, nothing, nothing, integrator.k, tsit5tab, nothing,
            Val{4}, differential_vars) * dt^4 / 24
        z[6] = zero(cache.z[6])
        fill!(dts, dt)
        perform_predict!(cache)
        cache.Δ = integrator.u - integrator.uprev
        update_nordsieck_vector!(cache)
        if integrator.opts.adaptive && integrator.EEst >= one(integrator.EEst)
            cache.order = 1
        end
    else
        # Reset time
        tmp = dts[6]
        for i in 5:-1:1
            dts[i + 1] = dts[i]
        end
        dts[1] = dt
        dt != dts[2] && nordsieck_rescale!(cache)
        integrator.k[1] = z[2] / dt
        # Perform 5th order Adams method in Nordsieck form
        perform_predict!(cache)
        calc_coeff!(cache)
        isucceed = nlsolve_functional!(integrator, cache)
        if !isucceed
            # rewind Nordsieck vector
            integrator.force_stepfail = true
            nordsieck_rewind!(cache)
            return nothing
        end

        ################################### Error estimation

        if integrator.opts.adaptive
            atmp = calculate_residuals(
                cache.Δ, uprev, integrator.u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t)
            integrator.EEst = integrator.opts.internalnorm(atmp, t) * cache.c_LTE
            if integrator.EEst > one(integrator.EEst)
                for i in 1:5
                    dts[i] = dts[i + 1]
                end
                dts[6] = tmp
            end
        end

        # Correct Nordsieck vector
        cache.order = 5
        update_nordsieck_vector!(cache)

        ################################### Finalize

        integrator.k[2] = cache.z[2] / dt
    end
    return nothing
end

function initialize!(integrator, cache::AN5Cache)
    integrator.kshortsize = 7
    resize!(integrator.k, integrator.kshortsize)
    # Setup k pointers
    integrator.k[1] = cache.tsit5cache.k1
    integrator.k[2] = cache.tsit5cache.k2
    integrator.k[3] = cache.tsit5cache.k3
    integrator.k[4] = cache.tsit5cache.k4
    integrator.k[5] = cache.tsit5cache.k5
    integrator.k[6] = cache.tsit5cache.k6
    integrator.k[7] = cache.tsit5cache.k7
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::AN5Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p, uprev2, differential_vars) = integrator
    (; z, l, m, c_LTE, dts, tmp, ratetmp, atmp, tsit5cache) = cache
    # handle callbacks, rewind back to order one.
    if integrator.u_modified
        cache.order = 1
    end
    # Nordsieck form needs to build the history vector
    if cache.order == 1
        ## Start the Nordsieck vector in two shots!
        perform_step!(integrator, tsit5cache, repeat_step)
        copyto!(tmp, integrator.u)
        cache.order = 4
        @.. broadcast=false z[1]=integrator.uprev
        @.. broadcast=false z[2]=integrator.k[1] * dt
        ode_interpolant!(z[3], t, dt, nothing, nothing, integrator.k, tsit5cache, nothing,
            Val{2}, differential_vars)
        ode_interpolant!(z[4], t, dt, nothing, nothing, integrator.k, tsit5cache, nothing,
            Val{3}, differential_vars)
        ode_interpolant!(z[5], t, dt, nothing, nothing, integrator.k, tsit5cache, nothing,
            Val{4}, differential_vars)
        @.. broadcast=false z[3]=z[3] * dt^2 / 2
        @.. broadcast=false z[4]=z[4] * dt^3 / 6
        @.. broadcast=false z[5]=z[5] * dt^4 / 24
        fill!(z[6], 0)
        fill!(dts, dt)
        perform_predict!(cache)
        @.. broadcast=false cache.Δ=integrator.u - integrator.uprev
        update_nordsieck_vector!(cache)
        if integrator.opts.adaptive && integrator.EEst >= one(integrator.EEst)
            cache.order = 1
        end
    else
        # Reset time
        tmp = dts[6]
        for i in 5:-1:1
            dts[i + 1] = dts[i]
        end
        dts[1] = dt
        # Rescale
        dt != dts[2] && nordsieck_rescale!(cache)
        @.. broadcast=false integrator.k[1]=z[2] / dt
        # Perform 5th order Adams method in Nordsieck form
        perform_predict!(cache)
        calc_coeff!(cache)
        isucceed = nlsolve_functional!(integrator, cache)
        if !isucceed
            integrator.force_stepfail = true
            # rewind Nordsieck vector
            nordsieck_rewind!(cache)
            return nothing
        end

        ################################### Error estimation

        if integrator.opts.adaptive
            calculate_residuals!(
                atmp, cache.Δ, uprev, integrator.u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t)
            integrator.EEst = integrator.opts.internalnorm(atmp, t) * cache.c_LTE
            if integrator.EEst > one(integrator.EEst)
                for i in 1:5
                    dts[i] = dts[i + 1]
                end
                dts[6] = tmp
            end
        end

        # Correct Nordsieck vector
        cache.order = 5
        update_nordsieck_vector!(cache)

        ################################### Finalize

        @.. broadcast=false integrator.k[2]=cache.z[2] / dt
    end
    return nothing
end

function initialize!(integrator, cache::JVODEConstantCache)
    integrator.kshortsize = 7
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::JVODEConstantCache, repeat_step = false)
    (; t, dt, uprev, u, f, p, differential_vars) = integrator
    (; z, l, m, c_LTE, dts, tsit5tab) = cache
    # handle callbacks, rewind back to order one.
    if integrator.u_modified || integrator.iter == 1
        cache.order = 1
        z[1] = integrator.uprev
        z[2] = f(uprev, p, t) * dt
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        dts[1] = dt
    end
    # Reset time
    tmp = dts[13]
    for i in 12:-1:1
        dts[i + 1] = dts[i]
    end
    dts[1] = dt
    dt != dts[2] && nordsieck_adjust!(integrator, cache)
    integrator.k[1] = z[2] / dt

    perform_predict!(cache)
    calc_coeff!(cache)
    isucceed = nlsolve_functional!(integrator, cache)
    if !isucceed
        # rewind Nordsieck vector
        integrator.force_stepfail = true
        nordsieck_rewind!(cache)
        return nothing
    end

    ################################### Error estimation
    if integrator.opts.adaptive
        atmp = calculate_residuals(cache.Δ, uprev, integrator.u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t) * cache.c_LTE
        if integrator.EEst > one(integrator.EEst)
            for i in 1:12
                dts[i] = dts[i + 1]
            end
            dts[13] = tmp
        end
    end

    ################################### Finalize
    nordsieck_finalize!(integrator, cache)
    nordsieck_prepare_next!(integrator, cache)
    integrator.k[2] = cache.z[2] / dt
    return nothing
end

function initialize!(integrator, cache::JVODECache)
    integrator.kshortsize = 7
    resize!(integrator.k, integrator.kshortsize)
    # Setup k pointers
    integrator.k[1] = cache.tsit5cache.k1
    integrator.k[2] = cache.tsit5cache.k2
    integrator.k[3] = cache.tsit5cache.k3
    integrator.k[4] = cache.tsit5cache.k4
    integrator.k[5] = cache.tsit5cache.k5
    integrator.k[6] = cache.tsit5cache.k6
    integrator.k[7] = cache.tsit5cache.k7
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

@muladd function perform_step!(integrator, cache::JVODECache, repeat_step = false)
    (; t, dt, uprev, u, f, p, uprev2, differential_vars) = integrator
    (; z, l, m, c_LTE, dts, tmp, ratetmp, atmp, tsit5cache) = cache
    # handle callbacks, rewind back to order one.
    if integrator.u_modified || integrator.iter == 1
        cache.order = 1
        @.. broadcast=false z[1]=integrator.uprev
        f(z[2], uprev, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast=false z[2]=z[2] * dt
        dts[1] = dt
    end
    # Reset time
    tmp = dts[13]
    for i in 12:-1:1
        dts[i + 1] = dts[i]
    end
    dts[1] = dt
    # Rescale
    dt != dts[2] && nordsieck_adjust!(integrator, cache)
    @.. broadcast=false integrator.k[1]=z[2] / dt

    perform_predict!(cache)
    calc_coeff!(cache)
    isucceed = nlsolve_functional!(integrator, cache)
    # TODO: Handle NLsolve better
    if !isucceed
        integrator.force_stepfail = true
        # rewind Nordsieck vector
        nordsieck_rewind!(cache)
        return nothing
    end

    ################################### Error estimation

    if integrator.opts.adaptive
        calculate_residuals!(atmp, cache.Δ, uprev, integrator.u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t) * cache.c_LTE
        if integrator.EEst > one(integrator.EEst)
            for i in 1:12
                dts[i] = dts[i + 1]
            end
            dts[13] = tmp
        end
    end

    ################################### Finalize

    nordsieck_finalize!(integrator, cache)
    nordsieck_prepare_next!(integrator, cache)
    @.. broadcast=false integrator.k[2]=cache.z[2] / dt
    return nothing
end
