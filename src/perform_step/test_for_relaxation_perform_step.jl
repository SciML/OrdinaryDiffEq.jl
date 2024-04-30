function initialize!(integrator, cache::Tsit5ConstantCache_for_relaxation)
    integrator.kshortsize = 7
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.stats.nf += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:(integrator.kshortsize - 1)
        integrator.k[i] = zero(integrator.fsalfirst)
    end
    integrator.k[integrator.kshortsize] = integrator.fsallast
end


function perform_step!(integrator, cache::Tsit5ConstantCache_for_relaxation, repeat_step = false)

    # Variable to know if dt has changed during perform_step
    integrator.dt_has_changed = false

    # computations! will only contain the mathematical scheme
    # i.e the computations of the u(t+dt)
    # the result is store not in integrator.u but integrator.u_propose
    computations!(integrator, cache, repeat_step)

    # modif_step! enables to modify the step like when we want to perform a relaxation
    # for this we give a new struture that can be defined either by us for already known
    # modification we want to do or by a user (see below)
    modif_step!(integrator)

    # finalize_step! will do staff related to the solver like integrator.stats, register integrator.fsal
    # and register integrator.u
    finalize_step!(integrator, cache)
end


@muladd function computations!(integrator, cache::Tsit5ConstantCache_for_relaxation, repeat_step = false)
    @unpack t, dt, uprev, u, f, p = integrator
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract Tsit5ConstantCacheActual T T2
    k1 = integrator.fsalfirst
    a = dt * a21
    k2 = f(uprev + a * k1, p, t + c1 * dt)
    k3 = f(uprev + dt * (a31 * k1 + a32 * k2), p, t + c2 * dt)
    k4 = f(uprev + dt * (a41 * k1 + a42 * k2 + a43 * k3), p, t + c3 * dt)
    k5 = f(uprev + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4), p, t + c4 * dt)
    g6 = uprev + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = f(g6, p, t + dt)
    u = uprev + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
    k7 = f(u, p, t + dt)
    integrator.k[1] = k1
    integrator.k[2] = k2
    integrator.k[3] = k3
    integrator.k[4] = k4
    integrator.k[5] = k5
    integrator.k[6] = k6
    integrator.k[7] = k7
    integrator.u_propose = u
end


function modif_step!(integrator)
    
    # Perform the modifications
    #integrator.modif(integrator)

    # Here we check the validity of chaging dt if it has changed
    # if it is valid integrator.changed_valid will be true, if not it will be false
    changed_valid = true
    if integrator.dt_has_changed
        # check dt in [dtmin, dtmax]
        # things related to tstops
        # surely other things
        if changed_valid
            integrator.u_propose = integrator.u_changed
            integrator.dt = integrator.dt_changed
        else
            # print error or warning
        end
    end
end


function finalize_step!(integrator, cache::Tsit5ConstantCache_for_relaxation)
    @unpack t, dt, uprev, u, u_propose, f, p = integrator
    integrator.u = u_propose
    integrator.fsallast = f(u_propose, p, t + dt)

    # We ask the Tableau again but this might be improve by option to ask only the wanted part of the tableau
    T = constvalue(recursive_unitless_bottom_eltype(u))
    T2 = constvalue(typeof(one(t)))
    @OnDemandTableauExtract Tsit5ConstantCacheActual T T2

    if integrator.opts.adaptive
        utilde = dt *
                 (btilde1 * integrator.k[1] + btilde2 * integrator.k[2] + btilde3 * integrator.k[3] + btilde4 * integrator.k[4] + btilde5 * integrator.k[5] +
                  btilde6 * integrator.k[6] + btilde7 * integrator.k[7])
        atmp = calculate_residuals(utilde, uprev, u_propose, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    integrator.stats.nf += 6
end