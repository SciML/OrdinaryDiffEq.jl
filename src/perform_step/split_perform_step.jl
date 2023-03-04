function initialize!(integrator, cache::SplitEulerConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f.f1(integrator.uprev, integrator.p, integrator.t) +
                           integrator.f.f2(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.destats.nf += 1
    integrator.destats.nf2 += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::SplitEulerConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, uprev, u, f, p) = integrator
    else
        @unpack t, dt, uprev, u, f, p = integrator
    end
    u = @.. broadcast=false uprev+dt * integrator.fsalfirst
    integrator.fsallast = f.f1(u, p, t + dt) + f.f2(u, p, t + dt)  # For the interpolation, needs k at the updated point
    integrator.destats.nf += 1
    integrator.destats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

function initialize!(integrator, cache::SplitEulerCache)
    integrator.kshortsize = 2
    @static if VERSION >= 1.8
        (; k, fsalfirst) = cache
    else
        @unpack k, fsalfirst = cache
    end
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    integrator.f.f2(cache.tmp, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    integrator.destats.nf += 1
    integrator.destats.nf2 += 1
    integrator.fsalfirst .+= cache.tmp
end

@muladd function perform_step!(integrator, cache::SplitEulerCache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, uprev, u, f, p) = integrator
    else
        @unpack t, dt, uprev, u, f, p = integrator
    end
    @.. broadcast=false u=uprev + dt * integrator.fsalfirst
    f.f1(integrator.fsallast, u, p, t + dt) # For the interpolation, needs k at the updated point
    f.f2(cache.tmp, u, p, t + dt) # For the interpolation, needs k at the updated point
    integrator.destats.nf2 += 1
    integrator.destats.nf += 1
    integrator.fsallast .+= cache.tmp
end
