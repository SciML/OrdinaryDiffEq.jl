function initialize!(integrator, cache::SplitEulerConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f.f1(integrator.uprev, integrator.p, integrator.t) +
        integrator.f.f2(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(
        integrator, cache::SplitEulerConstantCache,
        repeat_step = false
    )
    (; t, dt, uprev, u, f, p) = integrator
    u = @.. broadcast = false uprev + dt * integrator.fsalfirst
    integrator.fsallast = f.f1(u, p, t + dt) + f.f2(u, p, t + dt)  # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

get_fsalfirstlast(cache::SplitEulerCache, u) = (cache.fsalfirst, cache.k)
function initialize!(integrator, cache::SplitEulerCache)
    integrator.kshortsize = 2
    (; k, fsalfirst) = cache
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.f1(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    integrator.f.f2(cache.tmp, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    return integrator.fsalfirst .+= cache.tmp
end

@muladd function perform_step!(integrator, cache::SplitEulerCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    @.. broadcast = false u = uprev + dt * integrator.fsalfirst
    f.f1(integrator.fsallast, u, p, t + dt) # For the interpolation, needs k at the updated point
    f.f2(cache.tmp, u, p, t + dt) # For the interpolation, needs k at the updated point
    integrator.stats.nf2 += 1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast .+= cache.tmp
end
