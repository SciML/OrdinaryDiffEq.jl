function verify_f2(f, p, q, pa, t, integrator, ::BAOABConstantCache)
    res = f(p, q, pa, t)
    return res != p && throwex(integrator)
end
function verify_f2(f, res, p, q, pa, t, integrator, ::BAOABCache)
    f(res, p, q, pa, t)
    return res != p && throwex(integrator)
end
function throwex(integrator)
    algn = typeof(integrator.alg)
    throw(ArgumentError("Algorithm $algn is not applicable if f2(p, q, t) != p"))
end

function initialize!(integrator, cache::BAOABConstantCache)
    (; t, dt, uprev, u, p, W) = integrator
    du1 = integrator.uprev.x[1]
    u1 = integrator.uprev.x[2]

    verify_f2(integrator.f.f2, du1, u1, p, t, integrator, cache)
    return cache.k = integrator.f.f1(du1, u1, p, t)
end

function initialize!(integrator, cache::BAOABCache)
    (; t, dt, uprev, u, p, W) = integrator
    du1 = integrator.uprev.x[1]
    u1 = integrator.uprev.x[2]

    verify_f2(integrator.f.f2, cache.k, du1, u1, p, t, integrator, cache)
    return integrator.f.f1(cache.k, du1, u1, p, t)
end

@muladd function perform_step!(integrator, cache::BAOABConstantCache)
    (; t, dt, sqdt, uprev, u, p, W, f) = integrator
    (; half, c1, c2) = cache
    du1 = uprev.x[1]
    u1 = uprev.x[2]

    # B
    du2 = du1 + half * dt * cache.k

    # A
    u2 = u1 + half * dt * du2

    # O
    noise = integrator.f.g(u2, p, t + dt * half) .* W.dW
    du3 = c1 * du2 + c2 * noise

    # A
    u = u2 + half * dt * du3

    # B
    cache.k = f.f1(du3, u, p, t + dt)
    du = du3 + half * dt * cache.k

    integrator.u = ArrayPartition((du, u))
end

@muladd function perform_step!(integrator, cache::BAOABCache)
    (; t, dt, sqdt, uprev, u, p, W, f) = integrator
    (; utmp, dutmp, k, gtmp, noise, half, c1, c2) = cache
    du1 = uprev.x[1]
    u1 = uprev.x[2]

    # B
    @.. dutmp = du1 + half * dt * k

    # A
    @.. utmp = u1 + half * dt * dutmp

    # O
    integrator.f.g(gtmp, utmp, p, t + dt * half)
    @.. noise = gtmp * W.dW
    @.. dutmp = c1 * dutmp + c2 * noise

    # A
    @.. u.x[2] = utmp + half * dt * dutmp

    # B
    f.f1(k, dutmp, u.x[2], p, t + dt)
    @.. u.x[1] = dutmp + half * dt * k
end
