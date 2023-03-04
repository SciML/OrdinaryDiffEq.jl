# http://www.chimica.unipd.it/antonino.polimeno/pubblica/downloads/JChemPhys_101_4062.pdf

function initialize!(integrator, cache::SymplecticEulerConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    # Do the calculation pre
    # So that way FSAL interpolation
    duprev, uprev = integrator.uprev.x
    du, u = integrator.u.x
    kdu = integrator.f.f1(duprev, uprev, integrator.p, integrator.t)
    kuprev = integrator.f.f2(duprev, uprev, integrator.p, integrator.t)
    @muladd du = duprev + integrator.dt * kdu
    ku = integrator.f.f2(du, uprev, integrator.p, integrator.t)
    integrator.destats.nf2 += 1
    integrator.destats.nf += 2
    integrator.fsalfirst = ArrayPartition((kdu, kuprev))
    integrator.fsallast = ArrayPartition((zero(kdu), ku))
end

@muladd function perform_step!(integrator, cache::SymplecticEulerConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    duprev, uprev = integrator.uprev.x
    kuprev = integrator.fsalfirst.x[2]
    u = uprev + dt * kuprev
    # Now actually compute the step
    # Do it at the end for interpolations!
    kdu = f.f1(duprev, u, p, t)
    du = duprev + dt * kdu

    ku = f.f2(du, u, p, t)
    integrator.destats.nf2 += 1
    integrator.destats.nf += 1

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((kdu, ku))
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::SymplecticEulerCache)
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
    # Do the calculation pre
    # So that way FSAL interpolation
    duprev, uprev = integrator.uprev.x
    du, u = integrator.u.x
    kuprev = integrator.fsalfirst.x[2]
    kdu, ku = integrator.fsallast.x
    integrator.f.f1(kdu, duprev, uprev, integrator.p, integrator.t)
    integrator.f.f2(kuprev, duprev, uprev, integrator.p, integrator.t)
    @muladd @.. broadcast=false du=duprev + integrator.dt * kdu
    integrator.f.f2(ku, du, uprev, integrator.p, integrator.t)
    integrator.destats.nf += 1
    integrator.destats.nf2 += 2
end

@muladd function perform_step!(integrator, cache::SymplecticEulerCache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    duprev, uprev = integrator.uprev.x
    du, u = integrator.u.x
    kuprev = integrator.fsalfirst.x[2]
    kdu, ku = integrator.fsallast.x
    @.. broadcast=false u=uprev + dt * kuprev
    # Now actually compute the step
    # Do it at the end for interpolations!
    integrator.destats.nf2 += 1
    integrator.destats.nf += 1
    f.f1(kdu, duprev, u, p, t)
    @.. broadcast=false du=duprev + dt * kdu
    f.f2(ku, du, u, p, t)
end

const MutableCachesHamilton = Union{Symplectic2Cache, Symplectic3Cache,
                                    Symplectic4Cache, Symplectic45Cache, Symplectic5Cache,
                                    Symplectic6Cache, Symplectic62Cache,
                                    McAte8Cache, KahanLi8Cache, SofSpa10Cache}
const MutableCachesNewton = Union{VelocityVerletCache}

const ConstantCachesHamilton = Union{Symplectic2ConstantCache, Symplectic3ConstantCache,
                                     Symplectic4ConstantCache, Symplectic45ConstantCache,
                                     Symplectic5ConstantCache,
                                     Symplectic6ConstantCache, Symplectic62ConstantCache,
                                     McAte8ConstantCache, KahanLi8ConstantCache,
                                     SofSpa10ConstantCache}
const ConstantCachesNewton = Union{VelocityVerletConstantCache}

# some of the algorithms are designed only for the case
# f.f2(p, q, pa, t) = p which is the Newton/Lagrange equations
# If called with different functions (which are possible in the Hamiltonian case)
# an exception is thrown to avoid silently calculate wrong results.
verify_f2(f, p, q, pa, t, ::Any, ::C) where {C <: ConstantCachesHamilton} = f(p, q, pa, t)
function verify_f2(f, res, p, q, pa, t, ::Any, ::C) where {C <: MutableCachesHamilton}
    f(res, p, q, pa, t)
end

function verify_f2(f, p, q, pa, t, integrator, ::C) where {C <: ConstantCachesNewton}
    res = f(p, q, pa, t)
    res == p ? p : throwex(integrator)
end
function verify_f2(f, res, p, q, pa, t, integrator, ::C) where {C <: MutableCachesNewton}
    f(res, p, q, pa, t)
    res == p ? res : throwex(integrator)
end
function throwex(integrator)
    algn = typeof(integrator.alg)
    throw(ArgumentError("Algorithm $algn is not applicable if f2(p, q, t) != p"))
end

# provide the mutable uninitialized objects to keep state and derivative in case of mutable caches
# no such objects are required for constant caches
function alloc_symp_state(integrator)
    (integrator.u.x..., integrator.cache.tmp.x...)
end

# load state and derivatives at begin of symplectic iteration steps
function load_symp_state(integrator)
    (integrator.uprev.x..., integrator.fsallast.x...)
end

# store state and derivatives at the end of symplectic iteration steps
function store_symp_state!(integrator, ::OrdinaryDiffEqConstantCache, du, u, kdu, ku)
    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((kdu, ku))
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    nothing
end

function store_symp_state!(integrator, ::OrdinaryDiffEqMutableCache, kdu, ku)
    copyto!(integrator.k[1].x[1], integrator.k[2].x[1])
    copyto!(integrator.k[1].x[2], integrator.k[2].x[2])
    copyto!(integrator.k[2].x[2], ku)
    copyto!(integrator.k[2].x[1], kdu)
    nothing
end

function initialize!(integrator,
                     cache::C) where {C <:
                                      Union{MutableCachesHamilton, MutableCachesNewton}}
    integrator.fsalfirst = cache.fsalfirst
    integrator.fsallast = cache.k

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    duprev, uprev = integrator.uprev.x
    integrator.f.f1(integrator.k[2].x[1], duprev, uprev, integrator.p, integrator.t)
    verify_f2(integrator.f.f2, integrator.k[2].x[2], duprev, uprev, integrator.p,
              integrator.t, integrator, cache)
    integrator.destats.nf += 1
    integrator.destats.nf2 += 1
end

function initialize!(integrator,
                     cache::C) where {
                                      C <:
                                      Union{ConstantCachesHamilton, ConstantCachesNewton}}
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    duprev, uprev = integrator.uprev.x
    kdu = integrator.f.f1(duprev, uprev, integrator.p, integrator.t)
    ku = verify_f2(integrator.f.f2, duprev, uprev, integrator.p, integrator.t, integrator,
                   cache)
    integrator.destats.nf += 1
    integrator.destats.nf2 += 1
    integrator.fsallast = ArrayPartition((kdu, ku))
    integrator.k[2] = integrator.fsallast
    integrator.fsalfirst = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::VelocityVerletConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    duprev, uprev = load_symp_state(integrator)

    # x(t+Δt) = x(t) + v(t)*Δt + 1/2*a(t)*Δt^2
    ku = integrator.fsallast.x[1]
    dtsq = dt^2
    half = cache.half
    u = uprev + dt * duprev + dtsq * (half * ku)
    kdu = f.f1(duprev, u, p, t + dt)
    # v(t+Δt) = v(t) + 1/2*(a(t)+a(t+Δt))*Δt
    du = duprev + dt * (half * ku + half * kdu)

    integrator.destats.nf += 2
    store_symp_state!(integrator, cache, du, u, kdu, du)
end

@muladd function perform_step!(integrator, cache::VelocityVerletCache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    duprev, uprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # x(t+Δt) = x(t) + v(t)*Δt + 1/2*a(t)*Δt^2
    ku = integrator.fsallast.x[1]
    dtsq = dt^2
    half = cache.half
    @.. broadcast=false u=uprev + dt * duprev + dtsq * (half * ku)
    f.f1(kdu, duprev, u, p, t + dt)
    integrator.destats.nf += 2
    # v(t+Δt) = v(t) + 1/2*(a(t)+a(t+Δt))*Δt
    @.. broadcast=false du=duprev + dt * (half * ku + half * kdu)

    store_symp_state!(integrator, cache, kdu, du)
end

@muladd function perform_step!(integrator, cache::Symplectic2ConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, b1, b2) = cache
    else
        @unpack a1, a2, b1, b2 = cache
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)

    # update position
    u = uprev + dt * b1 * kuprev
    # update velocity
    kdu = f.f1(duprev, u, p, t)
    du = duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b2 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a2 * kdu
    kdu = f.f1(du, u, p, tnew)
    ku = f.f2(du, u, p, tnew)

    integrator.destats.nf += 3
    integrator.destats.nf2 += 2
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic2Cache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, b1, b2) = cache.tab
    else
        @unpack a1, a2, b1, b2 = cache.tab
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # update position
    @.. broadcast=false u=uprev + dt * b1 * kuprev
    # update velocity
    f.f1(kdu, duprev, u, p, t)
    @.. broadcast=false du=duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b2 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a2 * kdu
    f.f1(kdu, du, u, p, tnew)
    f.f2(ku, du, u, p, tnew)

    integrator.destats.nf += 3
    integrator.destats.nf2 += 2
    store_symp_state!(integrator, cache, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic3ConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, a3, b1, b2, b3) = cache
    else
        @unpack a1, a2, a3, b1, b2, b3 = cache
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)

    # update position
    u = uprev + dt * b1 * kuprev
    # update velocity
    kdu = f.f1(duprev, u, p, integrator.t)
    du = duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b2 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b3 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a3 * kdu
    kdu = f.f1(du, u, p, tnew)
    ku = f.f2(du, u, p, tnew)

    integrator.destats.nf += 4
    integrator.destats.nf2 += 3
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic3Cache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, a3, b1, b2, b3) = cache.tab
    else
        @unpack a1, a2, a3, b1, b2, b3 = cache.tab
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # update position
    @.. broadcast=false u=uprev + dt * b1 * kuprev
    # update velocity
    f.f1(kdu, duprev, u, p, integrator.t)
    @.. broadcast=false du=duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b2 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b3 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a3 * kdu
    f.f1(kdu, du, u, p, tnew)
    f.f2(ku, du, u, p, tnew)

    integrator.destats.nf += 4
    integrator.destats.nf2 += 3
    store_symp_state!(integrator, cache, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic4ConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, a3, a4, b1, b2, b3, b4) = cache
    else
        @unpack a1, a2, a3, a4, b1, b2, b3, b4 = cache
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)

    # update position
    u = uprev + dt * b1 * kuprev
    # update velocity
    kdu = f.f1(duprev, u, p, t)
    du = duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b2 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b3 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b4 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a4 * kdu
    kdu = f.f1(du, u, p, tnew)
    ku = f.f2(du, u, p, tnew)

    integrator.destats.nf += 5
    integrator.destats.nf2 += 4
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic4Cache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, a3, a4, b1, b2, b3, b4) = cache.tab
    else
        @unpack a1, a2, a3, a4, b1, b2, b3, b4 = cache.tab
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # update position
    @.. broadcast=false u=uprev + dt * b1 * kuprev
    # update velocity
    f.f1(kdu, duprev, u, p, t)
    @.. broadcast=false du=duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b2 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b3 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b4 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a4 * kdu
    f.f1(kdu, du, u, p, tnew)
    f.f2(ku, du, u, p, tnew)

    integrator.destats.nf += 5
    integrator.destats.nf2 += 4
    store_symp_state!(integrator, cache, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic45ConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    alg = unwrap_alg(integrator, false)
    @static if VERSION >= 1.8
        (; a1, a2, a3, a4, a5, b1, b2, b3, b4, b5) = cache
    else
        @unpack a1, a2, a3, a4, a5, b1, b2, b3, b4, b5 = cache
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)

    # update position
    u = uprev + dt * b1 * kuprev
    # update velocity
    kdu = f.f1(duprev, u, p, t)
    du = duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b2 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b3 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b4 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b5 * ku

    kdu = f.f1(du, u, p, tnew)
    if typeof(alg) <: McAte42
        du = du + dt * a5 * kdu
        kdu = f.f1(du, u, p, tnew)
        integrator.destats.nf += 1
    end
    ku = f.f2(du, u, p, tnew)

    integrator.destats.nf += 5
    integrator.destats.nf2 += 5
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic45Cache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    alg = unwrap_alg(integrator, false)
    @static if VERSION >= 1.8
        (; a1, a2, a3, a4, a5, b1, b2, b3, b4, b5) = cache.tab
    else
        @unpack a1, a2, a3, a4, a5, b1, b2, b3, b4, b5 = cache.tab
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # update position
    @.. broadcast=false u=uprev + dt * b1 * kuprev
    # update velocity
    f.f1(kdu, duprev, u, p, t)
    @.. broadcast=false du=duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b2 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b3 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b4 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b5 * ku

    f.f1(kdu, du, u, p, tnew)
    if typeof(alg) <: McAte42
        @.. broadcast=false du=du + dt * a5 * kdu
        f.f1(kdu, du, u, p, tnew)
        integrator.destats.nf += 1
    end
    f.f2(ku, du, u, p, tnew)

    integrator.destats.nf += 5
    integrator.destats.nf2 += 5
    store_symp_state!(integrator, cache, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic5ConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6) = cache
    else
        @unpack a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6 = cache
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)

    # update position
    u = uprev + dt * b1 * kuprev
    # update velocity
    kdu = f.f1(duprev, u, p, integrator.t)
    du = duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b2 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b3 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b4 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b5 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a5 * kdu

    tnew = tnew + t + a5 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b6 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a6 * kdu
    kdu = f.f1(du, u, p, tnew)
    ku = f.f2(du, u, p, tnew)

    integrator.destats.nf += 7
    integrator.destats.nf2 += 6
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic5Cache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6) = cache.tab
    else
        @unpack a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6 = cache.tab
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # update position
    @.. broadcast=false u=uprev + dt * b1 * kuprev
    # update velocity
    f.f1(kdu, duprev, u, p, integrator.t)
    @.. broadcast=false du=duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b2 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b3 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b4 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b5 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a5 * kdu

    tnew = tnew + t + a5 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b6 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a6 * kdu
    f.f1(kdu, du, u, p, tnew)
    f.f2(ku, du, u, p, tnew)

    integrator.destats.nf += 7
    integrator.destats.nf2 += 6
    store_symp_state!(integrator, cache, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic6ConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8) = cache
    else
        @unpack a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8 = cache
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)

    # update position
    u = uprev + dt * b1 * kuprev
    # update velocity
    kdu = f.f1(duprev, u, p, integrator.t)
    du = duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b2 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b3 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b4 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b5 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a5 * kdu

    tnew = tnew + a5 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b6 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a6 * kdu

    tnew = tnew + a6 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b7 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a7 * kdu

    tnew = tnew + a7 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b8 * ku

    kdu = f.f1(du, u, p, tnew)
    # @.. broadcast=false du = du + dt*a8*kdu
    ku = f.f2(du, u, p, tnew)

    integrator.destats.nf += 8
    integrator.destats.nf2 += 8
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic6Cache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8) = cache.tab
    else
        @unpack a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8 = cache.tab
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # update position
    @.. broadcast=false u=uprev + dt * b1 * kuprev
    # update velocity
    f.f1(kdu, duprev, u, p, integrator.t)
    @.. broadcast=false du=duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b2 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b3 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b4 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b5 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a5 * kdu

    tnew = tnew + a5 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b6 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a6 * kdu

    tnew = tnew + a6 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b7 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a7 * kdu

    tnew = tnew + a7 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b8 * ku

    f.f1(kdu, du, u, p, tnew)
    # @.. broadcast=false du = du + dt*a8*kdu
    f.f2(ku, du, u, p, tnew)

    integrator.destats.nf += 8
    integrator.destats.nf2 += 8
    store_symp_state!(integrator, cache, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic62ConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = cache
    else
        @unpack a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10 = cache
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)

    # update position
    u = uprev + dt * b1 * kuprev
    # update velocity
    kdu = f.f1(duprev, u, p, integrator.t)
    du = duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b2 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b3 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b4 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b5 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a5 * kdu

    tnew = tnew + a5 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b6 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a6 * kdu

    tnew = tnew + a6 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b7 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a7 * kdu

    tnew = tnew + a7 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b8 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a8 * kdu

    tnew = tnew + a8 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b9 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a9 * kdu

    tnew = tnew + a9 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b10 * ku

    kdu = f.f1(du, u, p, tnew)
    # @.. broadcast=false du = du + dt*a10*kdu
    ku = f.f2(du, u, p, tnew)

    integrator.destats.nf += 10
    integrator.destats.nf2 += 10
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(integrator, cache::Symplectic62Cache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @static if VERSION >= 1.8
        (; a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10) = cache.tab
    else
        @unpack a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10 = cache.tab
    end
    duprev, uprev, _, kuprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # update position
    @.. broadcast=false u=uprev + dt * b1 * kuprev
    # update velocity
    f.f1(kdu, duprev, u, p, integrator.t)
    @.. broadcast=false du=duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b2 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b3 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b4 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b5 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a5 * kdu

    tnew = tnew + a5 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b6 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a6 * kdu

    tnew = tnew + a6 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b7 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a7 * kdu

    tnew = tnew + a7 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b8 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a8 * kdu

    tnew = tnew + a8 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b9 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a9 * kdu

    tnew = tnew + a9 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b10 * ku

    f.f1(kdu, du, u, p, tnew)
    # @.. broadcast=false du = du + dt*a10*kdu
    f.f2(ku, du, u, p, tnew)

    integrator.destats.nf += 10
    integrator.destats.nf2 += 10
    store_symp_state!(integrator, cache, kdu, ku)
end

@muladd function perform_step!(integrator, cache::McAte8ConstantCache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @unpack a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16,
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16 = cache
    duprev, uprev, _, kuprev = load_symp_state(integrator)

    # update position
    u = uprev + dt * b1 * kuprev
    # update velocity
    kdu = f.f1(duprev, u, p, integrator.t)
    du = duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b2 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b3 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b4 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b5 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a5 * kdu

    tnew = tnew + a5 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b6 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a6 * kdu

    tnew = tnew + a6 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b7 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a7 * kdu

    tnew = tnew + a7 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b8 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a8 * kdu

    tnew = tnew + a8 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b9 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a9 * kdu

    tnew = tnew + a9 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b10 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a10 * kdu

    tnew = tnew + a10 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b11 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a11 * kdu

    tnew = tnew + a11 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b12 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a12 * kdu

    tnew = tnew + a12 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b13 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a13 * kdu

    tnew = tnew + a13 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b14 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a14 * kdu

    tnew = tnew + a14 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b15 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a15 * kdu

    tnew = tnew + a15 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b16 * ku

    kdu = f.f1(du, u, p, tnew)
    # @.. broadcast=false du = du + dt*a16*kdu
    ku = f.f2(du, u, p, tnew)

    integrator.destats.nf += 16
    integrator.destats.nf2 += 16
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(integrator, cache::McAte8Cache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @unpack a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16,
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16 = cache.tab
    duprev, uprev, _, kuprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # update position
    @.. broadcast=false u=uprev + dt * b1 * kuprev
    # update velocity
    f.f1(kdu, duprev, u, p, integrator.t)
    @.. broadcast=false du=duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b2 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b3 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b4 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b5 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a5 * kdu

    tnew = tnew + a5 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b6 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a6 * kdu

    tnew = tnew + a6 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b7 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a7 * kdu

    tnew = tnew + a7 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b8 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a8 * kdu

    tnew = tnew + a8 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b9 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a9 * kdu

    tnew = tnew + a9 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b10 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a10 * kdu

    tnew = tnew + a10 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b11 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a11 * kdu

    tnew = tnew + a11 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b12 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a12 * kdu

    tnew = tnew + a12 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b13 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a13 * kdu

    tnew = tnew + a13 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b14 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a14 * kdu

    tnew = tnew + a14 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b15 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a15 * kdu

    tnew = tnew + a15 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b16 * ku

    f.f1(kdu, du, u, p, tnew)
    # @.. broadcast=false du = du + dt*a16*kdu
    f.f2(ku, du, u, p, tnew)

    integrator.destats.nf += 16
    integrator.destats.nf2 += 16
    store_symp_state!(integrator, cache, kdu, ku)
end

@muladd function perform_step!(integrator, cache::KahanLi8ConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @unpack a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18,
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18 = cache
    duprev, uprev, _, kuprev = load_symp_state(integrator)

    # update position
    u = uprev + dt * b1 * kuprev
    # update velocity
    kdu = f.f1(duprev, u, p, integrator.t)
    du = duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b2 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b3 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b4 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b5 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a5 * kdu

    tnew = tnew + a5 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b6 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a6 * kdu

    tnew = tnew + a6 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b7 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a7 * kdu

    tnew = tnew + a7 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b8 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a8 * kdu

    tnew = tnew + a8 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b9 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a9 * kdu

    tnew = tnew + a9 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b10 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a10 * kdu

    tnew = tnew + a10 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b11 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a11 * kdu

    tnew = tnew + a11 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b12 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a12 * kdu

    tnew = tnew + a12 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b13 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a13 * kdu

    tnew = tnew + a13 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b14 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a14 * kdu

    tnew = tnew + a14 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b15 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a15 * kdu

    tnew = tnew + a15 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b16 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a16 * kdu

    tnew = tnew + a16 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b17 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a17 * kdu

    tnew = tnew + a17 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b18 * ku

    kdu = f.f1(du, u, p, tnew)
    # @.. broadcast=false du = du + dt*a18*kdu
    ku = f.f2(du, u, p, tnew)

    integrator.destats.nf += 18
    integrator.destats.nf2 += 18
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(integrator, cache::KahanLi8Cache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @unpack a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18,
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18 = cache.tab
    duprev, uprev, _, kuprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # update position
    @.. broadcast=false u=uprev + dt * b1 * kuprev
    # update velocity
    f.f1(kdu, duprev, u, p, integrator.t)
    @.. broadcast=false du=duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b2 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b3 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b4 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b5 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a5 * kdu

    tnew = tnew + a5 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b6 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a6 * kdu

    tnew = tnew + a6 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b7 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a7 * kdu

    tnew = tnew + a7 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b8 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a8 * kdu

    tnew = tnew + a8 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b9 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a9 * kdu

    tnew = tnew + a9 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b10 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a10 * kdu

    tnew = tnew + a10 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b11 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a11 * kdu

    tnew = tnew + a11 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b12 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a12 * kdu

    tnew = tnew + a12 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b13 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a13 * kdu

    tnew = tnew + a13 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b14 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a14 * kdu

    tnew = tnew + a14 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b15 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a15 * kdu

    tnew = tnew + a15 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b16 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a16 * kdu

    tnew = tnew + a16 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b17 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a17 * kdu

    tnew = tnew + a17 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b18 * ku

    f.f1(kdu, du, u, p, tnew)
    # @.. broadcast=false du = du + dt*a18*kdu
    f.f2(ku, du, u, p, tnew)

    integrator.destats.nf += 18
    integrator.destats.nf2 += 18
    store_symp_state!(integrator, cache, kdu, ku)
end

@muladd function perform_step!(integrator, cache::SofSpa10ConstantCache,
                               repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @unpack a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18,
    a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34,
    a35, a36,
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18,
    b19, b20, b21, b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33, b34,
    b35, b36 = cache
    duprev, uprev, _, kuprev = load_symp_state(integrator)

    # update position
    u = uprev + dt * b1 * kuprev
    # update velocity
    kdu = f.f1(duprev, u, p, integrator.t)
    du = duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b2 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b3 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b4 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b5 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a5 * kdu

    tnew = tnew + a5 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b6 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a6 * kdu

    tnew = tnew + a6 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b7 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a7 * kdu

    tnew = tnew + a7 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b8 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a8 * kdu

    tnew = tnew + a8 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b9 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a9 * kdu

    tnew = tnew + a9 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b10 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a10 * kdu

    tnew = tnew + a10 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b11 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a11 * kdu

    tnew = tnew + a11 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b12 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a12 * kdu

    tnew = tnew + a12 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b13 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a13 * kdu

    tnew = tnew + a13 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b14 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a14 * kdu

    tnew = tnew + a14 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b15 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a15 * kdu

    tnew = tnew + a15 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b16 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a16 * kdu

    tnew = tnew + a16 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b17 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a17 * kdu

    tnew = tnew + a17 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b18 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a18 * kdu

    tnew = tnew + a18 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b19 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a19 * kdu

    tnew = tnew + a19 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b20 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a20 * kdu

    tnew = tnew + a20 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b21 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a21 * kdu

    tnew = tnew + a21 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b22 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a22 * kdu

    tnew = tnew + a22 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b23 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a23 * kdu

    tnew = tnew + a23 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b24 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a24 * kdu

    tnew = tnew + a24 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b25 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a25 * kdu

    tnew = tnew + a25 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b26 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a26 * kdu

    tnew = tnew + a26 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b27 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a27 * kdu

    tnew = tnew + a27 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b28 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a28 * kdu

    tnew = tnew + a28 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b29 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a29 * kdu

    tnew = tnew + a29 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b30 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a30 * kdu

    tnew = tnew + a30 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b31 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a31 * kdu

    tnew = tnew + a31 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b32 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a32 * kdu

    tnew = tnew + a32 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b33 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a33 * kdu

    tnew = tnew + a33 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b34 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a34 * kdu

    tnew = tnew + a34 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b35 * ku

    kdu = f.f1(du, u, p, tnew)
    du = du + dt * a35 * kdu

    tnew = tnew + a35 * dt
    ku = f.f2(du, u, p, tnew)
    u = u + dt * b36 * ku

    kdu = f.f1(du, u, p, tnew)
    # @.. broadcast=false du = du + dt*a30*kdu
    ku = f.f2(du, u, p, tnew)

    integrator.destats.nf += 36
    integrator.destats.nf2 += 36
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(integrator, cache::SofSpa10Cache, repeat_step = false)
    @static if VERSION >= 1.8
        (; t, dt, f, p) = integrator
    else
        @unpack t, dt, f, p = integrator
    end
    @unpack a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18,
    a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34,
    a35, a36,
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18,
    b19, b20, b21, b22, b23, b24, b25, b26, b27, b28, b29, b30, b31, b32, b33, b34,
    b35, b36 = cache.tab
    duprev, uprev, _, kuprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # update position
    @.. broadcast=false u=uprev + dt * b1 * kuprev
    # update velocity
    f.f1(kdu, duprev, u, p, integrator.t)
    @.. broadcast=false du=duprev + dt * a1 * kdu
    # update position & velocity
    tnew = t + a1 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b2 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a2 * kdu

    # update position & velocity
    tnew = tnew + a2 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b3 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a3 * kdu

    # update position & velocity
    tnew = tnew + a3 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b4 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a4 * kdu

    # update position & velocity
    tnew = tnew + a4 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b5 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a5 * kdu

    tnew = tnew + a5 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b6 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a6 * kdu

    tnew = tnew + a6 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b7 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a7 * kdu

    tnew = tnew + a7 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b8 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a8 * kdu

    tnew = tnew + a8 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b9 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a9 * kdu

    tnew = tnew + a9 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b10 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a10 * kdu

    tnew = tnew + a10 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b11 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a11 * kdu

    tnew = tnew + a11 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b12 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a12 * kdu

    tnew = tnew + a12 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b13 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a13 * kdu

    tnew = tnew + a13 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b14 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a14 * kdu

    tnew = tnew + a14 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b15 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a15 * kdu

    tnew = tnew + a15 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b16 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a16 * kdu

    tnew = tnew + a16 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b17 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a17 * kdu

    tnew = tnew + a17 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b18 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a18 * kdu

    tnew = tnew + a18 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b19 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a19 * kdu

    tnew = tnew + a19 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b20 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a20 * kdu

    tnew = tnew + a20 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b21 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a21 * kdu

    tnew = tnew + a21 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b22 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a22 * kdu

    tnew = tnew + a22 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b23 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a23 * kdu

    tnew = tnew + a23 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b24 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a24 * kdu

    tnew = tnew + a24 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b25 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a25 * kdu

    tnew = tnew + a25 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b26 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a26 * kdu

    tnew = tnew + a26 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b27 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a27 * kdu

    tnew = tnew + a27 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b28 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a28 * kdu

    tnew = tnew + a28 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b29 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a29 * kdu

    tnew = tnew + a29 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b30 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a30 * kdu

    tnew = tnew + a30 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b31 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a31 * kdu

    tnew = tnew + a31 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b32 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a32 * kdu

    tnew = tnew + a32 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b33 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a33 * kdu

    tnew = tnew + a33 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b34 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a34 * kdu

    tnew = tnew + a34 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b35 * ku

    f.f1(kdu, du, u, p, tnew)
    @.. broadcast=false du=du + dt * a35 * kdu

    tnew = tnew + a35 * dt
    f.f2(ku, du, u, p, tnew)
    @.. broadcast=false u=u + dt * b36 * ku

    f.f1(kdu, du, u, p, tnew)
    # @.. broadcast=false du = du + dt*a30*kdu
    f.f2(ku, du, u, p, tnew)

    integrator.destats.nf += 36
    integrator.destats.nf2 += 36
    store_symp_state!(integrator, cache, kdu, ku)
end
