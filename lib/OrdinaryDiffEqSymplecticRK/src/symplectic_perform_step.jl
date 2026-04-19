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
    integrator.stats.nf2 += 1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    integrator.fsalfirst = ArrayPartition((kdu, kuprev))
    return integrator.fsallast = ArrayPartition((zero(kdu), ku))
end

@muladd function perform_step!(
        integrator, cache::SymplecticEulerConstantCache,
        repeat_step = false
    )
    (; t, dt, f, p) = integrator
    duprev, uprev = integrator.uprev.x
    kuprev = integrator.fsalfirst.x[2]
    u = uprev + dt * kuprev
    # Now actually compute the step
    # Do it at the end for interpolations!
    kdu = f.f1(duprev, u, p, t)
    du = duprev + dt * kdu

    ku = f.f2(du, u, p, t)
    integrator.stats.nf2 += 1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((kdu, ku))
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
end

function initialize!(integrator, cache::SymplecticEulerCache)
    integrator.kshortsize = 2
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
    @muladd @.. broadcast = false du = duprev + integrator.dt * kdu
    integrator.f.f2(ku, du, uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return integrator.stats.nf2 += 2
end

@muladd function perform_step!(integrator, cache::SymplecticEulerCache, repeat_step = false)
    (; t, dt, f, p) = integrator
    duprev, uprev = integrator.uprev.x
    du, u = integrator.u.x
    kuprev = integrator.fsalfirst.x[2]
    kdu, ku = integrator.fsallast.x
    @.. broadcast = false u = uprev + dt * kuprev
    # Now actually compute the step
    # Do it at the end for interpolations!
    integrator.stats.nf2 += 1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    f.f1(kdu, duprev, u, p, t)
    @.. broadcast = false du = duprev + dt * kdu
    f.f2(ku, du, u, p, t)
end

# some of the algorithms are designed only for the case
# f.f2(p, q, pa, t) = p which is the Newton/Lagrange equations
# If called with different functions (which are possible in the Hamiltonian case)
# an exception is thrown to avoid silently calculate wrong results.
function verify_f2(
        f::F, p, q, pa, t, ::Any,
        ::C
    ) where {
        F, C <: Union{
            HamiltonConstantCache, VerletLeapfrogConstantCache,
            LeapfrogDriftKickDriftConstantCache,
        },
    }
    return f(p, q, pa, t)
end
function verify_f2(
        f::F, res, p, q, pa, t, ::Any,
        ::C
    ) where {
        F, C <: Union{
            HamiltonMutableCache, VerletLeapfrogCache,
            LeapfrogDriftKickDriftCache,
        },
    }
    return f(res, p, q, pa, t)
end

function verify_f2(f::F, p, q, pa, t, integrator, ::C) where {F, C <: VelocityVerletConstantCache}
    res = f(p, q, pa, t)
    return res == p ? p : throwex(integrator)
end
function verify_f2(f::F, res, p, q, pa, t, integrator, ::C) where {F, C <: VelocityVerletCache}
    f(res, p, q, pa, t)
    return res == p ? res : throwex(integrator)
end
function throwex(integrator)
    algn = typeof(integrator.alg)
    throw(ArgumentError("Algorithm $algn is not applicable if f2(p, q, t) != p"))
end

# provide the mutable uninitialized objects to keep state and derivative in case of mutable caches
# no such objects are required for constant caches
function alloc_symp_state(integrator)
    return (integrator.u.x..., integrator.cache.tmp.x...)
end

# load state and derivatives at begin of symplectic iteration steps
function load_symp_state(integrator)
    return (integrator.uprev.x..., integrator.fsallast.x...)
end

# store state and derivatives at the end of symplectic iteration steps
function store_symp_state!(integrator, ::OrdinaryDiffEqConstantCache, du, u, kdu, ku)
    integrator.u = ArrayPartition((du, u))
    integrator.fsallast = ArrayPartition((kdu, ku))
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return nothing
end

function store_symp_state!(integrator, ::OrdinaryDiffEqMutableCache, kdu, ku)
    copyto!(integrator.k[1].x[1], integrator.k[2].x[1])
    copyto!(integrator.k[1].x[2], integrator.k[2].x[2])
    copyto!(integrator.k[2].x[2], ku)
    copyto!(integrator.k[2].x[1], kdu)
    return nothing
end

function initialize!(
        integrator,
        cache::C
    ) where {
        C <: Union{
            HamiltonMutableCache, VelocityVerletCache,
            VerletLeapfrogCache, LeapfrogDriftKickDriftCache,
        },
    }
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    duprev, uprev = integrator.uprev.x
    integrator.f.f1(integrator.k[2].x[1], duprev, uprev, integrator.p, integrator.t)
    verify_f2(
        integrator.f.f2, integrator.k[2].x[2], duprev, uprev, integrator.p,
        integrator.t, integrator, cache
    )
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return integrator.stats.nf2 += 1
end

function initialize!(
        integrator,
        cache::C
    ) where {
        C <: Union{
            HamiltonConstantCache, VelocityVerletConstantCache,
            VerletLeapfrogConstantCache, LeapfrogDriftKickDriftConstantCache,
        },
    }
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    duprev, uprev = integrator.uprev.x
    kdu = integrator.f.f1(duprev, uprev, integrator.p, integrator.t)
    ku = verify_f2(
        integrator.f.f2, duprev, uprev, integrator.p, integrator.t, integrator,
        cache
    )
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.fsallast = ArrayPartition((kdu, ku))
    integrator.k[2] = integrator.fsallast
    return integrator.fsalfirst = integrator.fsallast
end

@muladd function perform_step!(
        integrator, cache::VelocityVerletConstantCache,
        repeat_step = false
    )
    (; t, dt, f, p) = integrator
    duprev, uprev = load_symp_state(integrator)

    # x(t+Δt) = x(t) + v(t)*Δt + 1/2*a(t)*Δt^2
    ku = integrator.fsallast.x[1]
    dtsq = dt^2
    half = cache.half
    u = uprev + dt * duprev + dtsq * (half * ku)
    kdu = f.f1(duprev, u, p, t + dt)
    # v(t+Δt) = v(t) + 1/2*(a(t)+a(t+Δt))*Δt
    du = duprev + dt * (half * ku + half * kdu)

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    store_symp_state!(integrator, cache, du, u, kdu, du)
end

@muladd function perform_step!(integrator, cache::VelocityVerletCache, repeat_step = false)
    (; t, dt, f, p) = integrator
    duprev, uprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # x(t+Δt) = x(t) + v(t)*Δt + 1/2*a(t)*Δt^2
    ku = integrator.fsallast.x[1]
    dtsq = dt^2
    half = cache.half
    @.. broadcast = false u = uprev + dt * duprev + dtsq * (half * ku)
    f.f1(kdu, duprev, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # v(t+Δt) = v(t) + 1/2*(a(t)+a(t+Δt))*Δt
    @.. broadcast = false du = duprev + dt * (half * ku + half * kdu)

    store_symp_state!(integrator, cache, kdu, du)
end

@muladd function perform_step!(
        integrator, cache::VerletLeapfrogConstantCache,
        repeat_step = false
    )
    (; t, dt, f, p) = integrator
    duprev, uprev, kduprev, _ = load_symp_state(integrator)

    # kick-drift-kick scheme of the Leapfrog method:
    # update velocity
    half = cache.half
    du = duprev + dt * half * kduprev

    # update position
    ku = f.f2(du, uprev, p, t + half * dt)
    u = uprev + dt * ku

    # update velocity
    kdu = f.f1(du, u, p, t + dt)
    du = du + dt * half * kdu

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(integrator, cache::VerletLeapfrogCache, repeat_step = false)
    (; t, dt, f, p) = integrator
    duprev, uprev, kduprev, _ = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # kick-drift-kick scheme of the Leapfrog method:
    # update velocity
    half = cache.half
    @.. broadcast = false du = duprev + dt * half * kduprev

    # update position
    f.f2(ku, du, uprev, p, t + half * dt)
    @.. broadcast = false u = uprev + dt * ku

    # update velocity
    f.f1(kdu, du, u, p, t + dt)
    @.. broadcast = false du = du + dt * half * kdu

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    store_symp_state!(integrator, cache, kdu, ku)
end

@muladd function perform_step!(
        integrator, cache::LeapfrogDriftKickDriftConstantCache,
        repeat_step = false
    )
    (; t, dt, f, p) = integrator
    duprev, uprev, _, _ = load_symp_state(integrator)

    # drift-kick-drift scheme of the Leapfrog method, allowing for f1 to depend on v:
    # update position half step
    half = cache.half
    ku = f.f2(duprev, uprev, p, t)
    u = uprev + dt * half * ku

    # update velocity half step
    kdu = f.f1(duprev, uprev, p, t)
    du = duprev + dt * half * kdu

    # update velocity (add to previous full step velocity)
    # note that this extra step is only necessary if f1 depends on v/du (or t)
    kdu = f.f1(du, u, p, t + half * dt)
    du = duprev + dt * kdu

    # update position (add to half step position)
    ku = f.f2(du, u, p, t + dt)
    u = u + dt * half * ku

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    integrator.stats.nf2 += 2
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

@muladd function perform_step!(
        integrator, cache::LeapfrogDriftKickDriftCache,
        repeat_step = false
    )
    (; t, dt, f, p) = integrator
    duprev, uprev, _, _ = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # drift-kick-drift scheme of the Leapfrog method, allowing for f1 to depend on v:
    # update position half step
    half = cache.half
    f.f2(ku, duprev, uprev, p, t)
    @.. broadcast = false u = uprev + dt * half * ku

    # update velocity half step
    f.f1(kdu, duprev, uprev, p, t)
    @.. broadcast = false du = duprev + dt * half * kdu

    # update velocity (add to previous full step velocity)
    # note that this extra step is only necessary if f1 depends on v/du (or t)
    f.f1(kdu, du, u, p, t + half * dt)
    @.. broadcast = false du = duprev + dt * kdu

    # update position (add to half step position)
    f.f2(ku, du, u, p, t + dt)
    @.. broadcast = false u = u + dt * half * ku

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
    integrator.stats.nf2 += 2
    store_symp_state!(integrator, cache, kdu, ku)
end

# ===== Generic symplectic drift-kick perform_step! =====

# Out-of-place (OOP) generic symplectic integrator
@muladd function perform_step!(
        integrator, cache::SymplecticTableau,
        repeat_step = false
    )
    (; t, dt, f, p) = integrator
    a = cache.a
    b = cache.b
    N = length(a)
    duprev, uprev, _, kuprev = load_symp_state(integrator)

    # Stage 1: drift using kuprev from FSAL
    u = uprev + dt * b[1] * kuprev
    kdu = f.f1(duprev, u, p, t)
    du = duprev + dt * a[1] * kdu
    tnew = t + a[1] * dt

    # Stages 2 through N
    for i in 2:N
        ku = f.f2(du, u, p, tnew)
        u = u + dt * b[i] * ku
        kdu = f.f1(du, u, p, tnew)
        du = du + dt * a[i] * kdu
        if i < N
            tnew = tnew + a[i] * dt
        end
    end

    # FSAL finalization: if a[N] != 0, du changed so we need fresh kdu
    if !iszero(a[N])
        kdu = f.f1(du, u, p, tnew)
    end
    ku = f.f2(du, u, p, tnew)

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, N + (!iszero(a[N]) ? 1 : 0))
    integrator.stats.nf2 += N
    store_symp_state!(integrator, cache, du, u, kdu, ku)
end

# In-place (IIP) generic symplectic integrator
@muladd function perform_step!(integrator, cache::SymplecticGenericCache, repeat_step = false)
    (; t, dt, f, p) = integrator
    tab = cache.tab
    a = tab.a
    b = tab.b
    N = length(a)
    duprev, uprev, _, kuprev = load_symp_state(integrator)
    du, u, kdu, ku = alloc_symp_state(integrator)

    # Stage 1: drift using kuprev from FSAL
    @.. broadcast = false u = uprev + dt * b[1] * kuprev
    f.f1(kdu, duprev, u, p, t)
    @.. broadcast = false du = duprev + dt * a[1] * kdu
    tnew = t + a[1] * dt

    # Stages 2 through N
    for i in 2:N
        f.f2(ku, du, u, p, tnew)
        @.. broadcast = false u = u + dt * b[i] * ku
        f.f1(kdu, du, u, p, tnew)
        @.. broadcast = false du = du + dt * a[i] * kdu
        if i < N
            tnew = tnew + a[i] * dt
        end
    end

    # FSAL finalization: if a[N] != 0, du changed so we need fresh kdu
    if !iszero(a[N])
        f.f1(kdu, du, u, p, tnew)
    end
    f.f2(ku, du, u, p, tnew)

    OrdinaryDiffEqCore.increment_nf!(integrator.stats, N + (!iszero(a[N]) ? 1 : 0))
    integrator.stats.nf2 += N
    store_symp_state!(integrator, cache, kdu, ku)
end
