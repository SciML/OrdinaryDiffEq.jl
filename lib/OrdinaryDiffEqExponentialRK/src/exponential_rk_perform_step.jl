# Helper function to compute the G_nj factors for the classical ExpRK methods
@inline _compute_nl(f::SplitFunction, u, p, t, A) = f.f2(u, p, t)
@inline _compute_nl(f, u, p, t, A) = f(u, p, t) - A * u
@inline _compute_nl!(G, f::SplitFunction, u, p, t, A, Au_cache) = f.f2(G, u, p, t)
@inline function _compute_nl!(G, f, u, p, t, A, Au_cache)
    f(G, u, p, t)
    mul!(Au_cache, A, u)
    return G .-= Au_cache
end

# Error estimate of an embedded ExpRK pair. Called where `utilde` is formed, so
# that it is never read on a path that does not define it.
@inline function _expRK_EEst!(integrator, utilde, u)
    atmp = calculate_residuals(
        utilde, integrator.uprev, u, integrator.opts.abstol,
        integrator.opts.reltol, integrator.opts.internalnorm, integrator.t
    )
    return OrdinaryDiffEqCore.set_EEst!(
        integrator, integrator.opts.internalnorm(atmp, integrator.t)
    )
end

##########################################
# Common initializers for ExpRK integrators
function initialize!(integrator, cache::ExpRKConstantCache)
    # Pre-start fsal
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)

    # Initialize interpolation derivatives
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end
function initialize!(integrator, cache::ExpRKCache)
    # Pre-start fsal
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Initialize interpolation derivatives
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

###########################################
# Classical ExpRK integrators
function perform_step!(integrator, cache::LawsonEulerConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = isa(f, SplitFunction) ? f.f1.f : calc_J(integrator, cache) # get linear operator
    alg = unwrap_alg(integrator, true)

    nl = _compute_nl(f, uprev, p, t, A)
    if isa(f, SplitFunction)
        integrator.stats.nf2 += 1
    else
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end
    @muladd v = uprev + dt * nl
    if alg.krylov
        u = expv(
            dt, A, v; m = min(alg.m, size(A, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop
        )
    else
        exphA = cache.ops
        u = exphA * v
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::LawsonEulerCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, G, J, exphA, KsCache) = cache
    A = isa(f, SplitFunction) ? f.f1.f : (calc_J!(J, integrator, cache); J) # get linear operator
    alg = unwrap_alg(integrator, true)

    _compute_nl!(G, f, uprev, p, t, A, rtmp)
    if isa(f, SplitFunction)
        integrator.stats.nf2 += 1
    else
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end
    @muladd @.. broadcast = false tmp = uprev + dt * G
    if alg.krylov
        Ks, expv_cache = KsCache
        arnoldi!(
            Ks, A, tmp; m = min(alg.m, size(A, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop
        )
        expv!(u, dt, Ks; cache = expv_cache)
    else
        mul!(u, exphA, tmp)
    end

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::NorsettEulerConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = isa(f, SplitFunction) ? f.f1.f : calc_J(integrator, cache) # get linear operator
    alg = unwrap_alg(integrator, true)

    if alg.krylov
        w = phiv(
            dt, A, integrator.fsalfirst, 1; m = min(alg.m, size(A, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop
        )
        u = uprev + dt * w[:, 2]
    else
        phihA = cache.ops
        u = uprev + dt * (phihA * integrator.fsalfirst)
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::NorsettEulerCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; rtmp, J, KsCache) = cache
    A = isa(f, SplitFunction) ? f.f1.f : (calc_J!(J, integrator, cache); J) # get linear operator
    alg = unwrap_alg(integrator, true)

    if alg.krylov
        Ks, phiv_cache, ws, herm = KsCache
        w = ws[1]
        arnoldi!(Ks, A, integrator.fsalfirst; _arnoldi_kwargs(alg, A, integrator, herm)...)
        phiv!(w, dt, Ks, 1; cache = phiv_cache)
        @muladd @.. broadcast = false u = uprev + dt * @view(w[:, 2])
    else
        mul!(rtmp, cache.phihA, integrator.fsalfirst)
        @muladd @.. broadcast = false u = uprev + dt * rtmp
    end

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::ETDRK2ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = isa(f, SplitFunction) ? f.f1.f : calc_J(integrator, cache) # get linear operator
    alg = unwrap_alg(integrator, true)

    if alg.krylov
        F1 = integrator.fsalfirst
        w1 = phiv(
            dt, A, F1, 2; m = min(alg.m, size(A, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop
        )
        U2 = uprev + dt * w1[:, 2]
        F2 = _compute_nl(f, U2, p, t + dt, A) + A * uprev
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 1
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
        w2 = phiv(
            dt, A, F2, 2; m = min(alg.m, size(A, 1)),
            opnorm = integrator.opts.internalopnorm, iop = alg.iop
        )
        u = uprev + dt * (w1[:, 2] - w1[:, 3] + w2[:, 3])
    else
        phi1, phi2 = cache.ops
        # The caching version uses a special formula to save computation
        G1 = f.f2(uprev, p, t)
        integrator.stats.nf2 += 1
        F1 = integrator.fsalfirst
        U2 = uprev + dt * (phi1 * F1)
        G2 = f.f2(U2, p, t + dt)
        integrator.stats.nf2 += 1
        u = U2 + dt * (phi2 * (G2 - G1))
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::ETDRK2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, F2, J, KsCache) = cache
    A = isa(f, SplitFunction) ? f.f1.f : (calc_J!(J, integrator, cache); J) # get linear operator
    alg = unwrap_alg(integrator, true)

    if alg.krylov
        F1 = integrator.fsalfirst
        Ks, phiv_cache, ws, herm = KsCache
        w1, w2 = ws
        kwargs = _arnoldi_kwargs(alg, A, integrator, herm)
        # Krylov for F1
        arnoldi!(Ks, A, F1; kwargs...)
        phiv!(w1, dt, Ks, 2; cache = phiv_cache)
        # Krylov for F2
        @muladd @.. broadcast = false tmp = uprev + dt * @view(w1[:, 2])
        _compute_nl!(F2, f, tmp, p, t + dt, A, rtmp)
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 1
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
        F2 .+= mul!(rtmp, A, uprev)
        arnoldi!(Ks, A, F2; kwargs...)
        phiv!(w2, dt, Ks, 2; cache = phiv_cache)
        # Update u
        u .= uprev
        axpy!(dt, @view(w1[:, 2]), u)
        axpy!(-dt, @view(w1[:, 3]), u)
        axpy!(dt, @view(w2[:, 3]), u)
    else
        phi1, phi2 = cache.ops
        F1 = integrator.fsalfirst
        # The caching version uses a special formula to save computation
        # Compute U2
        mul!(rtmp, phi1, F1)
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U2
        # Compute G2 - G1, storing result in the cache F2
        f.f2(rtmp, uprev, p, t)
        integrator.stats.nf2 += 1
        f.f2(F2, tmp, p, t + dt)
        integrator.stats.nf2 += 1
        F2 .-= rtmp # "F2" is G2 - G1
        # Update u
        u .= tmp
        axpy!(dt, mul!(rtmp, phi2, F2), u)
    end

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::ETDRK3ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = isa(f, SplitFunction) ? f.f1.f : calc_J(integrator, cache) # get linear operator
    alg = unwrap_alg(integrator, true)

    Au = A * uprev
    F1 = integrator.fsalfirst
    if alg.krylov
        kwargs = _arnoldi_kwargs(alg, A, integrator, nothing)
        # Krylov on F1 (first column)
        Ks = arnoldi(A, F1; kwargs...)
        w1_half = phiv(dt / 2, Ks, 1)
        w1 = phiv(dt, Ks, 3)
        U2 = uprev + dt / 2 * w1_half[:, 2]
        F2 = _compute_nl(f, U2, p, t + dt / 2, A) + Au
        # Krylov on F2 (second column)
        w2 = phiv(dt, A, F2, 3; kwargs...)
        U3 = uprev + dt * (2w2[:, 2] - w1[:, 2])
        F3 = _compute_nl(f, U3, p, t + dt, A) + Au
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 2
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
        end
        # Krylov on F3 (third column)
        w3 = phiv(dt, A, F3, 3; kwargs...)
        u = uprev +
            dt * (
            4w1[:, 4] - 3w1[:, 3] + w1[:, 2]
                -
                8w2[:, 4] + 4w2[:, 3]
                + 4w3[:, 4] - w3[:, 3]
        )
    else
        A21, A3, B1, B2, B3 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        U2 = uprev + dt * (A21 * F1)
        F2 = f.f2(U2, p, t + dt / 2) + Au
        integrator.stats.nf2 += 1
        # stage 3
        U3 = uprev + dt * (A3 * (2F2 - F1))
        F3 = f.f2(U3, p, t + dt) + Au
        integrator.stats.nf2 += 1
        # update u
        u = uprev + dt * (B1 * F1 + B2 * F2 + B3 * F3)
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::ETDRK3Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, Au, F2, F3, J, KsCache) = cache
    A = isa(f, SplitFunction) ? f.f1.f : (calc_J!(J, integrator, cache); J) # get linear operator
    alg = unwrap_alg(integrator, true)

    F1 = integrator.fsalfirst
    mul!(Au, A, uprev)
    halfdt = dt / 2
    if alg.krylov
        Ks, phiv_cache, ws, herm = KsCache
        w1_half, w1, w2, w3 = ws
        kwargs = _arnoldi_kwargs(alg, A, integrator, herm)
        # Krylov for F1 (first column)
        arnoldi!(Ks, A, F1; kwargs...)
        phiv!(w1_half, halfdt, Ks, 1; cache = phiv_cache)
        phiv!(w1, dt, Ks, 3; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev + halfdt * w1_half[:, 2] # tmp is U2
        _compute_nl!(F2, f, tmp, p, t + halfdt, A, rtmp)
        F2 .+= Au
        # Krylov for F2 (second column)
        arnoldi!(Ks, A, F2; kwargs...)
        phiv!(w2, dt, Ks, 3; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev + dt * (2 * w2[:, 2] - w1[:, 2]) # tmp is U3
        _compute_nl!(F3, f, tmp, p, t + dt, A, rtmp)
        F3 .+= Au
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 2
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)
        end
        # Krylov for F3 (third column)
        arnoldi!(Ks, A, F3; kwargs...)
        phiv!(w3, dt, Ks, 3; cache = phiv_cache)
        # Update u
        @views @.. broadcast = false rtmp = 4w1[:, 4] - 3w1[:, 3] + w1[:, 2] - 8w2[:, 4] +
            4w2[:, 3] + 4w3[:, 4] - w3[:, 3]
        @muladd @.. broadcast = false u = uprev + dt * rtmp
    else
        A21, A3, B1, B2, B3 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        mul!(rtmp, A21, F1)
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U2
        f.f2(F2, tmp, p, t + halfdt)
        F2 .+= Au
        integrator.stats.nf2 += 1
        # stage 3
        @muladd @.. broadcast = false F3 = 2 * F2 - F1 # use F3 temporarily as cache
        mul!(rtmp, A3, F3)
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U3
        f.f2(F3, tmp, p, t + dt)
        F3 .+= Au
        integrator.stats.nf2 += 1
        # update u
        u .= uprev
        axpy!(dt, mul!(rtmp, B1, F1), u)
        axpy!(dt, mul!(rtmp, B2, F2), u)
        axpy!(dt, mul!(rtmp, B3, F3), u)
    end

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::ETDRK4ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = isa(f, SplitFunction) ? f.f1.f : calc_J(integrator, cache) # get linear operator
    alg = unwrap_alg(integrator, true)

    Au = A * uprev
    F1 = integrator.fsalfirst
    halfdt = dt / 2
    if alg.krylov
        kwargs = _arnoldi_kwargs(alg, A, integrator, nothing)
        # Krylov on F1 (first column)
        Ks = arnoldi(A, F1; kwargs...)
        w1_half = phiv(halfdt, Ks, 1)
        w1 = phiv(dt, Ks, 3)
        U2 = uprev + halfdt * w1_half[:, 2]
        F2 = _compute_nl(f, U2, p, t + halfdt, A) + Au
        # Krylov on F2 (second column)
        Ks = arnoldi(A, F2; kwargs...)
        w2_half = phiv(halfdt, Ks, 1)
        w2 = phiv(dt, Ks, 3)
        U3 = uprev + halfdt * w2_half[:, 2]
        F3 = _compute_nl(f, U3, p, t + halfdt, A) + Au
        # Krylov on F3 (third column)
        w3 = phiv(dt, A, F3, 3; kwargs...)
        # Extra Krylov for computing F4
        rtmp = 2F3 - F1 - Au + A * U2
        wtmp = phiv(halfdt, A, rtmp, 1; kwargs...)
        U4 = U2 + halfdt * wtmp[:, 2]
        F4 = _compute_nl(f, U4, p, t + dt, A) + Au
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 3
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
        end
        # Krylov on F4 (fourth column)
        w4 = phiv(dt, A, F4, 3; kwargs...)
        # update u
        u = uprev +
            dt * (
            w1[:, 2] - 3w1[:, 3] + 4w1[:, 4] + 2w2[:, 3] - 4w2[:, 4] +
                2w3[:, 3] - 4w3[:, 4] + 4w4[:, 4] - w4[:, 3]
        )
    else
        A21, A41, A43, B1, B2, B4 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        U2 = uprev + dt * (A21 * F1)
        F2 = f.f2(U2, p, t + halfdt) + Au
        # stage 3
        U3 = uprev + dt * (A21 * F2) # A32 = A21
        F3 = f.f2(U3, p, t + halfdt) + Au
        # stage 4
        U4 = uprev + dt * (A41 * F1 + A43 * F3)
        F4 = f.f2(U4, p, t + dt) + Au
        integrator.stats.nf2 += 3
        # update u
        u = uprev + dt * (B1 * F1 + B2 * (F2 + F3) + B4 * F4) # B3 = B2
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::ETDRK4Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, Au, F2, F3, F4, J, KsCache) = cache
    A = isa(f, SplitFunction) ? f.f1.f : (calc_J!(J, integrator, cache); J) # get linear operator
    alg = unwrap_alg(integrator, true)

    F1 = integrator.fsalfirst
    mul!(Au, A, uprev)
    halfdt = dt / 2
    if alg.krylov
        Ks, phiv_cache, ws, herm = KsCache
        w1_half, w2_half, w1, w2, w3, w4 = ws
        kwargs = _arnoldi_kwargs(alg, A, integrator, herm)
        # Krylov for F1 (first column)
        arnoldi!(Ks, A, F1; kwargs...)
        phiv!(w1_half, halfdt, Ks, 1; cache = phiv_cache)
        phiv!(w1, dt, Ks, 3; cache = phiv_cache)
        U2 = u # temporarily use u to store U2 (used in the extra Krylov step)
        @muladd @.. broadcast = false @views U2 = uprev + halfdt * w1_half[:, 2]
        _compute_nl!(F2, f, U2, p, t + halfdt, A, rtmp)
        F2 .+= Au
        # Krylov for F2 (second column)
        arnoldi!(Ks, A, F2; kwargs...)
        phiv!(w2_half, halfdt, Ks, 1; cache = phiv_cache)
        phiv!(w2, dt, Ks, 3; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev + halfdt * w2_half[:, 2] # tmp is U3
        _compute_nl!(F3, f, tmp, p, t + halfdt, A, rtmp)
        F3 .+= Au
        # Krylov for F3 (third column)
        arnoldi!(Ks, A, F3; kwargs...)
        phiv!(w3, dt, Ks, 3; cache = phiv_cache)
        # Extra Krylov for computing F4
        # Compute rtmp = 2F3 - F1 - Au + A*U2
        mul!(rtmp, A, U2)
        @.. broadcast = false rtmp += 2F3 - F1 - Au
        arnoldi!(Ks, A, rtmp; kwargs...)
        phiv!(w1_half, halfdt, Ks, 1; cache = phiv_cache) # original w1_half is no longer needed
        @muladd @.. broadcast = false @views tmp = U2 + halfdt * w1_half[:, 2] # tmp is U4
        _compute_nl!(F4, f, tmp, p, t + dt, A, rtmp)
        F4 .+= Au
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 3
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
        end
        # Krylov for F4 (fourth column)
        arnoldi!(Ks, A, F4; kwargs...)
        phiv!(w4, dt, Ks, 3; cache = phiv_cache)
        # update u
        @views @.. broadcast = false rtmp = w1[:, 2] - 3w1[:, 3] + 4w1[:, 4] + 2w2[
            :, 3,
        ] -
            4w2[:, 4] +
            2w3[:, 3] - 4w3[:, 4] + 4w4[:, 4] - w4[:, 3]
        @muladd @.. broadcast = false u = uprev + dt * rtmp
    else
        A21, A41, A43, B1, B2, B4 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        mul!(rtmp, A21, F1)
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U2
        f.f2(F2, tmp, p, t + halfdt)
        F2 .+= Au
        # stage 3
        mul!(rtmp, A21, F2) # A32 = A21
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U3
        f.f2(F3, tmp, p, t + halfdt)
        F3 .+= Au
        # stage 4
        @.. broadcast = false tmp = uprev
        axpy!(dt, mul!(rtmp, A41, F1), tmp)
        axpy!(dt, mul!(rtmp, A43, F3), tmp) # tmp is U4
        f.f2(F4, tmp, p, t + dt)
        F4 .+= Au
        integrator.stats.nf2 += 3
        # update u
        u .= uprev
        axpy!(dt, mul!(rtmp, B1, F1), u)
        F2 .+= F3
        axpy!(dt, mul!(rtmp, B2, F2), u) # B3 = B2
        axpy!(dt, mul!(rtmp, B4, F4), u)
    end

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::HochOst4ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = isa(f, SplitFunction) ? f.f1.f : calc_J(integrator, cache) # get linear operator
    alg = unwrap_alg(integrator, true)

    Au = A * uprev
    F1 = integrator.fsalfirst
    halfdt = dt / 2
    if alg.krylov
        kwargs = _arnoldi_kwargs(alg, A, integrator, nothing)
        # Krylov on F1 (first column)
        Ks = arnoldi(A, F1; kwargs...)
        w1_half = phiv(halfdt, Ks, 3)
        w1 = phiv(dt, Ks, 3)
        U2 = uprev + halfdt * w1_half[:, 2]
        F2 = _compute_nl(f, U2, p, t + halfdt, A) + Au
        # Krylov on F2 (second column)
        Ks = arnoldi(A, F2; kwargs...)
        w2_half = phiv(halfdt, Ks, 3)
        w2 = phiv(dt, Ks, 3)
        U3 = uprev + dt * (0.5w1_half[:, 2] - w1_half[:, 3] + w2_half[:, 3])
        F3 = _compute_nl(f, U3, p, t + halfdt, A) + Au
        # Krylov on F3 (third column)
        Ks = arnoldi(A, F3; kwargs...)
        w3_half = phiv(halfdt, Ks, 3)
        w3 = phiv(dt, Ks, 3)
        U4 = uprev + dt * (w1[:, 2] - 2w1[:, 3] + w2[:, 3] + w3[:, 3])
        F4 = _compute_nl(f, U4, p, t + dt, A) + Au
        # Krylov on F4 (fourth column)
        Ks = arnoldi(A, F4; kwargs...)
        w4_half = phiv(halfdt, Ks, 3)
        w4 = phiv(dt, Ks, 3)
        U5 = uprev +
            dt * (
            0.5w1_half[:, 2] - 0.75w1_half[:, 3] + 0.5w1_half[:, 4] + w1[:, 4] -
                0.25w1[:, 3] +
                0.5w2_half[:, 3] - w2[:, 4] + 0.25w2[:, 3] - 0.5w2_half[:, 4] +
                0.5w3_half[:, 3] - w2[:, 4] + 0.25w3[:, 3] - 0.5w3_half[:, 4] +
                w4[:, 4] - 0.25w4[:, 3] - 0.25w4_half[:, 3] + 0.5w4_half[:, 4]
        )
        F5 = _compute_nl(f, U5, p, t + halfdt, A) + Au
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 4
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
        end
        # Krylov on F5 (fifth column)
        w5 = phiv(dt, A, F5, 3; kwargs...)
        # update u
        u = uprev +
            dt * (
            w1[:, 2] - 3w1[:, 3] + 4w1[:, 4] - w4[:, 3] + 4w4[:, 4] + 4w5[:, 3] -
                8w5[:, 4]
        )
    else
        A21, A31, A32, A41, A42, A51, A52, A54, B1, B4, B5 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        U2 = uprev + dt * (A21 * F1)
        F2 = f.f2(U2, p, t + halfdt) + Au
        # stage 3
        U3 = uprev + dt * (A31 * F1 + A32 * F2)
        F3 = f.f2(U3, p, t + halfdt) + Au
        # stage 4
        U4 = uprev + dt * (A41 * F1 + A42 * (F2 + F3))
        F4 = f.f2(U4, p, t + dt) + Au
        # stage 5
        U5 = uprev + dt * (A51 * F1 + A52 * (F2 + F3) + A54 * F4)
        F5 = f.f2(U5, p, t + halfdt) + Au
        integrator.stats.nf2 += 4
        # update u
        u = uprev + dt * (B1 * F1 + B4 * F4 + B5 * F5)
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::HochOst4Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, rtmp2, Au, F2, F3, F4, F5, J, KsCache) = cache
    A = isa(f, SplitFunction) ? f.f1.f : (calc_J!(J, integrator, cache); J) # get linear operator
    alg = unwrap_alg(integrator, true)

    F1 = integrator.fsalfirst
    mul!(Au, A, uprev)
    halfdt = dt / 2
    if alg.krylov
        Ks, phiv_cache, ws, herm = KsCache
        w1_half, w2_half, w3_half, w4_half, w1, w2, w3, w4, w5 = ws
        kwargs = _arnoldi_kwargs(alg, A, integrator, herm)
        # Krylov on F1 (first column)
        arnoldi!(Ks, A, F1; kwargs...)
        phiv!(w1_half, halfdt, Ks, 3; cache = phiv_cache)
        phiv!(w1, dt, Ks, 3; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev + halfdt * w1_half[:, 2] # tmp is U2
        _compute_nl!(F2, f, tmp, p, t + halfdt, A, rtmp)
        F2 .+= Au
        # Krylov on F2 (second column)
        arnoldi!(Ks, A, F2; kwargs...)
        phiv!(w2_half, halfdt, Ks, 3; cache = phiv_cache)
        phiv!(w2, dt, Ks, 3; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev +
            dt * (
            0.5w1_half[:, 2] - w1_half[:, 3] +
                w2_half[:, 3]
        ) # tmp is U3
        _compute_nl!(F3, f, tmp, p, t + halfdt, A, rtmp)
        F3 .+= Au
        # Krylov on F3 (third column)
        arnoldi!(Ks, A, F3; kwargs...)
        phiv!(w3_half, halfdt, Ks, 3; cache = phiv_cache)
        phiv!(w3, dt, Ks, 3; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev +
            dt * (
            w1[:, 2] - 2w1[:, 3] + w2[:, 3] +
                w3[:, 3]
        ) # tmp is U4
        _compute_nl!(F4, f, tmp, p, t + dt, A, rtmp)
        F4 .+= Au
        # Krylov on F4 (fourth column)
        arnoldi!(Ks, A, F4; kwargs...)
        phiv!(w4_half, halfdt, Ks, 3; cache = phiv_cache)
        phiv!(w4, dt, Ks, 3; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev +
            dt *
            (
            0.5w1_half[:, 2] - 0.75w1_half[:, 3] +
                0.5w1_half[:, 4] + w1[:, 4] -
                0.25w1[:, 3] +
                0.5w2_half[:, 3] - w2[:, 4] +
                0.25w2[:, 3] - 0.5w2_half[:, 4] +
                0.5w3_half[:, 3] - w2[:, 4] +
                0.25w3[:, 3] - 0.5w3_half[:, 4] +
                w4[:, 4] - 0.25w4[:, 3] -
                0.25w4_half[:, 3] + 0.5w4_half[:, 4]
        ) # tmp is U5
        _compute_nl!(F5, f, tmp, p, t + dt, A, rtmp)
        F5 .+= Au
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 4
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
        end
        # Krylov on F5 (fifth column)
        arnoldi!(Ks, A, F5; kwargs...)
        phiv!(w5, dt, Ks, 3; cache = phiv_cache)
        # update u
        @muladd @.. broadcast = false @views rtmp = w1[:, 2] - 3w1[:, 3] + 4w1[:, 4] -
            w4[:, 3] + 4w4[:, 4] + 4w5[:, 3] -
            8w5[:, 4]
        @muladd @.. broadcast = false u = uprev + dt * rtmp
    else
        A21, A31, A32, A41, A42, A51, A52, A54, B1, B4, B5 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        mul!(rtmp, A21, F1)
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U2
        f.f2(F2, tmp, p, t + halfdt)
        F2 .+= Au
        # stage 3
        mul!(rtmp, A31, F1)
        mul!(rtmp2, A32, F2)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U3
        f.f2(F3, tmp, p, t + halfdt)
        F3 .+= Au
        # stage 4
        F2 .+= F3 # F2 now stores F2 + F3
        mul!(rtmp, A41, F1)
        mul!(rtmp2, A42, F2)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U4
        f.f2(F4, tmp, p, t + dt)
        F4 .+= Au
        # stage 5
        mul!(rtmp, A51, F1)
        mul!(rtmp2, A52, F2)
        rtmp .+= rtmp2
        mul!(rtmp2, A54, F4)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U5
        f.f2(F5, tmp, p, t + halfdt)
        F5 .+= Au
        integrator.stats.nf2 += 4
        # update u
        mul!(rtmp, B1, F1)
        mul!(rtmp2, B4, F4)
        rtmp .+= rtmp2
        mul!(rtmp2, B5, F5)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false u = uprev + dt * rtmp
    end

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::FriedliConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = isa(f, SplitFunction) ? f.f1.f : calc_J(integrator, cache) # get linear operator
    alg = unwrap_alg(integrator, true)

    Au = A * uprev
    F1 = integrator.fsalfirst
    halfdt = dt / 2
    if alg.krylov
        kwargs = (
            m = min(alg.m, size(A, 1)), opnorm = integrator.opts.internalopnorm,
            iop = alg.iop,
        )
        # Krylov on F1 (first column)
        Ks = arnoldi(A, F1; kwargs...)
        w1_half = phiv(halfdt, Ks, 2)
        w1 = phiv(dt, Ks, 3)
        U2 = uprev + halfdt * w1_half[:, 2]
        F2 = _compute_nl(f, U2, p, t + halfdt, A) + Au
        # Krylov on F2 (second column)
        Ks = arnoldi(A, F2; kwargs...)
        w2_half = phiv(halfdt, Ks, 2)
        w2 = phiv(dt, Ks, 2)
        U3 = uprev + halfdt * (w1_half[:, 2] - w1_half[:, 3] + w2_half[:, 3])
        F3 = _compute_nl(f, U3, p, t + halfdt, A) + Au
        # Krylov on F3 (third column)
        w3 = phiv(dt, A, F3, 3; kwargs...)
        U4 = uprev +
            dt * (
            w1[:, 2] - 2w1[:, 3] - (26 // 25) * w2[:, 2] + (2 // 25) * w2[:, 3] +
                (26 // 25) * w3[:, 2] + (48 // 25) * w3[:, 3]
        )
        F4 = _compute_nl(f, U4, p, t + dt, A) + Au
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 3
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
        end
        # Krylov on F4 (fourth column)
        w4 = phiv(dt, A, F4, 3; kwargs...)
        # update u
        u = uprev +
            dt * (
            w1[:, 2] - 3w1[:, 3] + 4w1[:, 4] + 4w3[:, 3] - 8w3[:, 4] - w4[:, 3] +
                4w4[:, 4]
        )
    else
        A21, A31, A32, A41, A42, A43, B1, B3, B4 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        U2 = uprev + dt * (A21 * F1)
        F2 = f.f2(U2, p, t + halfdt) + Au
        # stage 3
        U3 = uprev + dt * (A31 * F1 + A32 * F2)
        F3 = f.f2(U3, p, t + halfdt) + Au
        # stage 4
        U4 = uprev + dt * (A41 * F1 + A42 * F2 + A43 * F3)
        F4 = f.f2(U4, p, t + dt) + Au
        integrator.stats.nf2 += 3
        # update u
        u = uprev + dt * (B1 * F1 + B3 * F3 + B4 * F4)
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::FriedliCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, rtmp2, Au, F2, F3, F4, J, KsCache) = cache
    A = isa(f, SplitFunction) ? f.f1.f : (calc_J!(J, integrator, cache); J) # get linear operator
    alg = unwrap_alg(integrator, true)

    F1 = integrator.fsalfirst
    mul!(Au, A, uprev)
    halfdt = dt / 2
    if alg.krylov
        Ks, phiv_cache, ws = KsCache
        w1_half, w2_half, w1, w2, w3, w4 = ws
        kwargs = (
            m = min(alg.m, size(A, 1)), opnorm = integrator.opts.internalopnorm,
            iop = alg.iop,
        )
        # Krylov on F1 (first column)
        arnoldi!(Ks, A, F1; kwargs...)
        phiv!(w1_half, halfdt, Ks, 2; cache = phiv_cache)
        phiv!(w1, dt, Ks, 3; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev + halfdt * w1_half[:, 2] # tmp is U2
        _compute_nl!(F2, f, tmp, p, t + halfdt, A, rtmp)
        F2 .+= Au
        # Krylov on F2 (second column)
        arnoldi!(Ks, A, F2; kwargs...)
        phiv!(w2_half, halfdt, Ks, 2; cache = phiv_cache)
        phiv!(w2, dt, Ks, 2; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev +
            halfdt * (w1_half[:, 2] - w1_half[:, 3] + w2_half[:, 3]) # tmp is U3
        _compute_nl!(F3, f, tmp, p, t + halfdt, A, rtmp)
        F3 .+= Au
        # Krylov on F3 (third column)
        arnoldi!(Ks, A, F3; kwargs...)
        phiv!(w3, dt, Ks, 3; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev +
            dt * (
            w1[:, 2] - 2w1[:, 3] - (26 // 25) * w2[:, 2] + (2 // 25) * w2[:, 3] +
                (26 // 25) * w3[:, 2] + (48 // 25) * w3[:, 3]
        ) # tmp is U4
        _compute_nl!(F4, f, tmp, p, t + dt, A, rtmp)
        F4 .+= Au
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 3
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 3)
        end
        # Krylov on F4 (fourth column)
        arnoldi!(Ks, A, F4; kwargs...)
        phiv!(w4, dt, Ks, 3; cache = phiv_cache)
        # update u
        @muladd @.. broadcast = false @views rtmp = w1[:, 2] - 3w1[:, 3] + 4w1[:, 4] +
            4w3[:, 3] - 8w3[:, 4] - w4[:, 3] +
            4w4[:, 4]
        @muladd @.. broadcast = false u = uprev + dt * rtmp
    else
        A21, A31, A32, A41, A42, A43, B1, B3, B4 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        mul!(rtmp, A21, F1)
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U2
        f.f2(F2, tmp, p, t + halfdt)
        F2 .+= Au
        # stage 3
        mul!(rtmp, A31, F1)
        mul!(rtmp2, A32, F2)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U3
        f.f2(F3, tmp, p, t + halfdt)
        F3 .+= Au
        # stage 4
        mul!(rtmp, A41, F1)
        mul!(rtmp2, A42, F2)
        rtmp .+= rtmp2
        mul!(rtmp2, A43, F3)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U4
        f.f2(F4, tmp, p, t + dt)
        F4 .+= Au
        integrator.stats.nf2 += 3
        # update u
        mul!(rtmp, B1, F1)
        mul!(rtmp2, B3, F3)
        rtmp .+= rtmp2
        mul!(rtmp2, B4, F4)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false u = uprev + dt * rtmp
    end

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

#############################################
# Embedded ExpRK pair integrators

function perform_step!(integrator, cache::ERK32ZBConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = isa(f, SplitFunction) ? f.f1.f : calc_J(integrator, cache) # get linear operator
    alg = unwrap_alg(integrator, true)
    update_expRK_operators!(cache, alg, f, dt)

    Au = A * uprev
    F1 = integrator.fsalfirst
    halfdt = dt / 2
    threequarterdt = 3dt / 4
    nf2 = 2
    if alg.krylov
        kwargs = _arnoldi_kwargs(alg, A, integrator, nothing)
        # Krylov on F1 (first column)
        Ks = arnoldi(A, F1; kwargs...)
        w1 = phiv(dt, Ks, 3)
        w1h = phiv(halfdt, Ks, 3)
        w1q = phiv(threequarterdt, Ks, 2)
        U2 = uprev + dt * ((1 // 2) * w1h[:, 2])
        F2 = _compute_nl(f, U2, p, t + halfdt, A) + Au
        # Krylov on F2 (second column)
        Ks = arnoldi(A, F2; kwargs...)
        w2 = phiv(dt, Ks, 3)
        w2h = phiv(halfdt, Ks, 3)
        w2q = phiv(threequarterdt, Ks, 2)
        U3 = uprev +
            dt * (
            -(3 // 8) * w1h[:, 3] + (3 // 4) * w1q[:, 2] - (9 // 8) * w1q[:, 3] +
                (3 // 8) * w2h[:, 3] + (9 // 8) * w2q[:, 3]
        )
        F3 = _compute_nl(f, U3, p, t + threequarterdt, A) + Au
        # Krylov on F3 (third column)
        Ks = arnoldi(A, F3; kwargs...)
        w3 = phiv(dt, Ks, 3)
        # update u with the third order estimate
        u = uprev +
            dt * (
            w1[:, 2] - (19 // 12) * w1[:, 3] + (1 // 12) * w1[:, 4] +
                (3 // 4) * w2[:, 3] - (1 // 4) * w2[:, 4] +
                (5 // 6) * w3[:, 3] + (1 // 6) * w3[:, 4]
        )
        if integrator.opts.adaptive
            # the second order estimate uses the nonlinearity at the new solution
            w3h = phiv(halfdt, Ks, 3)
            w3q = phiv(threequarterdt, Ks, 2)
            F4 = _compute_nl(f, u, p, t + dt, A) + Au
            nf2 += 1
            Ks = arnoldi(A, F4; kwargs...)
            w4 = phiv(dt, Ks, 3)
            w4h = phiv(halfdt, Ks, 3)
            w4q = phiv(threequarterdt, Ks, 2)
            utilde = dt * (
                -(11 // 18) * w1[:, 2] - (7 // 3) * w1[:, 3] +
                    (1353 // 70) * w1[:, 4] - (9 // 14) * w1h[:, 2] -
                    (1 // 12) * w1h[:, 3] - (11 // 30) * w1h[:, 4] -
                    (7 // 6) * w1q[:, 2] - (2 // 7) * w1q[:, 3] +
                    (1 // 9) * w2[:, 2] + (5 // 4) * w2[:, 3] -
                    (5 // 12) * w2[:, 4] + (1 // 3) * w2h[:, 3] -
                    (1 // 6) * w2h[:, 4] + (1 // 6) * w2q[:, 2] +
                    (1 // 7) * w2q[:, 3] -
                    (2 // 3) * w3[:, 2] + (1 // 2) * w3[:, 3] +
                    (1 // 6) * w3[:, 4] + (1 // 7) * w3h[:, 2] +
                    (1 // 5) * w3h[:, 4] + (1 // 2) * w3q[:, 2] +
                    (1 // 7) * w3q[:, 3] +
                    (7 // 6) * w4[:, 2] + (7 // 12) * w4[:, 3] -
                    (2671 // 140) * w4[:, 4] + (1 // 2) * w4h[:, 2] -
                    (1 // 4) * w4h[:, 3] + (1 // 3) * w4h[:, 4] +
                    (1 // 2) * w4q[:, 2]
            )
            _expRK_EEst!(integrator, utilde, u)
        end
    else
        A21, A31, A32, B1, B2, B3, E1, E2, E3, E4 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        U2 = uprev + dt * (A21 * F1)
        F2 = f.f2(U2, p, t + halfdt) + Au
        # stage 3
        U3 = uprev + dt * (A31 * F1 + A32 * F2)
        F3 = f.f2(U3, p, t + threequarterdt) + Au
        # update u with the third order estimate
        u = uprev + dt * (B1 * F1 + B2 * F2 + B3 * F3)
        if integrator.opts.adaptive
            # the second order estimate uses the nonlinearity at the new solution
            F4 = f.f2(u, p, t + dt) + Au
            nf2 += 1
            utilde = dt * (E1 * F1 + E2 * F2 + E3 * F3 + E4 * F4)
            _expRK_EEst!(integrator, utilde, u)
        end
    end
    if isa(f, SplitFunction)
        integrator.stats.nf2 += nf2
    else
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, nf2)
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::ERK32ZBCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; utilde, tmp, rtmp, rtmp2, Au, F2, F3, F4, J, KsCache) = cache
    A = isa(f, SplitFunction) ? f.f1.f : (calc_J!(J, integrator, cache); J) # get linear operator
    alg = unwrap_alg(integrator, true)
    update_expRK_operators!(cache, alg, f, dt)

    F1 = integrator.fsalfirst
    mul!(Au, A, uprev)
    halfdt = dt / 2
    threequarterdt = 3dt / 4
    nf2 = 2
    if alg.krylov
        Ks, phiv_cache, ws, herm = KsCache
        w1, w1h, w1q, w2, w2h, w2q, w3, w3h, w3q, w4, w4h, w4q = ws
        kwargs = _arnoldi_kwargs(alg, A, integrator, herm)
        # Krylov on F1 (first column)
        arnoldi!(Ks, A, F1; kwargs...)
        phiv!(w1, dt, Ks, 3; cache = phiv_cache)
        phiv!(w1h, halfdt, Ks, 3; cache = phiv_cache)
        phiv!(w1q, threequarterdt, Ks, 2; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev + dt * ((1 // 2) * w1h[:, 2]) # tmp is U2
        _compute_nl!(F2, f, tmp, p, t + halfdt, A, rtmp)
        F2 .+= Au
        # Krylov on F2 (second column)
        arnoldi!(Ks, A, F2; kwargs...)
        phiv!(w2, dt, Ks, 3; cache = phiv_cache)
        phiv!(w2h, halfdt, Ks, 3; cache = phiv_cache)
        phiv!(w2q, threequarterdt, Ks, 2; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev +
            dt * (
            -(3 // 8) * w1h[:, 3] + (3 // 4) * w1q[:, 2] - (9 // 8) * w1q[:, 3] +
                (3 // 8) * w2h[:, 3] + (9 // 8) * w2q[:, 3]
        ) # tmp is U3
        _compute_nl!(F3, f, tmp, p, t + threequarterdt, A, rtmp)
        F3 .+= Au
        # Krylov on F3 (third column)
        arnoldi!(Ks, A, F3; kwargs...)
        phiv!(w3, dt, Ks, 3; cache = phiv_cache)
        # update u with the third order estimate
        @muladd @.. broadcast = false @views u = uprev +
            dt * (
            w1[:, 2] - (19 // 12) * w1[:, 3] + (1 // 12) * w1[:, 4] +
                (3 // 4) * w2[:, 3] - (1 // 4) * w2[:, 4] +
                (5 // 6) * w3[:, 3] + (1 // 6) * w3[:, 4]
        )
        if integrator.opts.adaptive
            # the second order estimate uses the nonlinearity at the new solution
            phiv!(w3h, halfdt, Ks, 3; cache = phiv_cache)
            phiv!(w3q, threequarterdt, Ks, 2; cache = phiv_cache)
            _compute_nl!(F4, f, u, p, t + dt, A, rtmp)
            F4 .+= Au
            nf2 += 1
            arnoldi!(Ks, A, F4; kwargs...)
            phiv!(w4, dt, Ks, 3; cache = phiv_cache)
            phiv!(w4h, halfdt, Ks, 3; cache = phiv_cache)
            phiv!(w4q, threequarterdt, Ks, 2; cache = phiv_cache)
            @muladd @.. broadcast = false @views utilde = dt * (
                -(11 // 18) * w1[:, 2] - (7 // 3) * w1[:, 3] +
                    (1353 // 70) * w1[:, 4] - (9 // 14) * w1h[:, 2] -
                    (1 // 12) * w1h[:, 3] - (11 // 30) * w1h[:, 4] -
                    (7 // 6) * w1q[:, 2] - (2 // 7) * w1q[:, 3] +
                    (1 // 9) * w2[:, 2] + (5 // 4) * w2[:, 3] -
                    (5 // 12) * w2[:, 4] + (1 // 3) * w2h[:, 3] -
                    (1 // 6) * w2h[:, 4] + (1 // 6) * w2q[:, 2] +
                    (1 // 7) * w2q[:, 3] -
                    (2 // 3) * w3[:, 2] + (1 // 2) * w3[:, 3] +
                    (1 // 6) * w3[:, 4] + (1 // 7) * w3h[:, 2] +
                    (1 // 5) * w3h[:, 4] + (1 // 2) * w3q[:, 2] +
                    (1 // 7) * w3q[:, 3] +
                    (7 // 6) * w4[:, 2] + (7 // 12) * w4[:, 3] -
                    (2671 // 140) * w4[:, 4] + (1 // 2) * w4h[:, 2] -
                    (1 // 4) * w4h[:, 3] + (1 // 3) * w4h[:, 4] +
                    (1 // 2) * w4q[:, 2]
            )
        end
    else
        A21, A31, A32, B1, B2, B3, E1, E2, E3, E4 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        mul!(rtmp, A21, F1)
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U2
        f.f2(F2, tmp, p, t + halfdt)
        F2 .+= Au
        # stage 3
        mul!(rtmp, A31, F1)
        mul!(rtmp2, A32, F2)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U3
        f.f2(F3, tmp, p, t + threequarterdt)
        F3 .+= Au
        # update u with the third order estimate
        mul!(rtmp, B1, F1)
        mul!(rtmp2, B2, F2)
        rtmp .+= rtmp2
        mul!(rtmp2, B3, F3)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false u = uprev + dt * rtmp
        if integrator.opts.adaptive
            # the second order estimate uses the nonlinearity at the new solution
            f.f2(F4, u, p, t + dt)
            F4 .+= Au
            nf2 += 1
            mul!(rtmp, E1, F1)
            mul!(rtmp2, E2, F2)
            rtmp .+= rtmp2
            mul!(rtmp2, E3, F3)
            rtmp .+= rtmp2
            mul!(rtmp2, E4, F4)
            rtmp .+= rtmp2
            @.. broadcast = false utilde = dt * rtmp
        end
    end
    if isa(f, SplitFunction)
        integrator.stats.nf2 += nf2
    else
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, nf2)
    end

    if integrator.opts.adaptive
        calculate_residuals!(
            tmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(tmp, t))
    end

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::ERK43ZBConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = isa(f, SplitFunction) ? f.f1.f : calc_J(integrator, cache) # get linear operator
    alg = unwrap_alg(integrator, true)
    update_expRK_operators!(cache, alg, f, dt)

    Au = A * uprev
    F1 = integrator.fsalfirst
    halfdt = dt / 2
    sixthdt = dt / 6
    if alg.krylov
        kwargs = _arnoldi_kwargs(alg, A, integrator, nothing)
        # Krylov on F1 (first column)
        Ks = arnoldi(A, F1; kwargs...)
        w1 = phiv(dt, Ks, 3)
        w1h = phiv(halfdt, Ks, 3)
        w1s = phiv(sixthdt, Ks, 2)
        U2 = uprev + dt * ((1 // 6) * w1s[:, 2])
        F2 = _compute_nl(f, U2, p, t + sixthdt, A) + Au
        # Krylov on F2 (second column)
        Ks = arnoldi(A, F2; kwargs...)
        w2 = phiv(dt, Ks, 3)
        w2h = phiv(halfdt, Ks, 3)
        w2s = phiv(sixthdt, Ks, 2)
        U3 = uprev +
            dt * (
            (1 // 2) * w1h[:, 2] - (3 // 2) * w1h[:, 3] - (1 // 2) * w1s[:, 3] +
                (3 // 2) * w2h[:, 3] + (1 // 2) * w2s[:, 3]
        )
        F3 = _compute_nl(f, U3, p, t + halfdt, A) + Au
        # Krylov on F3 (third column)
        Ks = arnoldi(A, F3; kwargs...)
        w3 = phiv(dt, Ks, 3)
        w3h = phiv(halfdt, Ks, 3)
        w3s = phiv(sixthdt, Ks, 2)
        U4 = uprev +
            dt * (
            -(19 // 90) * w1[:, 2] + (1 // 6) * w1h[:, 2] - (11 // 6) * w1h[:, 3] -
                (2 // 5) * w1h[:, 4] - (1 // 3) * w1s[:, 2] - (41 // 18) * w1s[:, 3] +
                (19 // 60) * w2[:, 2] + (1 // 2) * w2h[:, 2] + 2 * w2h[:, 3] +
                (3 // 5) * w2h[:, 4] + (1 // 2) * w2s[:, 2] + (13 // 6) * w2s[:, 3] -
                (19 // 180) * w3[:, 2] - (1 // 6) * w3h[:, 2] - (1 // 6) * w3h[:, 3] -
                (1 // 5) * w3h[:, 4] - (1 // 6) * w3s[:, 2] + (1 // 9) * w3s[:, 3]
        )
        F4 = _compute_nl(f, U4, p, t + halfdt, A) + Au
        # Krylov on F4 (fourth column)
        Ks = arnoldi(A, F4; kwargs...)
        w4 = phiv(dt, Ks, 3)
        w4h = phiv(halfdt, Ks, 3)
        # stage 5 is the third order estimate
        U5 = uprev +
            dt * (
            (71 // 90) * w1[:, 2] - 8 * w1[:, 3] + 24 * w1[:, 4] -
                (1 // 3) * w1h[:, 2] - (7 // 3) * w1h[:, 3] + (58 // 5) * w1h[:, 4] -
                (1 // 3) * w1s[:, 2] + (2 // 9) * w1s[:, 3] +
                (19 // 60) * w2[:, 2] + 9 * w2[:, 3] - 36 * w2[:, 4] +
                (1 // 2) * w2h[:, 2] + (7 // 2) * w2h[:, 3] - (87 // 5) * w2h[:, 4] +
                (1 // 2) * w2s[:, 2] - (1 // 3) * w2s[:, 3] -
                (19 // 180) * w3[:, 2] - 2 * w3[:, 3] + 18 * w3[:, 4] -
                (1 // 6) * w3h[:, 2] - (13 // 6) * w3h[:, 3] + (44 // 5) * w3h[:, 4] -
                (1 // 6) * w3s[:, 2] + (1 // 9) * w3s[:, 3] +
                w4[:, 3] - 6 * w4[:, 4] + w4h[:, 3] - 3 * w4h[:, 4]
        )
        F5 = _compute_nl(f, U5, p, t + dt, A) + Au
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 4
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
        end
        # Krylov on F5 (fifth column)
        Ks = arnoldi(A, F5; kwargs...)
        w5 = phiv(dt, Ks, 3)
        # update u with the fourth order estimate
        u = uprev +
            dt * (
            w1[:, 2] - (67 // 9) * w1[:, 3] + (52 // 3) * w1[:, 4] +
                8 * w2[:, 3] - 24 * w2[:, 4] -
                (11 // 9) * w3[:, 3] + (26 // 3) * w3[:, 4] +
                (7 // 9) * w4[:, 3] - (10 // 3) * w4[:, 4] -
                (1 // 9) * w5[:, 3] + (4 // 3) * w5[:, 4]
        )
        if integrator.opts.adaptive
            utilde = dt * (
                (19 // 90) * w1[:, 2] + (5 // 9) * w1[:, 3] - (20 // 3) * w1[:, 4] +
                    (1 // 3) * w1h[:, 2] + (7 // 3) * w1h[:, 3] -
                    (58 // 5) * w1h[:, 4] + (1 // 3) * w1s[:, 2] -
                    (2 // 9) * w1s[:, 3] -
                    (19 // 60) * w2[:, 2] - w2[:, 3] + 12 * w2[:, 4] -
                    (1 // 2) * w2h[:, 2] - (7 // 2) * w2h[:, 3] +
                    (87 // 5) * w2h[:, 4] - (1 // 2) * w2s[:, 2] +
                    (1 // 3) * w2s[:, 3] +
                    (19 // 180) * w3[:, 2] + (7 // 9) * w3[:, 3] -
                    (28 // 3) * w3[:, 4] + (1 // 6) * w3h[:, 2] +
                    (13 // 6) * w3h[:, 3] - (44 // 5) * w3h[:, 4] +
                    (1 // 6) * w3s[:, 2] - (1 // 9) * w3s[:, 3] -
                    (2 // 9) * w4[:, 3] + (8 // 3) * w4[:, 4] - w4h[:, 3] +
                    3 * w4h[:, 4] -
                    (1 // 9) * w5[:, 3] + (4 // 3) * w5[:, 4]
            )
            _expRK_EEst!(integrator, utilde, u)
        end
    else
        A21, A31, A32, A41, A42, A43, A51, A52, A53, A54,
            B1, B2, B3, B4, B5, E1, E2, E3, E4, E5 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        U2 = uprev + dt * (A21 * F1)
        F2 = f.f2(U2, p, t + sixthdt) + Au
        # stage 3
        U3 = uprev + dt * (A31 * F1 + A32 * F2)
        F3 = f.f2(U3, p, t + halfdt) + Au
        # stage 4
        U4 = uprev + dt * (A41 * F1 + A42 * F2 + A43 * F3)
        F4 = f.f2(U4, p, t + halfdt) + Au
        # stage 5 is the third order estimate
        U5 = uprev + dt * (A51 * F1 + A52 * F2 + A53 * F3 + A54 * F4)
        F5 = f.f2(U5, p, t + dt) + Au
        integrator.stats.nf2 += 4
        # update u with the fourth order estimate
        u = uprev + dt * (B1 * F1 + B2 * F2 + B3 * F3 + B4 * F4 + B5 * F5)
        if integrator.opts.adaptive
            utilde = dt * (E1 * F1 + E2 * F2 + E3 * F3 + E4 * F4 + E5 * F5)
            _expRK_EEst!(integrator, utilde, u)
        end
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::ERK43ZBCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; utilde, tmp, rtmp, rtmp2, Au, F2, F3, F4, F5, J, KsCache) = cache
    A = isa(f, SplitFunction) ? f.f1.f : (calc_J!(J, integrator, cache); J) # get linear operator
    alg = unwrap_alg(integrator, true)
    update_expRK_operators!(cache, alg, f, dt)

    F1 = integrator.fsalfirst
    mul!(Au, A, uprev)
    halfdt = dt / 2
    sixthdt = dt / 6
    if alg.krylov
        Ks, phiv_cache, ws, herm = KsCache
        w1, w1h, w1s, w2, w2h, w2s, w3, w3h, w3s, w4, w4h, w5 = ws
        kwargs = _arnoldi_kwargs(alg, A, integrator, herm)
        # Krylov on F1 (first column)
        arnoldi!(Ks, A, F1; kwargs...)
        phiv!(w1, dt, Ks, 3; cache = phiv_cache)
        phiv!(w1h, halfdt, Ks, 3; cache = phiv_cache)
        phiv!(w1s, sixthdt, Ks, 2; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev + dt * ((1 // 6) * w1s[:, 2]) # tmp is U2
        _compute_nl!(F2, f, tmp, p, t + sixthdt, A, rtmp)
        F2 .+= Au
        # Krylov on F2 (second column)
        arnoldi!(Ks, A, F2; kwargs...)
        phiv!(w2, dt, Ks, 3; cache = phiv_cache)
        phiv!(w2h, halfdt, Ks, 3; cache = phiv_cache)
        phiv!(w2s, sixthdt, Ks, 2; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev +
            dt * (
            (1 // 2) * w1h[:, 2] - (3 // 2) * w1h[:, 3] - (1 // 2) * w1s[:, 3] +
                (3 // 2) * w2h[:, 3] + (1 // 2) * w2s[:, 3]
        ) # tmp is U3
        _compute_nl!(F3, f, tmp, p, t + halfdt, A, rtmp)
        F3 .+= Au
        # Krylov on F3 (third column)
        arnoldi!(Ks, A, F3; kwargs...)
        phiv!(w3, dt, Ks, 3; cache = phiv_cache)
        phiv!(w3h, halfdt, Ks, 3; cache = phiv_cache)
        phiv!(w3s, sixthdt, Ks, 2; cache = phiv_cache)
        @muladd @.. broadcast = false @views tmp = uprev +
            dt * (
            -(19 // 90) * w1[:, 2] + (1 // 6) * w1h[:, 2] - (11 // 6) * w1h[:, 3] -
                (2 // 5) * w1h[:, 4] - (1 // 3) * w1s[:, 2] - (41 // 18) * w1s[:, 3] +
                (19 // 60) * w2[:, 2] + (1 // 2) * w2h[:, 2] + 2 * w2h[:, 3] +
                (3 // 5) * w2h[:, 4] + (1 // 2) * w2s[:, 2] + (13 // 6) * w2s[:, 3] -
                (19 // 180) * w3[:, 2] - (1 // 6) * w3h[:, 2] - (1 // 6) * w3h[:, 3] -
                (1 // 5) * w3h[:, 4] - (1 // 6) * w3s[:, 2] + (1 // 9) * w3s[:, 3]
        ) # tmp is U4
        _compute_nl!(F4, f, tmp, p, t + halfdt, A, rtmp)
        F4 .+= Au
        # Krylov on F4 (fourth column)
        arnoldi!(Ks, A, F4; kwargs...)
        phiv!(w4, dt, Ks, 3; cache = phiv_cache)
        phiv!(w4h, halfdt, Ks, 3; cache = phiv_cache)
        # stage 5 is the third order estimate
        @muladd @.. broadcast = false @views tmp = uprev +
            dt * (
            (71 // 90) * w1[:, 2] - 8 * w1[:, 3] + 24 * w1[:, 4] -
                (1 // 3) * w1h[:, 2] - (7 // 3) * w1h[:, 3] + (58 // 5) * w1h[:, 4] -
                (1 // 3) * w1s[:, 2] + (2 // 9) * w1s[:, 3] +
                (19 // 60) * w2[:, 2] + 9 * w2[:, 3] - 36 * w2[:, 4] +
                (1 // 2) * w2h[:, 2] + (7 // 2) * w2h[:, 3] - (87 // 5) * w2h[:, 4] +
                (1 // 2) * w2s[:, 2] - (1 // 3) * w2s[:, 3] -
                (19 // 180) * w3[:, 2] - 2 * w3[:, 3] + 18 * w3[:, 4] -
                (1 // 6) * w3h[:, 2] - (13 // 6) * w3h[:, 3] + (44 // 5) * w3h[:, 4] -
                (1 // 6) * w3s[:, 2] + (1 // 9) * w3s[:, 3] +
                w4[:, 3] - 6 * w4[:, 4] + w4h[:, 3] - 3 * w4h[:, 4]
        ) # tmp is U5
        _compute_nl!(F5, f, tmp, p, t + dt, A, rtmp)
        F5 .+= Au
        if isa(f, SplitFunction)
            integrator.stats.nf2 += 4
        else
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 4)
        end
        # Krylov on F5 (fifth column)
        arnoldi!(Ks, A, F5; kwargs...)
        phiv!(w5, dt, Ks, 3; cache = phiv_cache)
        # update u with the fourth order estimate
        @muladd @.. broadcast = false @views u = uprev +
            dt * (
            w1[:, 2] - (67 // 9) * w1[:, 3] + (52 // 3) * w1[:, 4] +
                8 * w2[:, 3] - 24 * w2[:, 4] -
                (11 // 9) * w3[:, 3] + (26 // 3) * w3[:, 4] +
                (7 // 9) * w4[:, 3] - (10 // 3) * w4[:, 4] -
                (1 // 9) * w5[:, 3] + (4 // 3) * w5[:, 4]
        )
        if integrator.opts.adaptive
            @muladd @.. broadcast = false @views utilde = dt * (
                (19 // 90) * w1[:, 2] + (5 // 9) * w1[:, 3] - (20 // 3) * w1[:, 4] +
                    (1 // 3) * w1h[:, 2] + (7 // 3) * w1h[:, 3] -
                    (58 // 5) * w1h[:, 4] + (1 // 3) * w1s[:, 2] -
                    (2 // 9) * w1s[:, 3] -
                    (19 // 60) * w2[:, 2] - w2[:, 3] + 12 * w2[:, 4] -
                    (1 // 2) * w2h[:, 2] - (7 // 2) * w2h[:, 3] +
                    (87 // 5) * w2h[:, 4] - (1 // 2) * w2s[:, 2] +
                    (1 // 3) * w2s[:, 3] +
                    (19 // 180) * w3[:, 2] + (7 // 9) * w3[:, 3] -
                    (28 // 3) * w3[:, 4] + (1 // 6) * w3h[:, 2] +
                    (13 // 6) * w3h[:, 3] - (44 // 5) * w3h[:, 4] +
                    (1 // 6) * w3s[:, 2] - (1 // 9) * w3s[:, 3] -
                    (2 // 9) * w4[:, 3] + (8 // 3) * w4[:, 4] - w4h[:, 3] +
                    3 * w4h[:, 4] -
                    (1 // 9) * w5[:, 3] + (4 // 3) * w5[:, 4]
            )
        end
    else
        A21, A31, A32, A41, A42, A43, A51, A52, A53, A54,
            B1, B2, B3, B4, B5, E1, E2, E3, E4, E5 = cache.ops
        # stage 1 (fsaled)
        # stage 2
        mul!(rtmp, A21, F1)
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U2
        f.f2(F2, tmp, p, t + sixthdt)
        F2 .+= Au
        # stage 3
        mul!(rtmp, A31, F1)
        mul!(rtmp2, A32, F2)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U3
        f.f2(F3, tmp, p, t + halfdt)
        F3 .+= Au
        # stage 4
        mul!(rtmp, A41, F1)
        mul!(rtmp2, A42, F2)
        rtmp .+= rtmp2
        mul!(rtmp2, A43, F3)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U4
        f.f2(F4, tmp, p, t + halfdt)
        F4 .+= Au
        # stage 5 is the third order estimate
        mul!(rtmp, A51, F1)
        mul!(rtmp2, A52, F2)
        rtmp .+= rtmp2
        mul!(rtmp2, A53, F3)
        rtmp .+= rtmp2
        mul!(rtmp2, A54, F4)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is U5
        f.f2(F5, tmp, p, t + dt)
        F5 .+= Au
        integrator.stats.nf2 += 4
        # update u with the fourth order estimate
        mul!(rtmp, B1, F1)
        mul!(rtmp2, B2, F2)
        rtmp .+= rtmp2
        mul!(rtmp2, B3, F3)
        rtmp .+= rtmp2
        mul!(rtmp2, B4, F4)
        rtmp .+= rtmp2
        mul!(rtmp2, B5, F5)
        rtmp .+= rtmp2
        @muladd @.. broadcast = false u = uprev + dt * rtmp
        if integrator.opts.adaptive
            mul!(rtmp, E1, F1)
            mul!(rtmp2, E2, F2)
            rtmp .+= rtmp2
            mul!(rtmp2, E3, F3)
            rtmp .+= rtmp2
            mul!(rtmp2, E4, F4)
            rtmp .+= rtmp2
            mul!(rtmp2, E5, F5)
            rtmp .+= rtmp2
            @.. broadcast = false utilde = dt * rtmp
        end
    end

    if integrator.opts.adaptive
        calculate_residuals!(
            tmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(tmp, t))
    end

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

#############################################
# EPIRK integrators
function perform_step!(integrator, cache::Exp4ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    J = calc_J(integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(uprev) is fsaled
    ts = [dt / 3, 2dt / 3, dt]
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)),
    )

    # Krylov for f(uprev)
    B1 = [zero(f0) f0]
    K1 = phiv_timestep(ts, J, B1; kwargs...) # tϕ(tA)f0
    @inbounds for i in 1:3
        @views K1[:, i] ./= ts[i]
    end
    w4 = K1 * [-7 / 300, 97 / 150, -37 / 300]
    u4 = uprev + dt * w4
    d4 = f(u4, p, t + dt) - f0 - dt * (J * w4) # TODO: what should be the time?
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # Krylov for the first remainder d4
    B2 = [zero(d4) d4]
    K2 = phiv_timestep(ts, J, B2; kwargs...)
    @inbounds for i in 1:3
        @views K2[:, i] ./= ts[i]
    end
    w7 = K1 * [59 / 300, -7 / 75, 269 / 300] + K2 * [2 / 3, 2 / 3, 2 / 3]
    u7 = uprev + dt * w7
    d7 = f(u7, p, t + dt) - f0 - dt * (J * w7)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # Krylov for the second remainder d7
    B3 = [zero(d7) d7]
    k7 = phiv_timestep(ts[1], J, B3; kwargs...)
    k7 ./= ts[1]
    # Update u
    u = uprev + dt * (K1[:, 3] + K2[:, 1] - 4 / 3 * K2[:, 2] + K2[:, 3] + k7 / 6)

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::Exp4Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, rtmp2, K, J, B, KsCache, ts) = cache
    calc_J!(J, integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(u0) is fsaled
    ts[1] = dt / 3
    ts[2] = 2dt / 3
    ts[3] = dt
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)), caches = KsCache,
    )

    # Krylov for f(uprev)
    B[:, 2] .= f0
    phiv_timestep!(K, ts, J, B; kwargs...)
    @inbounds for i in 1:3
        @views K[:, i] ./= ts[i]
    end
    @views @.. broadcast = false rtmp = (-7 / 300) * K[:, 1] + (97 / 150) * K[:, 2] +
        (-37 / 300) * K[:, 3] # rtmp is now w4
    @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is now u4
    mul!(rtmp2, J, rtmp)
    f(rtmp, tmp, p, t + dt) # TODO: what should be the time?
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @muladd @.. broadcast = false @view(B[:, 2]) = rtmp - f0 - dt * rtmp2 # B[:,2] is now d4
    # Partially update entities that use k1, k2, k3
    @views @.. broadcast = false rtmp = (59 / 300) * K[:, 1] + (-7 / 75) * K[:, 2] +
        (269 / 300) * K[:, 3] # rtmp is now w7
    @muladd @.. broadcast = false u = uprev + dt * @view(K[:, 3])
    # Krylov for the first remainder d4
    phiv_timestep!(K, ts, J, B; kwargs...)
    @inbounds for i in 1:3
        @views K[:, i] ./= ts[i]
    end
    @views @.. broadcast = false rtmp2 = (2 / 3) * (K[:, 1] + K[:, 2] + K[:, 3])
    rtmp .+= rtmp2 # w7 fully updated
    @muladd @.. broadcast = false tmp = uprev + dt * rtmp # tmp is now u7
    mul!(rtmp2, J, rtmp)
    f(rtmp, tmp, p, t + dt) # TODO: what should be the time?
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @muladd @.. broadcast = false @view(B[:, 2]) = rtmp - f0 - dt * rtmp2 # B[:,2] is now d7
    # Partially update entities that use k4, k5, k6
    @views @.. broadcast = false rtmp = K[:, 1] + (-4 / 3) * K[:, 2] + K[:, 3]
    axpy!(dt, rtmp, u)
    # Krylov for the second remainder d7
    k7 = @view(K[:, 1])
    phiv_timestep!(k7, ts[1], J, B; kwargs...)
    k7 ./= ts[1]
    axpy!(dt / 6, k7, u)

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::EPIRK4s3AConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    J = calc_J(integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(uprev) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)),
    )

    # Compute U2 and U3 vertically
    K = phiv_timestep([dt / 2, 2dt / 3], J, [zero(f0) f0]; kwargs...)
    U2 = uprev + K[:, 1]
    U3 = uprev + K[:, 2]
    R2 = f(U2, p, t + dt / 2) - f0 - J * K[:, 1] # remainder of U2
    R3 = f(U3, p, t + 2dt / 3) - f0 - J * K[:, 2] # remainder of U3
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)

    # Update u (horizontally)
    B = zeros(eltype(f0), length(f0), 5)
    B[:, 2] = f0
    B[:, 4] = (32R2 - 13.5R3) / dt^2
    B[:, 5] = (-144R2 + 81R3) / dt^3
    u = uprev + phiv_timestep(dt, J, B; kwargs...)

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::EPIRK4s3ACache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, rtmp2, K, J, B, KsCache, ts) = cache
    calc_J!(J, integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(u0) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)), caches = KsCache,
    )

    # Compute U2 and U3 vertically
    B[:, 2] .= f0
    ts[1] = dt / 2
    ts[2] = 2dt / 3
    phiv_timestep!(K, ts, J, @view(B[:, 1:2]); kwargs...)
    ## U2 and R2
    @.. broadcast = false tmp = uprev + @view(K[:, 1]) # tmp is now U2
    f(rtmp, tmp, p, t + dt / 2)
    mul!(rtmp2, J, @view(K[:, 1]))
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R2
    B[:, 4] .= (32 / dt^2) .* rtmp
    B[:, 5] .= (-144 / dt^3) .* rtmp
    ## U3 and R3
    @.. broadcast = false tmp = uprev + @view(K[:, 2]) # tmp is now U3
    f(rtmp, tmp, p, t + 2dt / 3)
    mul!(rtmp2, J, @view(K[:, 2]))
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R3
    @views B[:, 4] .-= (13.5 / dt^2) .* rtmp
    @views B[:, 5] .+= (81 / dt^3) .* rtmp
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)

    # Update u
    du = @view(K[:, 1])
    phiv_timestep!(du, dt, J, B; kwargs...)
    @.. broadcast = false u = uprev + du

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::EPIRK4s3BConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    J = calc_J(integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(uprev) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)),
    )

    # Compute U2 and U3 vertically
    K = phiv_timestep([dt / 2, 3dt / 4], J, [zero(f0) zero(f0) f0]; kwargs...)
    @views K[:, 1] .*= 8 / (3 * dt)
    @views K[:, 2] .*= 16 / (9 * dt)
    U2 = uprev + K[:, 1]
    U3 = uprev + K[:, 2]
    R2 = f(U2, p, t + dt / 2) - f0 - J * K[:, 1] # remainder of U2
    R3 = f(U3, p, t + 3dt / 4) - f0 - J * K[:, 2] # remainder of U3
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)

    # Update u (horizontally)
    B = zeros(eltype(f0), length(f0), 5)
    B[:, 2] = f0
    B[:, 4] = (54R2 - 16R3) / dt^2
    B[:, 5] = (-324R2 + 144R3) / dt^3
    u = uprev + phiv_timestep(dt, J, B; kwargs...)

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::EPIRK4s3BCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, rtmp2, K, J, B, KsCache, ts) = cache
    calc_J!(J, integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(u0) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)), caches = KsCache,
    )

    # Compute U2 and U3 vertically
    fill!(@view(B[:, 2]), zero(eltype(B)))
    B[:, 3] .= f0
    ts[1] = dt / 2
    ts[2] = 3dt / 4
    phiv_timestep!(K, ts, J, @view(B[:, 1:3]); kwargs...)
    @views K[:, 1] .*= 8 / (3 * dt)
    @views K[:, 2] .*= 16 / (9 * dt)
    ## U2 and R2
    @.. broadcast = false tmp = uprev + @view(K[:, 1]) # tmp is now U2
    f(rtmp, tmp, p, t + dt / 2)
    mul!(rtmp2, J, @view(K[:, 1]))
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R2
    B[:, 4] .= (54 / dt^2) .* rtmp
    B[:, 5] .= (-324 / dt^3) .* rtmp
    ## U3 and R3
    @.. broadcast = false tmp = uprev + @view(K[:, 2]) # tmp is now U3
    f(rtmp, tmp, p, t + 3dt / 4)
    mul!(rtmp2, J, @view(K[:, 2]))
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R3
    @views B[:, 4] .-= (16 / dt^2) .* rtmp
    @views B[:, 5] .+= (144 / dt^3) .* rtmp
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 2)

    # Update u
    fill!(@view(B[:, 3]), zero(eltype(B)))
    B[:, 2] .= f0
    du = @view(K[:, 1])
    phiv_timestep!(du, dt, J, B; kwargs...)
    @.. broadcast = false u = uprev + du

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::EPIRK5s3ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    J = calc_J(integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(uprev) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)),
    )

    # Compute U2 horizontally
    B = fill(zero(eltype(f0)), length(f0), 4)
    B[:, 3] = (55 / (8 * dt)) * f0
    B[:, 4] = (-3025 / (192 * dt^2)) * f0
    k = phiv_timestep(48dt / 55, J, B; kwargs...)
    U2 = uprev + k
    R2 = f(U2, p, t + 48dt / 55) - f0 - J * k # remainder of U2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Compute U3 horizontally
    B[:, 2] = (53 / 5) * f0
    B[:, 3] = (-648 / (5 * dt)) * f0
    B[:, 4] = (2916 / (5 * dt^2)) * f0 + (32065 / (1152 * dt^2)) * R2
    k = phiv_timestep(4dt / 9, J, B; kwargs...)
    U3 = uprev + k
    R3 = f(U3, p, t + 4dt / 9) - f0 - J * k # remainder of U3
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Update u (horizontally)
    B = fill(zero(eltype(f0)), length(f0), 5)
    B[:, 2] = f0
    B[:, 4] = (-166375 / (61056 * dt^2)) * R2 + (2187 / (106 * dt^2)) * R3
    B[:, 5] = (499125 / (27136 * dt^3)) * R2 - (2187 / (106 * dt^3)) * R3
    u = uprev + phiv_timestep(dt, J, B; kwargs...)

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::EPIRK5s3Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, k, rtmp, rtmp2, J, B, KsCache) = cache
    calc_J!(J, integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(u0) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)), caches = KsCache,
    )

    # Compute U2 horizontally
    fill!(@view(B[:, 2]), zero(eltype(B)))
    B[:, 3] .= (55 / (8 * dt)) .* f0
    B[:, 4] .= (-3025 / (192 * dt^2)) .* f0
    phiv_timestep!(k, 48dt / 55, J, @view(B[:, 1:4]); kwargs...)
    ## Compute R2
    @.. broadcast = false tmp = uprev + k # tmp is now U2
    f(rtmp, tmp, p, t + 48dt / 55)
    mul!(rtmp2, J, k)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R2

    # Compute U3 horizontally
    B[:, 2] .= (53 / 5) .* f0
    B[:, 3] .= (-648 / (5 * dt)) .* f0
    B[:, 4] .= (2916 / (5 * dt^2)) .* f0 .+ (32065 / (1152 * dt^2)) .* rtmp
    phiv_timestep!(k, 4dt / 9, J, @view(B[:, 1:4]); kwargs...)
    ## Update B matrix using R2
    B[:, 2] .= f0
    fill!(@view(B[:, 3]), zero(eltype(B)))
    B[:, 4] .= (-166375 / (61056 * dt^2)) .* rtmp
    B[:, 5] .= (499125 / (27136 * dt^3)) .* rtmp
    ## Compute R3 and update B
    @.. broadcast = false tmp = uprev + k # tmp is now U3
    f(rtmp, tmp, p, t + 4dt / 9)
    mul!(rtmp2, J, k)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R3
    @views B[:, 4] .+= (2187 / (106 * dt^2)) .* rtmp
    @views B[:, 5] .-= (2187 / (106 * dt^3)) .* rtmp

    # Update u
    phiv_timestep!(k, dt, J, B; kwargs...)
    @.. broadcast = false u = uprev + k

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::EXPRB53s3ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    J = calc_J(integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(uprev) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)),
    )

    # Compute the first group for U2 and U3
    B = [zero(f0) f0]
    K = phiv_timestep([dt / 2, 9dt / 10], J, B; kwargs...)
    U2 = uprev + K[:, 1]
    R2 = f(U2, p, t + dt / 2) - f0 - J * K[:, 1] # remainder of U2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    U3 = uprev + K[:, 2] # partially

    # Compute the second group for U3
    B = [zero(R2) zero(R2) zero(R2) R2]
    K = phiv_timestep([dt / 2, 9dt / 10], J, B; kwargs...)
    U3 .+= 216 / (25 * dt^2) .* K[:, 1] + 8 / dt^2 .* K[:, 2]
    R3 = f(U3, p, t + 9dt / 10) - f0 - J * (U3 - uprev) # remainder of U3
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Compute the third group for u
    B = fill(zero(eltype(f0)), length(f0), 5)
    B[:, 2] = f0
    B[:, 4] = (18 / dt^2) * R2 - (250 / (81 * dt^2)) * R3
    B[:, 5] = (-60 / dt^3) * R2 + (500 / (27 * dt^3)) * R3
    u = uprev + phiv_timestep(dt, J, B; kwargs...)

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::EXPRB53s3Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, rtmp2, K, J, B, KsCache, ts) = cache
    calc_J!(J, integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(u0) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)), caches = KsCache,
    )

    # Compute the first group for U2 and U3
    B[:, 2] .= f0
    ts[1] = dt / 2
    ts[2] = 9dt / 10
    phiv_timestep!(K, ts, J, @view(B[:, 1:2]); kwargs...)
    ## U2 and R2
    @.. broadcast = false tmp = uprev + @view(K[:, 1]) # tmp is now U2
    f(rtmp, tmp, p, t + dt / 2)
    mul!(rtmp2, J, @view(K[:, 1]))
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R2
    @.. broadcast = false tmp = uprev + @view(K[:, 2]) # tmp is now U3 (partially)

    # Compute the second group for U3
    fill!(@view(B[:, 2]), zero(eltype(B)))
    B[:, 4] .= rtmp
    phiv_timestep!(K, ts, J, @view(B[:, 1:4]); kwargs...)
    ## Update B using R2
    B[:, 2] .= f0
    B[:, 4] .= (18 / dt^2) .* rtmp
    B[:, 5] .= (-60 / dt^3) .* rtmp
    ## U3 and R3
    @views tmp .+= 216 / (25 * dt^2) .* K[:, 1] .+ 8 / dt^2 .* K[:, 2] # tmp is now U3
    f(rtmp, tmp, p, t + 9dt / 10)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    tmp .-= uprev
    mul!(rtmp2, J, tmp)
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R3
    ## Update B using R3
    @views B[:, 4] .-= (250 / (81 * dt^2)) .* rtmp
    @views B[:, 5] .+= (500 / (27 * dt^3)) .* rtmp

    # Update u
    du = @view(K[:, 1])
    phiv_timestep!(du, dt, J, B; kwargs...)
    @.. broadcast = false u = uprev + du

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::EPIRK5P1ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    J = calc_J(integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(uprev) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)),
    )

    # Coefficients (scaling factors absorbed)
    g11 = 0.35129592695058193092 * dt
    g21 = 0.84405472011657126298 * dt # g22 = g32
    g31 = dt
    g32 = 0.71111095364366870359 * dt
    g33 = 0.62378111953371494809 * dt
    a22 = 2.37739153404186675585
    b2 = 1.78975267532362337976
    b3 = 9.35854650579261718128 / dt^2

    # Compute the first column (f0)
    B = [zero(f0) f0]
    K1 = phiv_timestep([g11, g21, g31], J, B; kwargs...)
    ## U1 and R1
    U1 = uprev + K1[:, 1]
    R1 = f(U1, p, t + g11) - f0 - J * K1[:, 1] # remainder of U1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Compute the second column (R1)
    B = [zero(R1) R1]
    k2 = phiv_timestep(g32, J, B; kwargs...)
    ## U2 and R2
    U2 = uprev + K1[:, 2] + a22 * k2
    R2 = f(U2, p, t + g21) - f0 - J * (U2 - uprev) # remainder of U2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Compute the third column (dR = R2 - 2R1)
    B = fill(zero(eltype(f0)), length(f0), 4)
    B[:, 4] = R2 - 2R1
    k3 = phiv_timestep(g33, J, B; kwargs...)
    ## Update state
    u = uprev + K1[:, 3] + b2 * k2 + b3 * k3

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::EPIRK5P1Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, rtmp2, K, J, B, KsCache, ts) = cache
    calc_J!(J, integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(u0) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)), caches = KsCache,
    )

    # Coefficients (scaling factors absorbed)
    g11 = 0.35129592695058193092 * dt
    g21 = 0.84405472011657126298 * dt # g22 = g32
    g31 = dt
    g32 = 0.71111095364366870359 * dt
    g33 = 0.62378111953371494809 * dt
    a22 = 2.37739153404186675585
    b2 = 1.78975267532362337976
    b3 = 9.35854650579261718128 / dt^2

    # Compute the first column (f0)
    B[:, 2] .= f0
    ts[1] = g11
    ts[2] = g21
    ts[3] = g31
    phiv_timestep!(K, ts, J, @view(B[:, 1:2]); kwargs...)
    ## U1 and R1
    @.. broadcast = false tmp = uprev + @view(K[:, 1]) # tmp is now U1
    f(rtmp, tmp, p, t + g11)
    mul!(rtmp2, J, @view(K[:, 1]))
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R1
    @.. broadcast = false tmp = uprev + @view(K[:, 2]) # partially update U2 (stored tmp)
    @.. broadcast = false u = uprev + @view(K[:, 3]) # partially update u
    B[:, 2] .= rtmp
    @.. broadcast = false @view(B[:, 4]) = (-2) * rtmp

    # Compute the second column (R1)
    k = @view(K[:, 1])
    phiv_timestep!(k, g32, J, @view(B[:, 1:2]); kwargs...)
    ## U2 and R2
    axpy!(a22, k, tmp) # tmp is now U2
    f(rtmp, tmp, p, t + g21)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    tmp .-= uprev
    mul!(rtmp2, J, tmp)
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R2
    axpy!(b2, k, u) # partially update u
    @views B[:, 4] .+= rtmp # is now dR

    # Compute the third column (dR = R2 - 2R1)
    fill!(@view(B[:, 2]), zero(eltype(B)))
    phiv_timestep!(k, g33, J, B; kwargs...)
    axpy!(b3, k, u)

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::EPIRK5P2ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    J = calc_J(integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(uprev) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)),
    )

    # Coefficients (scaling factors absorbed)
    g11 = 0.46629408528088195806 * dt
    g21 = 0.8821791265336386514 * dt # g22 = g32
    g31 = dt
    g32 = 0.92074916488140031449 * dt
    g33 = 0.79791561832664517267 * dt
    a22 = 2.80620373289331259751 / dt
    b2 = 2.52806310256246280783 / dt
    b31 = -0.128486782657005566142
    b32 = -0.161028033172809183351 / dt
    b33 = 5.26726331616909606251 / dt^2

    # Compute the first column (f0)
    B = [zero(f0) f0]
    K1 = phiv_timestep([g11, g21, g31], J, B; kwargs...)
    ## U1 and R1
    U1 = uprev + K1[:, 1]
    R1 = f(U1, p, t + g11) - f0 - J * K1[:, 1] # remainder of U1
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Compute the second column (R1)
    B = [zero(R1) zero(R1) R1]
    k2 = phiv_timestep(g32, J, B; kwargs...)
    ## U2 and R2
    U2 = uprev + K1[:, 2] + a22 * k2
    R2 = f(U2, p, t + g21) - f0 - J * (U2 - uprev) # remainder of U2
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)

    # Compute the third column (dR = R2 - 2R1)
    dR = R2 - 2R1
    B = [zero(dR) b31 * dR b32 * dR b33 * dR]
    k3 = phiv_timestep(g33, J, B; kwargs...)
    ## Update state
    u = uprev + K1[:, 3] + b2 * k2 + k3

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::EPIRK5P2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; tmp, rtmp, rtmp2, dR, K, J, B, KsCache, ts) = cache
    calc_J!(J, integrator, cache)
    alg = unwrap_alg(integrator, true)
    f0 = integrator.fsalfirst # f(u0) is fsaled
    kwargs = (
        tol = integrator.opts.reltol, iop = alg.iop,
        opnorm = integrator.opts.internalopnorm,
        adaptive = alg.adaptive_krylov, tau = (alg.adaptive_krylov ? 0.0 : dt),
        m = min(alg.m, size(J, 1)), caches = KsCache,
    )

    # Coefficients (scaling factors absorbed)
    g11 = 0.46629408528088195806 * dt
    g21 = 0.8821791265336386514 * dt # g22 = g32
    g31 = dt
    g32 = 0.92074916488140031449 * dt
    g33 = 0.79791561832664517267 * dt
    a22 = 2.80620373289331259751 / dt
    b2 = 2.52806310256246280783 / dt
    b31 = -0.128486782657005566142
    b32 = -0.161028033172809183351 / dt
    b33 = 5.26726331616909606251 / dt^2

    # Compute the first column (f0)
    B[:, 2] .= f0
    ts[1] = g11
    ts[2] = g21
    ts[3] = g31
    phiv_timestep!(K, ts, J, @view(B[:, 1:2]); kwargs...)
    ## U1 and R1
    @.. broadcast = false tmp = uprev + @view(K[:, 1]) # tmp is now U1
    f(rtmp, tmp, p, t + g11)
    mul!(rtmp2, J, @view(K[:, 1]))
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R1
    @.. broadcast = false tmp = uprev + @view(K[:, 2]) # partially update U2 (stored in tmp)
    @.. broadcast = false u = uprev + @view(K[:, 3]) # partially update u
    @.. broadcast = false dR = -2rtmp # partially update dR

    # Compute the second column (R1)
    fill!(@view(B[:, 2]), zero(eltype(B)))
    B[:, 3] .= rtmp
    k = @view(K[:, 1])
    phiv_timestep!(k, g32, J, @view(B[:, 1:3]); kwargs...)
    ## U2 and R2
    axpy!(a22, k, tmp) # tmp is now U2
    f(rtmp, tmp, p, t + g21)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    tmp .-= uprev
    mul!(rtmp2, J, tmp)
    @.. broadcast = false rtmp = rtmp - f0 - rtmp2 # rtmp is now R2
    dR .+= rtmp # dR is now R2 - 2R1
    axpy!(b2, k, u) # partially update u

    # Compute the third column (dR = R2 - 2R1)
    @.. broadcast = false @view(B[:, 2]) = b31 * dR
    @.. broadcast = false @view(B[:, 3]) = b32 * dR
    @.. broadcast = false @view(B[:, 4]) = b33 * dR
    phiv_timestep!(k, g33, J, B; kwargs...)
    u .+= k

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

######################################################
# Adaptive exponential Rosenbrock integrators
function perform_step!(integrator, cache::Exprb32ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = calc_J(integrator, cache) # get linear operator
    alg = unwrap_alg(integrator, true)

    F1 = integrator.fsalfirst
    w1 = phiv(
        dt, A, F1, 3; m = min(alg.m, size(A, 1)),
        opnorm = integrator.opts.internalopnorm, iop = alg.iop
    )
    U2 = uprev + dt * w1[:, 2]
    F2 = _compute_nl(f, U2, p, t + dt, A) + A * uprev
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    w2 = phiv(
        dt, A, F2, 3; m = min(alg.m, size(A, 1)),
        opnorm = integrator.opts.internalopnorm, iop = alg.iop
    )
    u = uprev + dt * (w1[:, 2] - 2w1[:, 4] + 2w2[:, 4])
    if integrator.opts.adaptive
        # error estimator for the imbedded method
        utilde = 2dt * (-w1[:, 4] + w2[:, 4])
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::Exprb32Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; utilde, tmp, rtmp, F2, J, KsCache) = cache
    calc_J!(J, integrator, cache)
    alg = unwrap_alg(integrator, true)

    F1 = integrator.fsalfirst
    Ks, phiv_cache, ws, herm = KsCache
    w1, w2 = ws
    kwargs = _arnoldi_kwargs(alg, J, integrator, herm)
    # Krylov for F1
    arnoldi!(Ks, J, F1; kwargs...)
    phiv!(w1, dt, Ks, 3; cache = phiv_cache)
    # Krylov for F2
    @muladd @.. broadcast = false tmp = uprev + dt * @view(w1[:, 2])
    _compute_nl!(F2, f, tmp, p, t + dt, J, rtmp)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    F2 .+= mul!(rtmp, J, uprev)
    arnoldi!(Ks, J, F2; kwargs...)
    phiv!(w2, dt, Ks, 3; cache = phiv_cache)
    # Update u
    u .= uprev
    axpy!(dt, @view(w1[:, 2]), u)
    axpy!(-2dt, @view(w1[:, 4]), u)
    axpy!(2dt, @view(w2[:, 4]), u)
    if integrator.opts.adaptive
        # error estimator for the imbedded method
        @views @.. broadcast = false utilde = (2 * dt) * (-w1[:, 4] + w2[:, 4])
        calculate_residuals!(
            tmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(tmp, t))
    end

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

function perform_step!(integrator, cache::Exprb43ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = calc_J(integrator, cache) # get linear operator
    alg = unwrap_alg(integrator, true)

    Au = A * uprev
    F1 = integrator.fsalfirst
    kwargs = _arnoldi_kwargs(alg, A, integrator, nothing)
    # Krylov on F1 (first column)
    Ks = arnoldi(A, F1; kwargs...)
    w1_half = phiv(dt / 2, Ks, 1)
    w1 = phiv(dt, Ks, 4)
    U2 = uprev + dt / 2 * w1_half[:, 2]
    F2 = _compute_nl(f, U2, p, t + dt / 2, A) + Au
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # Krylov on F2 (second column)
    w2 = phiv(dt, A, F2, 4; kwargs...)
    U3 = uprev + dt * w2[:, 2]
    F3 = _compute_nl(f, U3, p, t + dt, A) + Au
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # Krylov on F3 (third column)
    w3 = phiv(dt, A, F3, 4; kwargs...)
    u = uprev +
        dt * (
        w1[:, 2] - 14w1[:, 4] + 36w1[:, 5] + 16w2[:, 4] - 48w2[:, 5] - 2w3[:, 4] +
            12w3[:, 5]
    )
    if integrator.opts.adaptive
        # error estimator for the imbedded method
        utilde = dt * (36w1[:, 5] - 48w2[:, 5] + 12w3[:, 5])
        atmp = calculate_residuals(
            utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(atmp, t))
    end

    # Update integrator state
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::Exprb43Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; utilde, tmp, rtmp, Au, F2, F3, J, KsCache) = cache
    calc_J!(J, integrator, cache)
    alg = unwrap_alg(integrator, true)

    F1 = integrator.fsalfirst
    mul!(Au, J, uprev)
    halfdt = dt / 2
    Ks, phiv_cache, ws, herm = KsCache
    w1_half, w1, w2, w3 = ws
    kwargs = _arnoldi_kwargs(alg, J, integrator, herm)
    # Krylov for F1
    arnoldi!(Ks, J, F1; kwargs...)
    phiv!(w1_half, halfdt, Ks, 1; cache = phiv_cache)
    phiv!(w1, dt, Ks, 4; cache = phiv_cache)
    @muladd @.. broadcast = false @views tmp = uprev + halfdt * w1_half[:, 2] # tmp is U2
    _compute_nl!(F2, f, tmp, p, t + halfdt, J, rtmp)
    F2 .+= Au
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # Krylov for F2
    arnoldi!(Ks, J, F2; kwargs...)
    phiv!(w2, dt, Ks, 4; cache = phiv_cache)
    @muladd @.. broadcast = false @views tmp = uprev + dt * w2[:, 2] # tmp is U3
    _compute_nl!(F3, f, tmp, p, t + dt, J, rtmp)
    F3 .+= Au
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # Krylov for F3 (third column)
    arnoldi!(Ks, J, F3; kwargs...)
    phiv!(w3, dt, Ks, 4; cache = phiv_cache)
    # Update u
    @views @.. broadcast = false rtmp = w1[:, 2] - 14w1[:, 4] + 36w1[:, 5] + 16w2[:, 4] -
        48w2[:, 5] - 2w3[:, 4] + 12w3[:, 5]
    @muladd @.. broadcast = false u = uprev + dt * rtmp
    if integrator.opts.adaptive
        @views @.. broadcast = false rtmp = 36w1[:, 5] - 48w2[:, 5] + 12w3[:, 5]
        @.. broadcast = false utilde = dt * rtmp
        calculate_residuals!(
            tmp, utilde, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        OrdinaryDiffEqCore.set_EEst!(integrator, integrator.opts.internalnorm(tmp, t))
    end

    # Update integrator state
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    # integrator.k is automatically set due to aliasing
end

######################################################
# Multistep exponential integrators
function initialize!(integrator, cache::ETD2ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

    # Pre-start fsal
    lin = integrator.f.f1(integrator.uprev, integrator.p, integrator.t)
    nl = integrator.f.f2(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    nlprev = zero(nl) # to be computed in the first iteration via ETD1
    integrator.fsalfirst = ETD2Fsal(lin, nl, nlprev)

    # Avoid undefined entries if k is an array of arrays
    rate_prototype = lin
    integrator.fsallast = ETD2Fsal(rate_prototype)
    integrator.k[1] = lin + nl
    return integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator, cache::ETD2ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    (; lin, nl, nlprev) = integrator.fsalfirst
    (; exphA, phihA, B1, B0) = cache
    integrator.k[1] = lin + nl

    if integrator.iter == 1 # ETD1 for initial step
        u = exphA * uprev + dt * (phihA * nl)
    else
        u = exphA * uprev + dt * (B1 * nl + B0 * nlprev)
    end
    integrator.u = u

    # Push the fsal at t+dt
    nlprev = nl
    lin = f.f1(u, p, t + dt)
    nl = f.f2(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    integrator.k[2] = lin + nl
    integrator.fsallast.lin = lin
    integrator.fsallast.nl = nl
    return integrator.fsallast.nlprev = nlprev
end

function initialize!(integrator, cache::ETD2Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    rate_prototype = cache.rtmp1

    # Pre-start fsal
    (; lin, nl) = integrator.fsalfirst
    integrator.f.f1(lin, integrator.uprev, integrator.p, integrator.t)
    integrator.f.f2(nl, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1

    # Avoid undefined entries if k is an array of arrays
    integrator.k[1] = lin + nl
    return integrator.k[2] = zero(rate_prototype)
end

function perform_step!(integrator, cache::ETD2Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; lin, nl, nlprev) = integrator.fsalfirst
    (; utmp, rtmp1, rtmp2, exphA, phihA, B1, B0) = cache
    @.. broadcast = false integrator.k[1] = lin + nl

    if integrator.iter == 1 # ETD1 for initial step
        mul!(utmp, exphA, uprev)
        mul!(rtmp1, phihA, nl)
        @muladd @.. broadcast = false u = utmp + dt * rtmp1
    else
        mul!(utmp, exphA, uprev)
        mul!(rtmp1, B1, nl)
        mul!(rtmp2, B0, nlprev)
        @muladd @.. broadcast = false u = utmp + dt * (rtmp1 + rtmp2)
    end

    # Push the fsal at t+dt
    fsallast = integrator.fsallast
    fsallast.nlprev .= nl
    f.f1(fsallast.lin, u, p, t + dt)
    f.f2(fsallast.nl, u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.stats.nf2 += 1
    return @.. broadcast = false integrator.k[2] = fsallast.lin + fsallast.nl
end
