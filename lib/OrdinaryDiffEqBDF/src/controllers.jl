"""
    BDFController(; qmin, qmax, qsteady_min, qsteady_max, gamma, qmax_first_step,
                  failfactor)

Step-size controller for the variable-order BDF family (`QNDF`, `FBDF`,
`DFBDF`). Composes the standard step-size knobs via [`CommonControllerOptions`](@ref);
the adaptive logic is integrated into the algorithm itself, so the cache
falls through to alg-level dispatch (the same way the legacy
`DummyControllerCache` did) but exposes the knobs as a real, settable
controller. Pass it explicitly with `solve(prob, alg; controller = BDFController(...))`,
or rely on the default constructed by `default_controller(QT, alg)`.
"""
struct BDFController{B <: Union{NamedTuple, CommonControllerOptions}} <: AbstractController
    basic::B
end

BDFController(; kwargs...) = BDFController(NamedTuple(kwargs))

BDFController(alg; kwargs...) = BDFController(Float64, alg; kwargs...)
BDFController(::Type{QT}, alg; kwargs...) where {QT} =
    BDFController(resolve_basic(NamedTuple(kwargs), alg, QT))

mutable struct BDFControllerCache{T, E, C, NLPType} <: AbstractControllerCache
    controller::BDFController{CommonControllerOptions{T, NLPType}}
    cache::C
    EEst::E
end

function setup_controller_cache(
        alg::Union{QNDF, FBDF, DFBDF, MOOSE234}, cache, controller::BDFController, ::Type{E}, disco_probs
    ) where {E}
    QT = _resolved_QT(controller.basic)
    basic = resolve_basic(controller.basic, alg, QT; disco_probs)
    resolved = BDFController(basic)
    return BDFControllerCache{QT, E, typeof(cache), eltype(disco_probs)}(resolved, cache, oneunit(E))
end

# The BDF stepsize/accept/reject logic lives at the algorithm level — the
# controller cache just delegates back, mirroring how DummyControllerCache did.
@inline OrdinaryDiffEqCore.stepsize_controller!(integrator, ::BDFControllerCache, alg) =
    stepsize_controller!(integrator, alg)
@inline OrdinaryDiffEqCore.step_accept_controller!(integrator, ::BDFControllerCache, alg, q) =
    step_accept_controller!(integrator, alg, q)
@inline OrdinaryDiffEqCore.step_reject_controller!(integrator, ::BDFControllerCache, alg) =
    step_reject_controller!(integrator, alg)
@inline OrdinaryDiffEqCore.post_newton_controller!(integrator, ::BDFControllerCache, alg) =
    post_newton_controller!(integrator, alg)
@inline OrdinaryDiffEqCore.accept_step_controller(
    integrator, cache::BDFControllerCache, alg,
) = get_EEst(cache) <= 1

# Per-algorithm defaults for `CommonControllerOptions` knobs that don't match the
# generic IController defaults. These match the historical alg-struct kwargs
# `QNDF(qmax=5//1, qsteady_min=9//10, qsteady_max=12//10)` etc. and the formerly
# hard-coded `zₛ = 1.2` step-size safety factor in the BDF stepsize logic.
qmax_default(::Union{QNDF, FBDF, DFBDF}) = 5 // 1
qsteady_min_default(::Union{QNDF, QNDF1, QNDF2, FBDF, DFBDF}) = 9 // 10
qsteady_max_default(::Union{QNDF, QNDF1, QNDF2, FBDF, DFBDF}) = 12 // 10
gamma_default(::Union{QNDF, FBDF, DFBDF}) = 12 // 10

function default_controller(QT, alg::Union{QNDF, FBDF, DFBDF})
    # Thread the alg-level kwargs through to the CommonControllerOptions so that
    # `QNDF(qmax = 20)` keeps working (qmax = 20 ends up on the controller).
    return BDFController(
        QT, alg;
        qmax = alg.qmax,
        qsteady_min = alg.qsteady_min,
        qsteady_max = alg.qsteady_max,
    )
end

# MOOSE234's struct carries no qmax/qsteady fields, so build the delegating
# BDFController from the trait-based defaults instead of threading alg kwargs.
default_controller(QT, alg::MOOSE234) = BDFController(QT, alg)

# QNDF
stepsize_controller!(integrator, alg::QNDF) = nothing

# this stepsize and order controller is taken from
# Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s paper
# E. Alberdi Celaya, J. J. Anza Aguirrezabala, and P. Chatzipantelidis

function step_accept_controller!(integrator, alg::QNDF, q)
    return step_accept_controller!(integrator, integrator.cache, alg, q)
end

function step_accept_controller!(integrator, cache::Union{QNDFCache, QNDFConstantCache}, alg::QNDF{max_order}, q) where {max_order}
    #step is accepted, reset count of consecutive failed steps
    cache.consfailcnt = 0
    cache.nconsteps += 1
    if iszero(OrdinaryDiffEqCore.get_EEst(integrator))
        return integrator.dt * get_current_qmax(integrator, get_qmax(integrator))
    else
        est = OrdinaryDiffEqCore.get_EEst(integrator)
        estₖ₋₁ = cache.EEst1
        estₖ₊₁ = cache.EEst2
        h = integrator.dt
        k = cache.order
        prefer_const_step = cache.nconsteps < cache.order + 2
        zₛ = get_gamma(integrator)
        zᵤ = 0.1
        Fᵤ = 10
        expo = 1 / (k + 1)
        z = zₛ * ((est)^expo)
        F = inv(z)
        hₙ = h
        kₙ = k
        if z <= zᵤ
            hₖ = Fᵤ * h
        else
            hₖ = F * h
        end
        hₖ₋₁ = 0.0
        hₖ₊₁ = 0.0

        if k > 1
            expo = 1 / k
            zₖ₋₁ = 1.3 * ((estₖ₋₁)^expo)
            Fₖ₋₁ = inv(zₖ₋₁)
            if zₖ₋₁ <= 0.1
                hₖ₋₁ = 10 * h
            elseif zₖ₋₁ <= 1.3
                hₖ₋₁ = Fₖ₋₁ * h
            end
            if hₖ₋₁ > hₖ
                hₙ = hₖ₋₁
                kₙ = k - 1
            else
                hₙ = hₖ
                kₙ = k
            end
        else
            hₙ = hₖ
            kₙ = k
        end

        if k < max_order
            expo = 1 / (k + 2)
            zₖ₊₁ = 1.4 * ((estₖ₊₁)^expo)
            Fₖ₊₁ = inv(zₖ₊₁)

            if zₖ₊₁ <= 0.1
                hₖ₊₁ = 10 * h
            elseif 0.1 < zₖ₊₁ <= 1.4
                hₖ₊₁ = Fₖ₊₁ * h
            end
            if hₖ₊₁ > hₙ
                hₙ = hₖ₊₁
                kₙ = k + 1
            end
        end
        cache.order = kₙ
        q = integrator.dt / hₙ
    end
    if prefer_const_step
        if q < 1.2 && q > 0.6
            return integrator.dt
        end
    end
    if q <= get_qsteady_max(integrator) && q >= get_qsteady_min(integrator)
        return integrator.dt
    end
    return integrator.dt / q
end

function bdf_step_reject_controller!(integrator, cache, EEst1)
    k = cache.order
    h = integrator.dt
    cache.consfailcnt += 1
    cache.nconsteps = 0

    controller_cache = integrator.controller_cache
    discontinuity_detection = if controller_cache isa OrdinaryDiffEqCore.CompositeControllerCache
        current_idx = integrator.cache.current
        controller_cache.caches[current_idx].controller.basic.discontinuity_detection
    else
        controller_cache.controller.basic.discontinuity_detection
    end

    if discontinuity_detection
        disco_dt = set_discontinuity(integrator.u, integrator.uprev, integrator)
        if disco_dt != -1
            integrator.dt = disco_dt
            return integrator.dt
        end
    end

    if cache.consfailcnt > 1
        h = h / 2
    end
    zₛ = get_gamma(integrator)
    expo = 1 / (k + 1)
    z = zₛ * ((OrdinaryDiffEqCore.get_EEst(integrator))^expo)
    F = inv(z)
    if z <= 10
        hₖ = F * h
    else # z > 10
        hₖ = 0.1 * h
    end
    hₙ = hₖ
    kₙ = k
    if k > 1
        expo = 1 / k
        zₖ₋₁ = 1.3 * (EEst1^expo)
        Fₖ₋₁ = inv(zₖ₋₁)
        if zₖ₋₁ <= 10
            hₖ₋₁ = Fₖ₋₁ * h
        else # zₖ₋₁ > 10
            hₖ₋₁ = 0.1 * h
        end
        if cache.consfailcnt > 2 || hₖ₋₁ > hₖ
            hₙ = min(h, hₖ₋₁)
            kₙ = k - 1
        end
    end
    # Restart BDf (clear history) when we failed repeatedly
    if kₙ == 1 && cache.consfailcnt > 3
        derivative_discontinuity!(integrator, true)
    end
    integrator.dt = hₙ
    return cache.order = kₙ
end

function step_reject_controller!(integrator, alg::QNDF)
    return step_reject_controller!(integrator, integrator.cache, alg)
end
function step_reject_controller!(integrator, cache::Union{QNDFCache, QNDFConstantCache}, ::QNDF)
    return bdf_step_reject_controller!(integrator, cache, cache.EEst1)
end

# Without this method, falling through to the generic 2-arg
# `post_newton_controller!(integrator, alg)` in OrdinaryDiffEqCore recurses
# back into the `BDFControllerCache` 3-arg dispatch above, which calls the
# 2-arg again — a StackOverflowError on the first Newton failure.
function post_newton_controller!(integrator, alg::QNDF)
    return post_newton_controller!(integrator, integrator.cache, alg)
end
function post_newton_controller!(integrator, cache::Union{QNDFCache, QNDFConstantCache}, ::QNDF)
    integrator.dt = integrator.dt / get_failfactor(integrator)
    return nothing
end

function step_reject_controller!(integrator, alg::FBDF)
    return step_reject_controller!(integrator, integrator.cache, alg)
end
function step_reject_controller!(integrator, cache::Union{FBDFCache, FBDFConstantCache}, ::FBDF)
    return bdf_step_reject_controller!(integrator, cache, cache.terkm1)
end

function post_newton_controller!(integrator, alg::FBDF)
    return post_newton_controller!(integrator, integrator.cache, alg)
end
function post_newton_controller!(integrator, cache::Union{FBDFCache, FBDFConstantCache}, alg::FBDF)
    if cache.order > 1 && cache.nlsolver.nfails >= 3
        cache.order -= 1
    end
    integrator.dt = integrator.dt / get_failfactor(integrator)
    cache.consfailcnt += 1
    cache.nconsteps = 0
    return nothing
end

function choose_order!(
        alg::FBDF, integrator,
        cache::FBDFCache,
        ::Val{max_order}
    ) where {max_order}
    (; t, dt, u, uprev) = integrator
    (; atmp, ts_tmp, terkm2, terkm1, terk, terkp1, terk_tmp, u_history, fd_weights) = cache
    k = cache.order
    # Use CVODE-style qwait countdown: only consider order increase when qwait reaches 0
    if k < max_order && cache.qwait == 0 &&
            (
            (k == 1 && terk > terkp1) ||
                (k == 2 && terkm1 > terk > terkp1) ||
                (k > 2 && terkm2 > terkm1 > terk > terkp1)
        )
        k += 1
        terk = terkp1
    else
        while !(terkm2 > terkm1 > terk > terkp1) && k > 2
            terkp1 = terk
            terk = terkm1
            terkm1 = terkm2
            calc_finite_difference_weights!(
                fd_weights, ts_tmp, t + dt, k - 2
            )
            @.. broadcast = false terk_tmp = fd_weights[k - 2, 1] * u
            for i in 2:(k - 2)
                @.. broadcast = false terk_tmp += fd_weights[i, k - 2] * u_history[i - 1]
            end
            @.. broadcast = false terk_tmp *= abs(dt^(k - 2))
            calculate_residuals!(
                atmp, terk_tmp, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            terkm2 = integrator.opts.internalnorm(atmp, t)
            k -= 1
        end
    end
    return k, terk
end

function choose_order!(
        alg::FBDF, integrator,
        cache::FBDFConstantCache,
        ::Val{max_order}
    ) where {max_order}
    (; t, dt, u, uprev) = integrator
    (; ts_tmp, terkm2, terkm1, terk, terkp1, u_history, fd_weights) = cache
    k = cache.order
    if k < max_order && cache.qwait == 0 &&
            (
            (k == 1 && terk > terkp1) ||
                (k == 2 && terkm1 > terk > terkp1) ||
                (k > 2 && terkm2 > terkm1 > terk > terkp1)
        )
        k += 1
        terk = terkp1
    else
        while !(terkm2 > terkm1 > terk > terkp1) && k > 2
            terkp1 = terk
            terk = terkm1
            terkm1 = terkm2
            calc_finite_difference_weights!(
                fd_weights, ts_tmp, t + dt, k - 2
            )
            local terk_tmp
            if u isa Number
                terk_tmp = fd_weights[k - 2, 1] * u
                for i in 2:(k - 2)
                    terk_tmp += fd_weights[i, k - 2] * u_history[i - 1]
                end
                terk_tmp *= abs(dt^(k - 2))
            else
                # we need terk_tmp to be mutable.
                # so it can be updated
                terk_tmp = fd_weights[k - 2, 1] * u
                for i in 2:(k - 2)
                    terk_tmp = @.. terk_tmp + fd_weights[i, k - 2] * u_history[i - 1]
                end
                terk_tmp = @.. terk_tmp * abs(dt^(k - 2))
            end
            atmp = calculate_residuals(
                terk_tmp, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            terkm2 = integrator.opts.internalnorm(atmp, t)
            k -= 1
        end
    end
    return k, terk
end

function stepsize_controller!(
        integrator,
        alg::FBDF
    )
    return stepsize_controller!(integrator, integrator.cache, alg)
end

function stepsize_controller!(
        integrator,
        cache::Union{FBDFCache, FBDFConstantCache},
        alg::FBDF{max_order}
    ) where {
        max_order,
    }
    cache.prev_order = cache.order

    # CVODE-style Stability Limit Detection (STALD)
    # Collect data and check BEFORE order selection, using the step's order and norms.
    # BDF orders 3-5 are only alpha-stable, so eigenvalues near the imaginary axis
    # can cause instability. STALD detects this and forces order reduction.
    step_order = cache.prev_order
    stald_reduce = false
    if step_order >= 3
        stald_collect_data!(
            cache.stald, step_order,
            cache.terkm2, cache.terkm1, cache.terk
        )
        stald_reduce = stald_check!(cache.stald, step_order)
    end

    k, terk = choose_order!(alg, integrator, cache, Val(max_order))
    if k != cache.order
        cache.nconsteps = 0
        cache.order = k
    end

    if stald_reduce
        # Stability violation detected at the step's order: constrain new order
        k = min(k, step_order - 1)
        cache.order = k
        cache.nconsteps = 0
        terk = cache.terkm1
    end

    if iszero(terk)
        q = inv(get_current_qmax(integrator, get_qmax(integrator)))
    else
        # CVODE-style step size formula: eta = 1 / (BIAS2 * dsm)^(1/(k+1))
        # where dsm = terk / (alpha0 * (k+1)) and alpha0 is the BDF leading coefficient.
        # FBDF uses fixed leading coefficients, so alpha0 = bdf_coeffs[k, 1].
        # BIAS2 = 6 matches CVODE (cvode_impl.h).
        alpha0 = cache.bdf_coeffs[k, 1]
        q = ((6 * terk / (alpha0 * (k + 1)))^(1 / (k + 1)))
    end
    return q
end

function step_accept_controller!(integrator, alg::FBDF, q)
    return step_accept_controller!(integrator, integrator.cache, alg, q)
end

function step_accept_controller!(
        integrator, cache::Union{FBDFCache, FBDFConstantCache}, alg::FBDF{max_order},
        q
    ) where {max_order}
    cache.consfailcnt = 0
    if q <= get_qsteady_max(integrator) && q >= get_qsteady_min(integrator)
        q = one(q)
    end
    cache.nconsteps += 1
    cache.iters_from_event += 1
    # CVODE-style qwait countdown for order change gating
    if cache.order != cache.prev_order
        cache.qwait = cache.order + 2 # reset after order change, matching nconsteps >= order + 2
    elseif cache.qwait > 0
        cache.qwait -= 1 # countdown
    end
    return integrator.dt / q
end

function step_reject_controller!(integrator, alg::DFBDF)
    return step_reject_controller!(integrator, integrator.cache, alg)
end
function step_reject_controller!(integrator, cache::Union{DFBDFCache, DFBDFConstantCache}, ::DFBDF)
    return bdf_step_reject_controller!(integrator, cache, cache.terkm1)
end

function post_newton_controller!(integrator, alg::DFBDF)
    return post_newton_controller!(integrator, integrator.cache, alg)
end
function post_newton_controller!(integrator, cache::Union{DFBDFCache, DFBDFConstantCache}, alg::DFBDF)
    if cache.order > 1 && cache.nlsolver.nfails >= 3
        cache.order -= 1
    end
    integrator.dt = integrator.dt / get_failfactor(integrator)
    cache.consfailcnt += 1
    cache.nconsteps = 0
    return nothing
end

function choose_order!(
        alg::DFBDF, integrator,
        cache::DFBDFCache,
        ::Val{max_order}
    ) where {max_order}
    (; t, dt, u, uprev) = integrator
    (; atmp, ts_tmp, terkm2, terkm1, terk, terkp1, terk_tmp, u_history, fd_weights) = cache
    k = cache.order
    # Use CVODE-style qwait countdown: only consider order increase when qwait reaches 0
    if k < max_order && cache.qwait == 0 &&
            (
            (k == 1 && terk > terkp1) ||
                (k == 2 && terkm1 > terk > terkp1) ||
                (k > 2 && terkm2 > terkm1 > terk > terkp1)
        )
        k += 1
        terk = terkp1
    else
        while !(terkm2 > terkm1 > terk > terkp1) && k > 2
            terkp1 = terk
            terk = terkm1
            terkm1 = terkm2
            calc_finite_difference_weights!(
                fd_weights, ts_tmp, t + dt, k - 2
            )
            @.. broadcast = false terk_tmp = fd_weights[k - 2, 1] * u
            for i in 2:(k - 2)
                @.. broadcast = false terk_tmp += fd_weights[i, k - 2] * u_history[i - 1]
            end
            @.. broadcast = false terk_tmp *= abs(dt^(k - 2))
            calculate_residuals!(
                atmp, terk_tmp, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            terkm2 = integrator.opts.internalnorm(atmp, t)
            k -= 1
        end
    end
    return k, terk
end

function choose_order!(
        alg::DFBDF, integrator,
        cache::DFBDFConstantCache,
        ::Val{max_order}
    ) where {max_order}
    (; t, dt, u, uprev) = integrator
    (; ts_tmp, terkm2, terkm1, terk, terkp1, u_history, fd_weights) = cache
    k = cache.order
    if k < max_order && cache.qwait == 0 &&
            (
            (k == 1 && terk > terkp1) ||
                (k == 2 && terkm1 > terk > terkp1) ||
                (k > 2 && terkm2 > terkm1 > terk > terkp1)
        )
        k += 1
        terk = terkp1
    else
        while !(terkm2 > terkm1 > terk > terkp1) && k > 2
            terkp1 = terk
            terk = terkm1
            terkm1 = terkm2
            calc_finite_difference_weights!(
                fd_weights, ts_tmp, t + dt, k - 2
            )
            terk_tmp = @.. broadcast = false fd_weights[k - 2, 1] * u
            if u isa Number
                for i in 2:(k - 2)
                    terk_tmp += fd_weights[i, k - 2] * u_history[i - 1]
                end
                terk_tmp *= abs(dt^(k - 2))
            else
                for i in 2:(k - 2)
                    terk_tmp = @.. terk_tmp + fd_weights[i, k - 2] * u_history[i - 1]
                end
                terk_tmp = @.. broadcast = false terk_tmp * abs(dt^(k - 2))
            end
            atmp = calculate_residuals(
                terk_tmp, uprev, u,
                integrator.opts.abstol, integrator.opts.reltol,
                integrator.opts.internalnorm, t
            )
            terkm2 = integrator.opts.internalnorm(atmp, t)
            k -= 1
        end
    end
    return k, terk
end

function stepsize_controller!(
        integrator,
        alg::DFBDF
    )
    return stepsize_controller!(integrator, integrator.cache, alg)
end

function stepsize_controller!(
        integrator,
        cache::Union{DFBDFCache, DFBDFConstantCache},
        alg::DFBDF{max_order}
    ) where {
        max_order,
    }
    cache.prev_order = cache.order
    k, terk = choose_order!(alg, integrator, cache, Val(max_order))
    if k != cache.order
        cache.nconsteps = 0
        cache.order = k
    end
    if iszero(terk)
        q = inv(get_current_qmax(integrator, get_qmax(integrator)))
    else
        # CVODE-style step size formula matching FBDF change
        alpha0 = cache.bdf_coeffs[k, 1]
        q = ((6 * terk / (alpha0 * (k + 1)))^(1 / (k + 1)))
    end
    return q
end

function step_accept_controller!(integrator, alg::DFBDF, q)
    return step_accept_controller!(integrator, integrator.cache, alg, q)
end

function step_accept_controller!(
        integrator, cache::Union{DFBDFCache, DFBDFConstantCache}, alg::DFBDF{max_order},
        q
    ) where {max_order}
    cache.consfailcnt = 0
    if q <= get_qsteady_max(integrator) && q >= get_qsteady_min(integrator)
        q = one(q)
    end
    cache.nconsteps += 1
    cache.iters_from_event += 1
    # CVODE-style qwait countdown for order change gating
    if cache.order != cache.prev_order
        cache.qwait = cache.order + 2 # reset after order change, matching nconsteps >= order + 2
    elseif cache.qwait > 0
        cache.qwait -= 1 # countdown
    end
    return integrator.dt / q
end

# ============================================================================
# MOOSE234 Controllers
#
# MOOSE234 error estimates (set in perform_step!):
#   cache.terkm1 = |Est2| — order 2 LTE norm (from BDF3-Stab filter: y² − y³)
#   cache.terkp1 = |Est3| — order 3 LTE norm (from FBDF4 filter: y⁴ − y³)
#   cache.terk   = LTE norm at the current active order (used for EEst)
#
# Est4 (order 4 LTE) would require an extra f-eval and is not computed;
# order-4 candidate step uses a geometric-extrapolation heuristic.
# ============================================================================

# --------------------------------------------------------------------------
# Order selection: compute candidate stepsizes for each available order
# (2, 3, 4) and pick the order allowing the largest next step.
#
#   hₖ = dt · (γ / |Estₖ|)^(1/(k+1)),   γ = 0.9
#
# Returns (kₙ, terk_new) where terk_new is the error norm at the chosen order.
# --------------------------------------------------------------------------
function choose_order_moose!(
        integrator,
        cache::Union{MOOSE234Cache, MOOSE234ConstantCache}
    )
    k = cache.order
    dt = integrator.dt
    n_hist = min(cache.iters_from_event + 1, 5)
    max_avail = n_hist >= 4 ? 4 : (n_hist >= 3 ? 3 : 2)

    est2 = cache.terkm1
    est3 = cache.terkp1

    # Candidate step for order 2:  h₂ = dt · (γ / est2)^(1/3)
    if est2 > zero(est2)
        h2 = dt * (0.9 / est2)^(1 / 3)
    else
        h2 = 10 * dt
    end

    # Candidate step for order 3:  h₃ = dt · (γ / est3)^(1/4)
    if max_avail >= 3 && est3 > zero(est3)
        h3 = dt * (0.9 / est3)^(1 / 4)
    else
        h3 = zero(dt)
    end

    # Candidate step for order 4: Est4 is not available, so approximate it.
    # When errors decrease geometrically (est2 > est3), extrapolate
    # est4 ≈ est3² / est2.
    h4 = zero(dt)
    if max_avail >= 4 && est2 > zero(est2) && est3 > zero(est3) && est3 < est2
        est4_approx = est3 * est3 / est2
        h4 = dt * (0.9 / est4_approx)^(1 / 5)
    end

    # Pick the order with the largest candidate step
    kₙ = 2
    hₙ = h2
    if max_avail >= 3 && h3 > hₙ
        kₙ = 3
        hₙ = h3
    end
    if max_avail >= 4 && h4 > hₙ
        kₙ = 4
        hₙ = h4
    end

    # Only allow order *increase* when the qwait countdown has reached 0
    if kₙ > k && cache.qwait > 0
        kₙ = k
    end

    # Error norm at the chosen order (drives the stepsize formula)
    if kₙ == 2
        terk_new = est2
    elseif kₙ == 3
        terk_new = est3
    else
        terk_new = (est2 > zero(est2) && est3 > zero(est3) && est3 < est2) ?
                   est3 * est3 / est2 : est3
    end

    return kₙ, terk_new
end

# --------------------------------------------------------------------------
# Stepsize controller: run order selection, compute q = dt / dt_new
# --------------------------------------------------------------------------
function stepsize_controller!(integrator, alg::MOOSE234)
    return stepsize_controller!(integrator, integrator.cache, alg)
end

function stepsize_controller!(
        integrator,
        cache::Union{MOOSE234Cache, MOOSE234ConstantCache},
        alg::MOOSE234
    )
    cache.prev_order = cache.order

    kₙ, terk = choose_order_moose!(integrator, cache)
    if kₙ != cache.order
        cache.nconsteps = 0
        cache.order = kₙ
    end

    if iszero(terk)
        q = inv(get_current_qmax(integrator, get_qmax(integrator)))
    else
        q = (terk / 0.9)^(1 / (kₙ + 1))
    end
    integrator.qold = q
    return q
end

# --------------------------------------------------------------------------
# Step accept controller
# --------------------------------------------------------------------------
function step_accept_controller!(integrator, alg::MOOSE234, q)
    return step_accept_controller!(integrator, integrator.cache, alg, q)
end

function step_accept_controller!(
        integrator,
        cache::Union{MOOSE234Cache, MOOSE234ConstantCache},
        alg::MOOSE234, q
    )
    cache.consfailcnt = 0
    if q <= get_qsteady_max(integrator) && q >= get_qsteady_min(integrator)
        q = one(q)
    end
    cache.nconsteps += 1
    cache.iters_from_event += 1
    # CVODE-style qwait countdown: order + 2 steps after an order change
    if cache.order != cache.prev_order
        cache.qwait = cache.order + 2
    elseif cache.qwait > 0
        cache.qwait -= 1
    end
    return integrator.dt / q
end

# --------------------------------------------------------------------------
# Step reject controller (MOOSE234-specific, min order = 2)
#
# On repeated failures: drop order, halve step size.
# On consfailcnt > 3 at min order: restart (clear history).
# --------------------------------------------------------------------------
function step_reject_controller!(integrator, ::MOOSE234)
    k = integrator.cache.order
    h = integrator.dt
    integrator.cache.consfailcnt += 1
    integrator.cache.nconsteps = 0

    if integrator.cache.consfailcnt > 1
        h = h / 2
    end

    # Candidate step at current order
    expo = 1 / (k + 1)
    z = 1.2 * integrator.EEst^expo
    hₖ = z <= 10 ? h / z : 0.1 * h

    hₙ = hₖ
    kₙ = k

    if k > 2
        # Error at one order below:
        #   k = 3 → order 2 error (terkm1)
        #   k = 4 → order 3 error (terkp1)
        est_lower = k == 4 ? integrator.cache.terkp1 : integrator.cache.terkm1
        expo_lower = 1 / k  # 1 / ((k-1) + 1)
        zₖ₋₁ = 1.3 * est_lower^expo_lower
        hₖ₋₁ = zₖ₋₁ <= 10 ? h / zₖ₋₁ : 0.1 * h

        if integrator.cache.consfailcnt > 2 || hₖ₋₁ > hₖ
            hₙ = min(h, hₖ₋₁)
            kₙ = k - 1
        end
    end

    # Restart from order 2 (clear history) when stuck at minimum order
    if kₙ == 2 && integrator.cache.consfailcnt > 3
        u_modified!(integrator, true)
    end

    integrator.dt = hₙ
    return integrator.cache.order = kₙ
end

# --------------------------------------------------------------------------
# Post-Newton failure controller
# --------------------------------------------------------------------------
function post_newton_controller!(integrator, alg::MOOSE234)
    (; cache) = integrator
    if cache.order > 2 && cache.nlsolver.nfails >= 3
        cache.order -= 1
    end
    integrator.dt = integrator.dt / get_failfactor(integrator)
    integrator.cache.consfailcnt += 1
    integrator.cache.nconsteps = 0
    return nothing
end
