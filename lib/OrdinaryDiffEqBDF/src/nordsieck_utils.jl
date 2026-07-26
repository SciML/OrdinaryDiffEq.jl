# Variable-order variable-step BDF on a propagated Nordsieck history array,
#
#     zn[j] = h^j / j! * y^(j)(t_n),   j = 0…q
#
# following SUNDIALS CVODE (`cvode.c`, `cvode_nls.c`). Unlike `FBDF`, which stores
# raw `(t_i, u_i)` history and rebuilds the predictor and the error/order estimates
# from it each step, every operation here is a cheap transform of one propagated
# array: the predictor is a Pascal-triangle shift, a step-size change is
# `zn[j] *= eta^j`, and the corrector update is rank-1 `zn[j] += l[j] * acor`.
#
# Nothing is reconstructed from stored points, so a nonlinear-solve residual is
# never amplified by an extrapolating predictor. That is what lets the corrector
# be solved loosely (NLSCOEF-style) without the step size collapsing.

# ---- CVODE constants (cvode_impl.h) ----
const NORD_ADDON = 1 // 1000000
const NORD_BIAS1 = 6
const NORD_BIAS2 = 6
const NORD_BIAS3 = 10
const NORD_ETA_MAX_FX = 3 // 2
const NORD_ETA_MAX_FS = 10000
const NORD_ETA_MAX_GS = 10
const NORD_ETA_MIN = 1 // 10
const NORD_ETA_MAX_EF = 1 // 5
const NORD_ETA_MIN_EF = 1 // 10
const NORD_ETA_CF = 1 // 4
const NORD_SMALL_NST = 10
const NORD_SMALL_NEF = 2
const NORD_MXNEF1 = 2
const NORD_LONG_WAIT = 10

# ---------------------------------------------------------------- coefficients
"""
    nordsieck_set_coeffs!(cache, dt)

CVODE `cvSetBDF` + `cvSetTqBDF`: the Nordsieck correction vector `l` and the test
quantities `tq` for the current order and step-size history.
"""
function nordsieck_set_coeffs!(cache, dt)
    l = cache.l
    tau = cache.tau
    tq = cache.tq
    q = cache.order
    T = eltype(l)
    l[1] = one(T)
    l[2] = one(T)
    xi_inv = one(T)
    xistar_inv = one(T)
    @inbounds for i in 2:q
        l[i + 1] = zero(T)
    end
    alpha0 = -one(T)
    alpha0_hat = -one(T)
    hsum = dt
    if q > 1
        for j in 2:(q - 1)
            hsum += tau[j - 1]
            xi_inv = dt / hsum
            alpha0 -= one(T) / j
            @inbounds for i in j:-1:1
                l[i + 1] += l[i] * xi_inv
            end
        end
        alpha0 -= one(T) / q
        xistar_inv = -l[2] - alpha0
        hsum += tau[q - 1]
        xi_inv = dt / hsum
        alpha0_hat = -l[2] - xi_inv
        @inbounds for i in q:-1:1
            l[i + 1] += l[i] * xistar_inv
        end
    end

    A1 = one(T) - alpha0_hat + alpha0
    A2 = one(T) + q * A1
    tq[2] = abs(A1 / (alpha0 * A2))
    tq[5] = abs(A2 * xistar_inv / (l[q + 1] * xi_inv))
    if cache.qwait == 1
        if q > 1
            C = xistar_inv / l[q + 1]
            A3 = alpha0 + one(T) / q
            A4 = alpha0_hat + xi_inv
            tq[1] = abs(C * (one(T) - A4 + A3) / A3)
        else
            tq[1] = one(T)
        end
        hsum += tau[q]
        xi_inv = dt / hsum
        A5 = alpha0 - one(T) / (q + 1)
        A6 = alpha0_hat - xi_inv
        tq[3] = abs((one(T) - A6 + A5) / A2 / (xi_inv * (q + 2) * A5))
    end
    return nothing
end

# ---------------------------------------------------------------- array ops
# Two methods each: in-place for mutable caches, rebinding for constant caches
# (where `u` may be a scalar or a `StaticArray`).
@inline function _nord_shift!(zn::Vector, q, ::Val{true})
    @inbounds for k in 1:q, j in q:-1:k
        zj = zn[j + 1]
        @.. broadcast = false zn[j] = zn[j] + zj
    end
    return nothing
end
@inline function _nord_shift!(zn::Vector, q, ::Val{false})
    @inbounds for k in 1:q, j in q:-1:k
        zn[j] = zn[j] + zn[j + 1]
    end
    return nothing
end

@inline function _nord_unshift!(zn::Vector, q, ::Val{true})
    @inbounds for k in 1:q, j in q:-1:k
        zj = zn[j + 1]
        @.. broadcast = false zn[j] = zn[j] - zj
    end
    return nothing
end
@inline function _nord_unshift!(zn::Vector, q, ::Val{false})
    @inbounds for k in 1:q, j in q:-1:k
        zn[j] = zn[j] - zn[j + 1]
    end
    return nothing
end

@inline function _nord_scale!(zn::Vector, q, eta, ::Val{true})
    factor = eta
    @inbounds for j in 1:q
        @.. broadcast = false zn[j + 1] = zn[j + 1] * factor
        factor *= eta
    end
    return nothing
end
@inline function _nord_scale!(zn::Vector, q, eta, ::Val{false})
    factor = eta
    @inbounds for j in 1:q
        zn[j + 1] = zn[j + 1] * factor
        factor *= eta
    end
    return nothing
end

@inline function _nord_rank1!(zn::Vector, l, q, acor, ::Val{true})
    @inbounds for j in 0:q
        lj = l[j + 1]
        @.. broadcast = false zn[j + 1] = zn[j + 1] + lj * acor
    end
    return nothing
end
@inline function _nord_rank1!(zn::Vector, l, q, acor, ::Val{false})
    @inbounds for j in 0:q
        zn[j + 1] = zn[j + 1] + l[j + 1] * acor
    end
    return nothing
end

# ---------------------------------------------------------------- history moves
"""
    nordsieck_predict!(cache, iip)

CVODE `cvPredict`: advance the Nordsieck array to the new time by repeated
addition (a Pascal-triangle shift). `zn[0]` becomes the predictor.
"""
function nordsieck_predict!(cache, iip)
    _nord_shift!(cache.zn, cache.order, iip)
    cache.predicted = true
    return nothing
end

"""
    nordsieck_restore!(cache, iip)

CVODE `cvRestore`: undo [`nordsieck_predict!`](@ref) after a rejected step.
"""
function nordsieck_restore!(cache, iip)
    cache.predicted || return nothing
    _nord_unshift!(cache.zn, cache.order, iip)
    cache.predicted = false
    return nothing
end

"""
    nordsieck_rescale!(cache, eta, iip)

CVODE `cvRescale`: a step-size change is exactly `zn[j] *= eta^j`. No history is
re-interpolated, which is why a step-size change costs nothing in accuracy here.
"""
function nordsieck_rescale!(cache, eta, iip)
    _nord_scale!(cache.zn, cache.order, eta, iip)
    cache.hscale = cache.hscale * eta
    return nothing
end

"""
    nordsieck_complete!(cache, dt, acor, iip)

CVODE `cvCompleteStep`: commit an accepted step with the rank-1 update
`zn[j] += l[j] * acor` and shift the step-size history `tau`.
"""
function nordsieck_complete!(cache, dt, acor, iip)
    q = cache.order
    tau = cache.tau
    @inbounds for i in q:-1:2
        tau[i] = tau[i - 1]
    end
    if q == 1 && cache.nst > 0
        tau[2] = tau[1]
    end
    tau[1] = dt
    _nord_rank1!(cache.zn, cache.l, q, acor, iip)
    cache.nst += 1
    cache.predicted = false
    cache.qwait -= 1
    if cache.qwait == 1 && q != cache.max_order_int
        _nord_setacor!(cache, acor, iip)
        cache.saved_tq5 = cache.tq[5]
        cache.indx_acor = cache.max_order_int
    end
    return nothing
end

@inline function _nord_setacor!(cache, acor, ::Val{true})
    @.. broadcast = false cache.zn[cache.max_order_int + 1] = acor
    return nothing
end
@inline function _nord_setacor!(cache, acor, ::Val{false})
    cache.zn[cache.max_order_int + 1] = acor
    return nothing
end

# ---------------------------------------------------------------- order changes
"""
    nordsieck_increase!(cache, iip)

CVODE `cvIncreaseBDF`: extend the array with a new highest column built from the
saved correction.
"""
function nordsieck_increase!(cache, iip)
    l = cache.l
    T = eltype(l)
    mo = cache.max_order_int
    @inbounds for i in 0:mo
        l[i + 1] = zero(T)
    end
    l[3] = one(T)
    alpha1 = one(T)
    prod = one(T)
    xiold = one(T)
    alpha0 = -one(T)
    hsum = cache.hscale
    q = cache.order
    if q > 1
        for j in 1:(q - 1)
            hsum += cache.tau[j + 1]
            xi = hsum / cache.hscale
            prod *= xi
            alpha0 -= one(T) / (j + 1)
            alpha1 += one(T) / xi
            @inbounds for i in (j + 2):-1:2
                l[i + 1] = l[i + 1] * xiold + l[i]
            end
            xiold = xi
        end
    end
    A1 = (-alpha0 - alpha1) / prod
    _nord_newcol!(cache, A1, iip)
    _nord_addcol!(cache, iip)
    return nothing
end

@inline function _nord_newcol!(cache, A1, ::Val{true})
    src = cache.zn[cache.indx_acor + 1]
    @.. broadcast = false cache.zn[cache.order + 2] = A1 * src
    return nothing
end
@inline function _nord_newcol!(cache, A1, ::Val{false})
    cache.zn[cache.order + 2] = A1 * cache.zn[cache.indx_acor + 1]
    return nothing
end

@inline function _nord_addcol!(cache, ::Val{true})
    dst = cache.zn[cache.order + 2]
    @inbounds for j in 2:(cache.order)
        lj = cache.l[j + 1]
        @.. broadcast = false cache.zn[j + 1] = cache.zn[j + 1] + lj * dst
    end
    return nothing
end
@inline function _nord_addcol!(cache, ::Val{false})
    dst = cache.zn[cache.order + 2]
    @inbounds for j in 2:(cache.order)
        cache.zn[j + 1] = cache.zn[j + 1] + cache.l[j + 1] * dst
    end
    return nothing
end

"""
    nordsieck_decrease!(cache, iip)

CVODE `cvDecreaseBDF`: drop the highest column.
"""
function nordsieck_decrease!(cache, iip)
    l = cache.l
    T = eltype(l)
    @inbounds for i in 0:(cache.max_order_int)
        l[i + 1] = zero(T)
    end
    l[3] = one(T)
    hsum = zero(T)
    q = cache.order
    for j in 1:(q - 2)
        hsum += cache.tau[j]
        xi = hsum / cache.hscale
        @inbounds for i in (j + 2):-1:2
            l[i + 1] = l[i + 1] * xi + l[i]
        end
    end
    _nord_subcol!(cache, iip)
    return nothing
end

@inline function _nord_subcol!(cache, ::Val{true})
    zq = cache.zn[cache.order + 1]
    @inbounds for j in 2:(cache.order - 1)
        lj = cache.l[j + 1]
        @.. broadcast = false cache.zn[j + 1] = cache.zn[j + 1] - lj * zq
    end
    return nothing
end
@inline function _nord_subcol!(cache, ::Val{false})
    zq = cache.zn[cache.order + 1]
    @inbounds for j in 2:(cache.order - 1)
        cache.zn[j + 1] = cache.zn[j + 1] - cache.l[j + 1] * zq
    end
    return nothing
end

function nordsieck_adjust_order!(cache, deltaq, iip)
    (cache.order == 2 && deltaq != 1) && return nothing
    deltaq == 1 && return nordsieck_increase!(cache, iip)
    deltaq == -1 && return nordsieck_decrease!(cache, iip)
    return nothing
end

"""
    nordsieck_prepare!(integrator, cache, iip)

CVODE `cvAdjustParams`: apply any pending order change and rescale the array to
the step size the controller selected. Runs at the top of each attempted step.
"""
function nordsieck_prepare!(integrator, cache, iip)
    if cache.qprime != cache.order
        nordsieck_adjust_order!(cache, cache.qprime - cache.order, iip)
        cache.order = cache.qprime
        cache.qwait = cache.order + 1
    end
    dt = integrator.dt
    if dt != cache.hscale && !iszero(cache.hscale)
        nordsieck_rescale!(cache, dt / cache.hscale, iip)
    end
    cache.hscale = dt
    return nothing
end

# ---------------------------------------------------------------- eta selection
function nordsieck_compute_etaqm1(cache, integrator, u, uprev)
    T = typeof(cache.eta)
    cache.order <= 1 && return zero(T)
    ddn = _nord_wrms(integrator, cache, cache.zn[cache.order + 1], uprev, u) * cache.tq[1]
    return inv((NORD_BIAS1 * ddn)^(one(T) / cache.order) + NORD_ADDON)
end

function nordsieck_compute_etaqp1(cache, integrator, u, uprev, acor, dt)
    T = typeof(cache.eta)
    cache.order == cache.max_order_int && return zero(T)
    iszero(cache.saved_tq5) && return zero(T)
    L = cache.order + 1
    cquot = (cache.tq[5] / cache.saved_tq5) * (dt / cache.tau[2])^L
    dup = _nord_wrms_diff(integrator, cache, acor, cquot, u, uprev) * cache.tq[3]
    return inv((NORD_BIAS3 * dup)^(one(T) / (L + 1)) + NORD_ADDON)
end

"""
    nordsieck_choose_eta!(cache, ...)

CVODE `cvChooseEta`: pick whichever of orders `q-1`, `q`, `q+1` permits the
largest step.
"""
function nordsieck_choose_eta!(cache, integrator, u, uprev, acor, dt, iip)
    etam = max(cache.etaqm1, max(cache.etaq, cache.etaqp1))
    if etam < NORD_ETA_MAX_FX
        cache.eta = one(cache.eta)
        cache.qprime = cache.order
        return nothing
    end
    if etam == cache.etaq
        cache.eta = cache.etaq
        cache.qprime = cache.order
    elseif etam == cache.etaqm1
        cache.eta = cache.etaqm1
        cache.qprime = cache.order - 1
    else
        cache.eta = cache.etaqp1
        cache.qprime = cache.order + 1
        # stash Delta_n for the order increase
        _nord_setacor!(cache, acor, iip)
    end
    return nothing
end

function nordsieck_set_eta!(cache, integrator)
    if cache.eta < NORD_ETA_MAX_FX
        # inside the fixed band: keep the step size. CVODE never shrinks after an
        # accepted step, which is what keeps W refactorizations rare.
        cache.eta = one(cache.eta)
    else
        cache.eta = min(cache.eta, cache.etamax)
    end
    return nothing
end

# ---------------------------------------------------------------- weighted norms
@inline function _nord_wrms(integrator, cache::OrdinaryDiffEqMutableCache, x, uprev, u)
    (; abstol, reltol, internalnorm) = integrator.opts
    calculate_residuals!(cache.atmp, x, uprev, u, abstol, reltol, internalnorm, integrator.t)
    return internalnorm(cache.atmp, integrator.t)
end
@inline function _nord_wrms(integrator, cache::OrdinaryDiffEqConstantCache, x, uprev, u)
    (; abstol, reltol, internalnorm) = integrator.opts
    atmp = calculate_residuals(x, uprev, u, abstol, reltol, internalnorm, integrator.t)
    return internalnorm(atmp, integrator.t)
end

# ‖acor - cquot * zn[qmax]‖ without allocating a second temporary
@inline function _nord_wrms_diff(
        integrator, cache::OrdinaryDiffEqMutableCache, acor, cquot, u, uprev
    )
    zq = cache.zn[cache.max_order_int + 1]
    @.. broadcast = false cache.tempv = acor - cquot * zq
    return _nord_wrms(integrator, cache, cache.tempv, uprev, u)
end
@inline function _nord_wrms_diff(
        integrator, cache::OrdinaryDiffEqConstantCache, acor, cquot, u, uprev
    )
    tempv = @.. acor - cquot * cache.zn[cache.max_order_int + 1]
    return _nord_wrms(integrator, cache, tempv, uprev, u)
end

function _nordsieck_finish_fixed!(integrator, cache, iip)
    integrator.opts.adaptive && return nothing
    nordsieck_complete!(cache, integrator.dt, cache.acor, iip)
    _nordsieck_store_k!(integrator, cache, iip)
    if cache.qwait <= 0
        cache.qwait = 2
        if cache.order < cache.max_order_int
            cache.qprime = cache.order + 1
        end
    end
    return nothing
end

# ================================================================= start / restart
"""
    nordsieck_start!(integrator, cache, iip)

Seed the array with the order-1 history `zn = [uprev, dt*f(uprev)]`. Used on the
first step and after a discontinuity, mirroring CVODE's cold start.
"""
function nordsieck_start!(integrator, cache, iip::Val{true})
    (; uprev, dt) = integrator
    zn = cache.zn
    copyto!(zn[1], uprev)
    @.. broadcast = false zn[2] = dt * integrator.fsalfirst
    @inbounds for j in 3:length(zn)
        fill!(zn[j], zero(eltype(uprev)))
    end
    return _nordsieck_start_common!(cache, dt)
end
function nordsieck_start!(integrator, cache, iip::Val{false})
    (; uprev, dt) = integrator
    zn = cache.zn
    zn[1] = uprev
    zn[2] = dt * integrator.fsalfirst
    @inbounds for j in 3:length(zn)
        zn[j] = zero(uprev)
    end
    return _nordsieck_start_common!(cache, dt)
end
function _nordsieck_start_common!(cache, dt)
    cache.order = 1
    cache.qprime = 1
    cache.qwait = 2
    cache.nef = 0
    cache.ncf = 0
    cache.hscale = dt
    cache.eta = one(cache.eta)
    cache.etamax = typeof(cache.etamax)(NORD_ETA_MAX_FS)
    cache.saved_tq5 = zero(cache.saved_tq5)
    cache.indx_acor = cache.max_order_int
    cache.predicted = false
    stald_reset!(cache.stald)
    fill!(cache.tau, zero(eltype(cache.tau)))
    return nothing
end

@inline function nordsieck_needs_start(integrator, cache)
    return cache.nst == 0 || integrator.derivative_discontinuity
end
