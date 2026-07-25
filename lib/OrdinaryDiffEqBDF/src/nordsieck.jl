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

# ================================================================= algorithms
@doc generic_solver_docstring(
    """
    An adaptive-order, adaptive-time BDF method on a propagated **Nordsieck**
    history array `zn[j] = h^j/j! * y^(j)(t_n)`, following SUNDIALS CVODE.

    Where `FBDF` stores raw `(t_i, u_i)` history and rebuilds the predictor and the
    error/order estimates from it every step, this method propagates one array:
    predicting is a Pascal-triangle shift, changing the step size is
    `zn[j] *= eta^j`, and accepting a step is a rank-1 update. Because nothing is
    reconstructed from stored points, a loose nonlinear solve cannot be amplified
    into a step-size collapse, so the corrector can be solved to a fraction of the
    local error budget (`nlsolve = NLNewton(κ = …)`, CVODE's NLSCOEF) instead of to
    full accuracy. On stiff benchmarks that is worth roughly a factor of two in
    f-evaluations relative to `FBDF`.

    Dense output is the Nordsieck polynomial itself and is free.
    """,
    "NordsieckBDF",
    "Multistep Method",
    """
    @article{byrne1975polyalgorithm,
    title={A polyalgorithm for the numerical solution of ordinary differential equations},
    author={Byrne, George D and Hindmarsh, Alan C},
    journal={ACM Transactions on Mathematical Software},
    volume={1}, number={1}, pages={71--96}, year={1975}}
    @article{hindmarsh2005sundials,
    title={{SUNDIALS}: Suite of nonlinear and differential/algebraic equation solvers},
    author={Hindmarsh, Alan C and Brown, Peter N and Grant, Keith E and Lee, Steven L
            and Serban, Radu and Shumaker, Dan E and Woodward, Carol S},
    journal={ACM Transactions on Mathematical Software},
    volume={31}, number={3}, pages={363--396}, year={2005}}""",
    """
    - `nlsolve`: nonlinear solver for the implicit stage. Its `κ` acts as CVODE's
      NLSCOEF, i.e. the fraction of the local error budget the corrector is allowed
      to consume, because the increment norm is scaled by the test quantity `tq[2]`.
    - `max_order`: maximum BDF order (1–5).
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`.
    """,
    """
    nlsolve = NLNewton(),
    extrapolant = :linear,
    max_order::Val{MO} = Val{5}(),
    step_limiter! = trivial_limiter!,
    """
)
struct NordsieckBDF{MO, AD, F, F2, T, StepLimiter, CJ, QT} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    tol::T
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
    stald::Bool
    qmax::QT
    qsteady_min::QT
    qsteady_max::QT
end

function NordsieckBDF(;
        max_order::Val{MO} = Val{5}(),
        autodiff = AutoForwardDiff(), concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(), tol = nothing,
        extrapolant = :linear, step_limiter! = trivial_limiter!, stald = false,
        qsteady_min = 1 // 1, qsteady_max = 1 // 1, qmax = 10 // 1
    ) where {MO}
    autodiff = _fixup_ad(autodiff)
    return NordsieckBDF(
        max_order, linsolve, nlsolve, tol, extrapolant, step_limiter!,
        autodiff, _unwrap_val(concrete_jac), stald, qmax, qsteady_min, qsteady_max
    )
end

@truncate_stacktrace NordsieckBDF

@doc generic_solver_docstring(
    """
    Fully implicit DAE solver: the [`NordsieckBDF`](@ref) method applied to
    `f(du, u, p, t) = 0`. The Nordsieck array supplies both the state predictor and
    the derivative `du = zn[1]/h`, so the corrector solves
    `f((zn[1] + l[1]*acor)/h, ypred + acor, p, t) = 0` with leading coefficient
    `cj = l[1]/h` — the same role IDA's `cj` plays.
    """,
    "DNordsieckBDF",
    "Fully Implicit Multistep Method",
    """
    @article{hindmarsh2005sundials,
    title={{SUNDIALS}: Suite of nonlinear and differential/algebraic equation solvers},
    author={Hindmarsh, Alan C and Brown, Peter N and Grant, Keith E and Lee, Steven L
            and Serban, Radu and Shumaker, Dan E and Woodward, Carol S},
    journal={ACM Transactions on Mathematical Software},
    volume={31}, number={3}, pages={363--396}, year={2005}}""",
    """
    - `nlsolve`: nonlinear solver for the implicit stage; its `κ` acts as NLSCOEF.
    - `max_order`: maximum BDF order (1–5).
    """,
    """
    nlsolve = NLNewton(),
    extrapolant = :linear,
    max_order::Val{MO} = Val{5}(),
    """
)
struct DNordsieckBDF{MO, AD, F, F2, T, CJ, QT} <: DAEAlgorithm
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    tol::T
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
    stald::Bool
    qmax::QT
    qsteady_min::QT
    qsteady_max::QT
end

function DNordsieckBDF(;
        max_order::Val{MO} = Val{5}(),
        autodiff = AutoForwardDiff(), concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(), tol = nothing,
        extrapolant = :linear, stald = false,
        qsteady_min = 1 // 1, qsteady_max = 1 // 1, qmax = 10 // 1
    ) where {MO}
    autodiff = _fixup_ad(autodiff)
    return DNordsieckBDF(
        max_order, linsolve, nlsolve, tol, extrapolant, autodiff,
        _unwrap_val(concrete_jac), stald, qmax, qsteady_min, qsteady_max
    )
end

@truncate_stacktrace DNordsieckBDF

const NordsieckBDFAlgs = Union{NordsieckBDF, DNordsieckBDF}

alg_order(alg::NordsieckBDF) = 1  # dummy: the running order lives in the cache
alg_order(alg::DNordsieckBDF) = 1
isadaptive(alg::DNordsieckBDF) = true
get_current_alg_order(alg::NordsieckBDFAlgs, cache) = cache.order
get_current_adaptive_order(alg::NordsieckBDFAlgs, cache) = cache.order
has_stiff_interpolation(::NordsieckBDFAlgs) = true

# The Newton increment norm is scaled by tq[2], which converts it into the units of
# the local error test. That makes `NLNewton(κ = ...)` mean CVODE's NLSCOEF: the
# fraction of the error-test budget the corrector may consume.
has_special_newton_error(alg::NordsieckBDFAlgs) = true
error_constant(integrator, alg::NordsieckBDFAlgs, k) = integrator.cache.tq[2]


# The step-size logic is CVODE's (`cvSetEta` keeps h unless eta >= 1.5), so the
# generic qsteady band must not also clamp it.
qmax_default(::NordsieckBDFAlgs) = 10 // 1
qsteady_min_default(::NordsieckBDFAlgs) = 1 // 1
qsteady_max_default(::NordsieckBDFAlgs) = 1 // 1
gamma_default(::NordsieckBDFAlgs) = 1 // 1

function default_controller(QT, alg::NordsieckBDFAlgs)
    return BDFController(
        QT, alg; qmax = alg.qmax,
        qsteady_min = alg.qsteady_min, qsteady_max = alg.qsteady_max
    )
end

function setup_controller_cache(
        alg::NordsieckBDFAlgs, cache, controller::BDFController, ::Type{E}, disco_probs
    ) where {E}
    QT = _resolved_QT(controller.basic)
    basic = resolve_basic(controller.basic, alg, QT; disco_probs)
    return BDFControllerCache{QT, E, typeof(cache), eltype(disco_probs)}(
        BDFController(basic), cache, oneunit(E)
    )
end

# ================================================================= caches
@cache mutable struct NordsieckBDFCache{
        MO, N, rateType, uNoUnitsType, uType, tType, coeffType, staldType, StepLimiter,
    } <: BDFMutableCache
    fsalfirst::rateType
    nlsolver::N
    zn::Vector{uType}
    ypred::uType
    acor::uType
    tempv::uType
    tmp::uType
    atmp::uNoUnitsType
    l::coeffType
    tau::coeffType
    tq::coeffType
    order::Int
    qprime::Int
    qwait::Int
    nst::Int
    nef::Int
    ncf::Int
    indx_acor::Int
    max_order::Val{MO}
    max_order_int::Int
    hscale::tType
    eta::tType
    etamax::tType
    etaq::tType
    etaqm1::tType
    etaqp1::tType
    saved_tq5::tType
    predicted::Bool
    stald::staldType
    step_limiter!::StepLimiter
end

@truncate_stacktrace NordsieckBDFCache 1

mutable struct NordsieckBDFConstantCache{MO, N, uType, tType, coeffType, staldType} <:
    OrdinaryDiffEqConstantCache
    nlsolver::N
    zn::Vector{uType}
    ypred::uType
    acor::uType
    l::coeffType
    tau::coeffType
    tq::coeffType
    order::Int
    qprime::Int
    qwait::Int
    nst::Int
    nef::Int
    ncf::Int
    indx_acor::Int
    max_order::Val{MO}
    max_order_int::Int
    hscale::tType
    eta::tType
    etamax::tType
    etaq::tType
    etaqm1::tType
    etaqp1::tType
    saved_tq5::tType
    predicted::Bool
    stald::staldType
end

# DAE caches carry `u₀` because `get_dae_uprev` uses it as the predictor the
# correction `z` is measured against.
@cache mutable struct DNordsieckBDFCache{
        MO, N, rateType, uNoUnitsType, uType, tType, coeffType, staldType,
    } <: BDFMutableCache
    fsalfirst::rateType
    nlsolver::N
    zn::Vector{uType}
    u₀::uType
    ypred::uType
    acor::uType
    tempv::uType
    tmp::uType
    atmp::uNoUnitsType
    l::coeffType
    tau::coeffType
    tq::coeffType
    order::Int
    qprime::Int
    qwait::Int
    nst::Int
    nef::Int
    ncf::Int
    indx_acor::Int
    max_order::Val{MO}
    max_order_int::Int
    hscale::tType
    eta::tType
    etamax::tType
    etaq::tType
    etaqm1::tType
    etaqp1::tType
    saved_tq5::tType
    predicted::Bool
    stald::staldType
end

@truncate_stacktrace DNordsieckBDFCache 1

mutable struct DNordsieckBDFConstantCache{MO, N, uType, tType, coeffType, staldType} <:
    OrdinaryDiffEqConstantCache
    nlsolver::N
    zn::Vector{uType}
    u₀::uType
    ypred::uType
    acor::uType
    l::coeffType
    tau::coeffType
    tq::coeffType
    order::Int
    qprime::Int
    qwait::Int
    nst::Int
    nef::Int
    ncf::Int
    indx_acor::Int
    max_order::Val{MO}
    max_order_int::Int
    hscale::tType
    eta::tType
    etamax::tType
    etaq::tType
    etaqm1::tType
    etaqp1::tType
    saved_tq5::tType
    predicted::Bool
    stald::staldType
end

const NordsieckCaches = Union{
    NordsieckBDFCache, NordsieckBDFConstantCache,
    DNordsieckBDFCache, DNordsieckBDFConstantCache,
}

function alg_cache(
        alg::NordsieckBDF{MO}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck, ::Val{true}, verbose
    ) where {MO, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = one(tTypeNoUnits), one(tTypeNoUnits)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    zn = [zero(u) for _ in 1:(MO + 1)]
    coeffs() = zeros(typeof(t), MO + 3)
    tq = zeros(typeof(t), 6)
    stald = StabilityLimitDetectionState(real(uBottomEltypeNoUnits); enabled = alg.stald)
    return NordsieckBDFCache{
        MO, typeof(nlsolver), typeof(rate_prototype),
        typeof(similar(u, uEltypeNoUnits)), typeof(u), typeof(t),
        typeof(coeffs()), typeof(stald), typeof(alg.step_limiter!),
    }(
        zero(rate_prototype), nlsolver, zn, zero(u), zero(u), zero(u), zero(u),
        similar(u, uEltypeNoUnits), coeffs(), coeffs(), tq,
        1, 1, 2, 0, 0, 0, MO, Val(MO), MO,
        zero(t), one(t), typeof(t)(NORD_ETA_MAX_FS), one(t), zero(t), zero(t),
        zero(t), false, stald, alg.step_limiter!
    )
end

function alg_cache(
        alg::NordsieckBDF{MO}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck, ::Val{false}, verbose
    ) where {MO, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = one(tTypeNoUnits), one(tTypeNoUnits)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    zn = [zero(u) for _ in 1:(MO + 1)]
    coeffs() = zeros(typeof(t), MO + 3)
    stald = StabilityLimitDetectionState(real(uBottomEltypeNoUnits); enabled = alg.stald)
    return NordsieckBDFConstantCache{
        MO, typeof(nlsolver), typeof(u), typeof(t), typeof(coeffs()), typeof(stald),
    }(
        nlsolver, zn, zero(u), zero(u), coeffs(), coeffs(), zeros(typeof(t), 6),
        1, 1, 2, 0, 0, 0, MO, Val(MO), MO,
        zero(t), one(t), typeof(t)(NORD_ETA_MAX_FS), one(t), zero(t), zero(t),
        zero(t), false, stald
    )
end

function alg_cache(
        alg::DNordsieckBDF{MO}, du, u, res_prototype, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck, ::Val{true}, verbose
    ) where {MO, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = one(tTypeNoUnits), one(tTypeNoUnits)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    zn = [zero(u) for _ in 1:(MO + 1)]
    coeffs() = zeros(typeof(t), MO + 3)
    stald = StabilityLimitDetectionState(real(uBottomEltypeNoUnits); enabled = alg.stald)
    return DNordsieckBDFCache{
        MO, typeof(nlsolver), typeof(rate_prototype),
        typeof(similar(u, uEltypeNoUnits)), typeof(u), typeof(t), typeof(coeffs()),
        typeof(stald),
    }(
        zero(rate_prototype), nlsolver, zn, zero(u), zero(u), zero(u), zero(u),
        zero(u), similar(u, uEltypeNoUnits), coeffs(), coeffs(), zeros(typeof(t), 6),
        1, 1, 2, 0, 0, 0, MO, Val(MO), MO,
        zero(t), one(t), typeof(t)(NORD_ETA_MAX_FS), one(t), zero(t), zero(t),
        zero(t), false, stald
    )
end

function alg_cache(
        alg::DNordsieckBDF{MO}, du, u, res_prototype, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck, ::Val{false}, verbose
    ) where {MO, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = one(tTypeNoUnits), one(tTypeNoUnits)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    zn = [zero(u) for _ in 1:(MO + 1)]
    coeffs() = zeros(typeof(t), MO + 3)
    stald = StabilityLimitDetectionState(real(uBottomEltypeNoUnits); enabled = alg.stald)
    return DNordsieckBDFConstantCache{
        MO, typeof(nlsolver), typeof(u), typeof(t), typeof(coeffs()), typeof(stald),
    }(
        nlsolver, zn, zero(u), zero(u), zero(u), coeffs(), coeffs(),
        zeros(typeof(t), 6), 1, 1, 2, 0, 0, 0, MO, Val(MO), MO,
        zero(t), one(t), typeof(t)(NORD_ETA_MAX_FS), one(t), zero(t), zero(t),
        zero(t), false, stald
    )
end

# With `adaptive = false` the integrator never calls `step_accept_controller!`, so
# the array has to be committed here instead. The order still ramps up to
# `max_order` through the usual `qwait` countdown, which keeps the saved
# correction column valid for the order increase.
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

# ================================================================= perform_step!
function initialize!(integrator, cache::NordsieckBDFCache)
    integrator.kshortsize = cache.max_order_int + 1
    resize!(integrator.k, integrator.kshortsize)
    @inbounds for i in 1:(integrator.kshortsize)
        integrator.k[i] = zero(integrator.u)
    end
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return nothing
end

function initialize!(integrator, cache::NordsieckBDFConstantCache)
    integrator.kshortsize = cache.max_order_int + 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    @inbounds for i in 1:(integrator.kshortsize)
        integrator.k[i] = zero(integrator.u)
    end
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return nothing
end

# Store the Nordsieck columns for dense output. Unused columns stay zero, so the
# interpolant can sum all of them without knowing the order at evaluation time.
@inline function _nordsieck_store_k!(integrator, cache, ::Val{true})
    @inbounds for j in 0:(cache.max_order_int)
        if j <= cache.order
            copyto!(integrator.k[j + 1], cache.zn[j + 1])
        else
            fill!(integrator.k[j + 1], zero(eltype(integrator.u)))
        end
    end
    return nothing
end
@inline function _nordsieck_store_k!(integrator, cache, ::Val{false})
    @inbounds for j in 0:(cache.max_order_int)
        integrator.k[j + 1] = j <= cache.order ? cache.zn[j + 1] : zero(integrator.u)
    end
    return nothing
end

function perform_step!(integrator, cache::NordsieckBDFCache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    iip = Val(true)
    nlsolver = cache.nlsolver
    nordsieck_needs_start(integrator, cache) && nordsieck_start!(integrator, cache, iip)

    nordsieck_prepare!(integrator, cache, iip)
    nordsieck_predict!(cache, iip)
    nordsieck_set_coeffs!(cache, dt)

    zn = cache.zn
    copyto!(cache.ypred, zn[1])
    l1 = cache.l[2]

    # (l1/dt) M u - f(u) = (l1*ypred - zn[1])/dt  =>  gamma = 1/l1, alpha = 1
    mass_matrix = f.mass_matrix
    @.. broadcast = false cache.tmp = (l1 * cache.ypred - zn[2]) / dt
    if mass_matrix === I
        copyto!(nlsolver.tmp, cache.tmp)
    else
        mul!(nlsolver.tmp, mass_matrix, cache.tmp)
    end
    markfirststage!(nlsolver)
    copyto!(nlsolver.z, cache.ypred)
    nlsolver.γ = inv(l1)
    nlsolver.α = one(l1)
    nlsolver.method = COEFFICIENT_MULTISTEP
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    @.. broadcast = false u = z
    cache.step_limiter!(u, integrator, p, t + dt)

    @.. broadcast = false cache.acor = u - cache.ypred
    if integrator.opts.adaptive
        OrdinaryDiffEqCore.set_EEst!(
            integrator,
            cache.tq[2] * _nord_wrms(integrator, cache, cache.acor, uprev, u)
        )
    end

    # zn[1] is h*f(u) exactly once the corrector has converged, so fsallast is free
    @.. broadcast = false integrator.fsallast = (zn[2] + l1 * cache.acor) / dt
    _nordsieck_finish_fixed!(integrator, cache, iip)
    return nothing
end

function perform_step!(integrator, cache::NordsieckBDFConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    iip = Val(false)
    nlsolver = cache.nlsolver
    nordsieck_needs_start(integrator, cache) && nordsieck_start!(integrator, cache, iip)

    nordsieck_prepare!(integrator, cache, iip)
    nordsieck_predict!(cache, iip)
    nordsieck_set_coeffs!(cache, dt)

    zn = cache.zn
    cache.ypred = zn[1]
    l1 = cache.l[2]

    mass_matrix = f.mass_matrix
    tmp = @.. (l1 * cache.ypred - zn[2]) / dt
    nlsolver.tmp = mass_matrix === I ? tmp : mass_matrix * tmp
    markfirststage!(nlsolver)
    nlsolver.z = cache.ypred
    nlsolver.γ = inv(l1)
    nlsolver.α = one(l1)
    nlsolver.method = COEFFICIENT_MULTISTEP
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return
    u = z

    cache.acor = @.. u - cache.ypred
    if integrator.opts.adaptive
        OrdinaryDiffEqCore.set_EEst!(
            integrator,
            cache.tq[2] * _nord_wrms(integrator, cache, cache.acor, uprev, u)
        )
    end
    integrator.fsallast = @.. (zn[2] + l1 * cache.acor) / dt
    integrator.u = u
    _nordsieck_finish_fixed!(integrator, cache, iip)
    return nothing
end

"""
    nordsieck_stald!(cache, integrator, u, uprev, dsm)

CVODE `cvBDFStab`: feed the stability-limit detector the scaled derivative norms

    sqm2 = (q-1)! * ‖zn[q-1]‖,  sqm1 = (q-1)! * q * ‖zn[q]‖,
    sq   = (q-1)! * q * (q+1) * acnrm / tq[5]

and return whether the order must be reduced. BDF orders 3–5 are only α-stable, so
eigenvalues near the imaginary axis can destabilise the integration at high order;
this detects that from the history and drops the order. Off by default, as in CVODE.
"""
function nordsieck_stald!(cache, integrator, u, uprev, dsm)
    cache.stald.enabled || return false
    q = cache.order
    q < 3 && return false
    T = typeof(cache.eta)
    fact = one(T)
    for i in 1:(q - 1)
        fact *= i
    end
    acnrm = dsm / max(cache.tq[2], eps(T))
    sq = fact * q * (q + 1) * acnrm / max(cache.tq[5], eps(T))
    sqm1 = fact * q * _nord_wrms(integrator, cache, cache.zn[q + 1], uprev, u)
    sqm2 = fact * _nord_wrms(integrator, cache, cache.zn[q], uprev, u)
    stald_collect_data!(cache.stald, q, sqm2, sqm1, sq)
    return stald_check!(cache.stald, q)
end

# ================================================================= controllers
# The controller hooks own the Nordsieck bookkeeping: `perform_step!` leaves the
# array in the *predicted* state, accepting commits it with the rank-1 update, and
# rejecting undoes the Pascal shift. `cache.predicted` makes both idempotent, so
# the hooks are safe regardless of the order the integrator calls them in.

_nordsieck_iip(::Union{NordsieckBDFCache, DNordsieckBDFCache}) = Val(true)
_nordsieck_iip(::Union{NordsieckBDFConstantCache, DNordsieckBDFConstantCache}) = Val(false)

stepsize_controller!(integrator, alg::NordsieckBDFAlgs) = nothing

function step_accept_controller!(integrator, alg::NordsieckBDFAlgs, q)
    cache = integrator.cache
    iip = _nordsieck_iip(cache)
    (; dt, u, uprev) = integrator
    T = typeof(cache.eta)
    dsm = OrdinaryDiffEqCore.get_EEst(integrator)
    acor = cache.acor

    # STALD inspects the step that was just taken, so it runs before the array is
    # advanced and before the new order is chosen.
    stald_reduce = nordsieck_stald!(cache, integrator, u, uprev, dsm)

    nordsieck_complete!(cache, dt, acor, iip)
    cache.nef = 0
    cache.ncf = 0
    # dense output data: the committed Nordsieck columns about t_{n+1}
    _nordsieck_store_k!(integrator, cache, iip)

    if cache.etamax == one(T)
        # a failure earlier in this step forbids growth (CVODE cvPrepareNextStep)
        cache.qwait = max(cache.qwait, 2)
        cache.qprime = cache.order
        cache.eta = one(T)
    else
        cache.etaq = inv((NORD_BIAS2 * dsm)^(one(T) / (cache.order + 1)) + NORD_ADDON)
        if cache.qwait != 0
            cache.eta = cache.etaq
            cache.qprime = cache.order
            nordsieck_set_eta!(cache, integrator)
        else
            cache.qwait = 2
            cache.etaqm1 = nordsieck_compute_etaqm1(cache, integrator, u, uprev)
            cache.etaqp1 = nordsieck_compute_etaqp1(cache, integrator, u, uprev, acor, dt)
            nordsieck_choose_eta!(cache, integrator, u, uprev, acor, dt, iip)
            nordsieck_set_eta!(cache, integrator)
        end
    end
    if stald_reduce && cache.qprime > 1
        # a stability violation overrides whatever order the error estimates chose
        cache.qprime = min(cache.qprime, cache.order - 1)
        cache.eta = min(cache.eta, one(T))
    end
    # after the first step the growth cap drops from ETA_MAX_FS to the steady value
    cache.etamax = T(NORD_ETA_MAX_GS)
    eta = min(cache.eta, get_current_qmax(integrator, get_qmax(integrator)))
    return dt * eta
end

function step_reject_controller!(integrator, alg::NordsieckBDFAlgs)
    cache = integrator.cache
    iip = _nordsieck_iip(cache)
    T = typeof(cache.eta)
    dsm = OrdinaryDiffEqCore.get_EEst(integrator)
    nordsieck_restore!(cache, iip)
    cache.nef += 1
    cache.etamax = one(T)

    if cache.nef <= NORD_MXNEF1
        eta = inv((NORD_BIAS2 * dsm)^(one(T) / (cache.order + 1)) + NORD_ADDON)
        eta = max(T(NORD_ETA_MIN_EF), eta)
        if cache.nef >= NORD_SMALL_NEF
            eta = min(eta, T(NORD_ETA_MAX_EF))
        end
        cache.eta = eta
    elseif cache.order > 1
        # after repeated failures drop the order and retry
        cache.eta = T(NORD_ETA_MIN_EF)
        nordsieck_adjust_order!(cache, -1, iip)
        cache.order -= 1
        cache.qprime = cache.order
        cache.qwait = cache.order + 1
    else
        # already at order 1: restart the history from scratch
        cache.eta = T(NORD_ETA_MIN_EF)
        cache.qwait = NORD_LONG_WAIT
        derivative_discontinuity!(integrator, true)
    end
    nordsieck_rescale!(cache, cache.eta, iip)
    integrator.dt = cache.hscale
    return integrator.dt
end

function post_newton_controller!(integrator, alg::NordsieckBDFAlgs)
    cache = integrator.cache
    iip = _nordsieck_iip(cache)
    T = typeof(cache.eta)
    nordsieck_restore!(cache, iip)
    cache.etamax = one(T)
    cache.ncf += 1
    if cache.ncf >= 3 && cache.order > 1
        # repeated corrector failures usually mean the high-order predictor is the
        # problem, so drop the order as well as the step size
        nordsieck_adjust_order!(cache, -1, iip)
        cache.order -= 1
        cache.qprime = cache.order
        cache.qwait = cache.order + 1
        cache.ncf = 0
    end
    cache.eta = T(NORD_ETA_CF)
    nordsieck_rescale!(cache, cache.eta, iip)
    integrator.dt = cache.hscale
    return nothing
end

# ================================================================= DAE stepping
# The corrector solves  f((zn[1] + l1*acor)/dt, ypred + acor, p, t+dt) = 0.
# `_compute_rhs!` for a `DAEFunction` forms `du = (tmp + α*z)*invγdt` and
# `u = cache.u₀ + z`, so with α = 1, γ = 1/l1 (invγdt = l1/dt) it is enough to set
# `tmp = zn[1]/l1` and `u₀ = ypred`; then `cj = l1/dt`, exactly IDA's leading
# coefficient.
function initialize!(integrator, cache::DNordsieckBDFCache)
    integrator.kshortsize = cache.max_order_int + 1
    resize!(integrator.k, integrator.kshortsize)
    @inbounds for i in 1:(integrator.kshortsize)
        integrator.k[i] = zero(integrator.u)
    end
    copyto!(integrator.fsalfirst, integrator.du)
    return nothing
end

function initialize!(integrator, cache::DNordsieckBDFConstantCache)
    integrator.kshortsize = cache.max_order_int + 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    @inbounds for i in 1:(integrator.kshortsize)
        integrator.k[i] = zero(integrator.u)
    end
    integrator.fsalfirst = integrator.du
    return nothing
end

# For a DAE the derivative history comes from `integrator.du`, not from f.
function nordsieck_start_dae!(integrator, cache, iip::Val{true})
    (; uprev, dt) = integrator
    zn = cache.zn
    copyto!(zn[1], uprev)
    @.. broadcast = false zn[2] = dt * integrator.du
    @inbounds for j in 3:length(zn)
        fill!(zn[j], zero(eltype(uprev)))
    end
    return _nordsieck_start_common!(cache, dt)
end
function nordsieck_start_dae!(integrator, cache, iip::Val{false})
    (; uprev, dt) = integrator
    zn = cache.zn
    zn[1] = uprev
    zn[2] = dt * integrator.du
    @inbounds for j in 3:length(zn)
        zn[j] = zero(uprev)
    end
    return _nordsieck_start_common!(cache, dt)
end

function perform_step!(integrator, cache::DNordsieckBDFCache, repeat_step = false)
    (; t, dt, uprev, u, p) = integrator
    iip = Val(true)
    nlsolver = cache.nlsolver
    nordsieck_needs_start(integrator, cache) &&
        nordsieck_start_dae!(integrator, cache, iip)

    nordsieck_prepare!(integrator, cache, iip)
    nordsieck_predict!(cache, iip)
    nordsieck_set_coeffs!(cache, dt)

    zn = cache.zn
    copyto!(cache.ypred, zn[1])
    copyto!(cache.u₀, cache.ypred)
    l1 = cache.l[2]

    @.. broadcast = false nlsolver.tmp = zn[2] / l1
    markfirststage!(nlsolver)
    fill!(nlsolver.z, zero(eltype(nlsolver.z)))
    nlsolver.γ = inv(l1)
    nlsolver.α = one(l1)
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    @.. broadcast = false cache.acor = z
    @.. broadcast = false u = cache.ypred + z
    if integrator.opts.adaptive
        OrdinaryDiffEqCore.set_EEst!(
            integrator,
            cache.tq[2] * _nord_wrms(integrator, cache, cache.acor, uprev, u)
        )
    end
    @.. broadcast = false integrator.fsallast = (zn[2] + l1 * cache.acor) / dt
    _nordsieck_finish_fixed!(integrator, cache, iip)
    return nothing
end

function perform_step!(integrator, cache::DNordsieckBDFConstantCache, repeat_step = false)
    (; t, dt, uprev, p) = integrator
    iip = Val(false)
    nlsolver = cache.nlsolver
    nordsieck_needs_start(integrator, cache) &&
        nordsieck_start_dae!(integrator, cache, iip)

    nordsieck_prepare!(integrator, cache, iip)
    nordsieck_predict!(cache, iip)
    nordsieck_set_coeffs!(cache, dt)

    zn = cache.zn
    cache.ypred = zn[1]
    cache.u₀ = cache.ypred
    l1 = cache.l[2]

    nlsolver.tmp = @.. zn[2] / l1
    markfirststage!(nlsolver)
    nlsolver.z = zero(cache.ypred)
    nlsolver.γ = inv(l1)
    nlsolver.α = one(l1)
    z = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    cache.acor = z
    u = @.. cache.ypred + z
    if integrator.opts.adaptive
        OrdinaryDiffEqCore.set_EEst!(
            integrator,
            cache.tq[2] * _nord_wrms(integrator, cache, cache.acor, uprev, u)
        )
    end
    integrator.fsallast = @.. (zn[2] + l1 * cache.acor) / dt
    integrator.u = u
    _nordsieck_finish_fixed!(integrator, cache, iip)
    return nothing
end

# ================================================================= interpolation
# Dense output is the Nordsieck polynomial itself:
#     y(t_n + s*h) = sum_j zn[j] * s^j,   s = Θ - 1
# `integrator.k` holds the columns (unused ones zeroed), so the order at
# evaluation time is not needed.
const NORDSIECK_CACHES = Union{
    NordsieckBDFCache, NordsieckBDFConstantCache,
    DNordsieckBDFCache, DNordsieckBDFConstantCache,
}

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::NORDSIECK_CACHES,
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    s = Θ - one(Θ)
    out = k[length(k)]
    @inbounds for j in (length(k) - 1):-1:1
        out = @.. out * s + k[j]
    end
    return out
end

@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k, cache::NORDSIECK_CACHES,
        idxs, T::Type{Val{0}}, differential_vars
    )
    s = Θ - one(Θ)
    out = k[length(k)][idxs]
    @inbounds for j in (length(k) - 1):-1:1
        out = @.. out * s + k[j][idxs]
    end
    return out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::NORDSIECK_CACHES,
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    s = Θ - one(Θ)
    copyto!(out, k[length(k)])
    @inbounds for j in (length(k) - 1):-1:1
        kj = k[j]
        @.. broadcast = false out = out * s + kj
    end
    return out
end

@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k, cache::NORDSIECK_CACHES,
        idxs, T::Type{Val{0}}, differential_vars
    )
    s = Θ - one(Θ)
    @views copyto!(out, k[length(k)][idxs])
    @inbounds for j in (length(k) - 1):-1:1
        kj = @view k[j][idxs]
        @.. broadcast = false out = out * s + kj
    end
    return out
end

# derivatives: differentiate the same polynomial, d^m/dt^m = (1/dt^m) d^m/ds^m
for (TV, m) in ((Val{1}, 1), (Val{2}, 2), (Val{3}, 3))
    @eval @muladd function _ode_interpolant(
            Θ, dt, y₀, y₁, k, cache::NORDSIECK_CACHES,
            idxs::Nothing, T::Type{$TV}, differential_vars
        )
        s = Θ - one(Θ)
        n = length(k)
        out = zero(k[1])
        @inbounds for j in n:-1:($m + 1)
            c = one(Θ)
            for r in 0:($m - 1)
                c *= (j - 1 - r)
            end
            out = @.. out * s + c * k[j]
        end
        return @.. out / dt^$m
    end

    @eval @muladd function _ode_interpolant!(
            out, Θ, dt, y₀, y₁, k, cache::NORDSIECK_CACHES,
            idxs::Nothing, T::Type{$TV}, differential_vars
        )
        s = Θ - one(Θ)
        n = length(k)
        fill!(out, zero(eltype(out)))
        @inbounds for j in n:-1:($m + 1)
            c = one(Θ)
            for r in 0:($m - 1)
                c *= (j - 1 - r)
            end
            kj = k[j]
            @.. broadcast = false out = out * s + c * kj
        end
        @.. broadcast = false out = out / dt^$m
        return out
    end
end

# `addsteps!` is a no-op: the interpolation data is written during the step.
function _ode_addsteps!(
        k, t, uprev, u, dt, f, p, cache::NORDSIECK_CACHES,
        always_calc_begin = false, allow_calc_end = true, force_calc_end = false
    )
    return nothing
end
