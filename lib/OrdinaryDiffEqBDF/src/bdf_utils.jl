@inline function U!(k, U)
    @inbounds for r in 1:k
        U[1, r] = -r
        for j in 2:k
            U[j, r] = U[j - 1, r] * ((j - 1) - r) / j
        end
    end
    return nothing
end

function R!(k, ρ, cache)
    (; R) = cache
    @inbounds for r in 1:k
        R[1, r] = -r * ρ
        for j in 2:k
            R[j, r] = R[j - 1, r] * ((j - 1) - r * ρ) / j
        end
    end
    return nothing
end

# This functions takes help of D2 array to create backward differences array D
# Ith row of D2 keeps Ith order backward differences (∇ⁱyₙ)

function backward_diff!(cache::OrdinaryDiffEqMutableCache, D, D2, k, flag = true)
    flag && copyto!(D[1], D2[1, 1])
    for i in 2:k
        for j in 1:(k - i + 1)
            @.. broadcast = false D2[i, j] = D2[i - 1, j] - D2[i - 1, j + 1]
        end
        flag && copyto!(D[i], D2[i, 1])
    end
    return
end

function backward_diff!(cache::OrdinaryDiffEqConstantCache, D, D2, k, flag = true)
    flag && (D[1] = D2[1, 1])
    for i in 2:k
        for j in 1:(k - i + 1)
            D2[i, j] = D2[i - 1, j] - D2[i - 1, j + 1]
        end
        flag && (D[i] = D2[i, 1])
    end
    return
end

# this function updates backward difference array D when stepsize gets change
# Formula -> D = D * (R * U)
# and it is taken from the paper -
# Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s paper
# E. Alberdi Celaya, J. J. Anza Aguirrezabala, and P. Chatzipantelidis
function reinterpolate_history!(cache::OrdinaryDiffEqMutableCache, D, R, k)
    (; tmp) = cache.nlsolver
    fill!(tmp, zero(eltype(D[1])))
    for j in 1:k
        for k in 1:k
            @. tmp += D[k] * R[k, j]
        end
        D[j] .= tmp
        fill!(tmp, zero(eltype(tmp)))
    end
    return
end

function reinterpolate_history!(cache::OrdinaryDiffEqConstantCache, D, R, k)
    tmp = zero(D[1])
    for j in 1:k
        for k in 1:k
            tmp += D[k] * R[k, j]
        end
        D[j] = tmp
    end
    return
end

function _make_bdf_coeffs_fbdf()
    Rat = Rational{Int64}
    return Rat[
        1 -1 0 0 0 0;
        3 // 2 -2 1 // 2 0 0 0;
        11 // 6 -3 3 // 2 -1 // 3 0 0;
        25 // 12 -4 3 -4 // 3 1 // 4 0;
        137 // 60 -5 5 -10 // 3 5 // 4 -1 // 5
    ]
end

function _make_bdf_coeffs_dfbdf()
    Rat = Rational{Int64}
    return Rat[
        1 -1 0 0 0 0;
        2 // 3 -4 // 3 1 // 3 0 0 0;
        6 // 11 -18 // 11 9 // 11 -2 // 11 0 0;
        12 // 25 -48 // 25 36 // 25 -16 // 25 3 // 25 0;
        60 // 137 -300 // 137 300 // 137 -200 // 137 75 // 137 -12 // 137
    ]
end

function calc_R!(R, ρ, k)
    fill!(R, zero(eltype(R)))
    @inbounds for r in 1:k
        R[1, r] = -r * ρ
        for j in 2:k
            R[j, r] = R[j - 1, r] * ((j - 1) - r * ρ) / j
        end
    end
    return R
end

function calc_R(ρ, k, ::Val{N}) where {N}
    R = zeros(typeof(ρ), N, N)
    return calc_R!(R, ρ, k)
end

function update_D!(D, dd, k)
    if dd isa AbstractArray && ArrayInterface.ismutable(dd)
        @.. broadcast = false D[k + 2] = dd - D[k + 1]
        @.. broadcast = false D[k + 1] = dd
        for i in k:-1:1
            @.. broadcast = false D[i] = D[i] + D[i + 1]
        end
    else
        D[k + 2] = dd - D[k + 1]
        D[k + 1] = dd
        for i in k:-1:1
            D[i] = D[i] + D[i + 1]
        end
    end
    return nothing
end

const γₖ = ntuple(k -> sum(Int64(1) // j for j in 1:k), 6)

function error_constant(integrator, alg::QNDF, k)
    (; γₖ) = integrator.cache
    κ = alg.kappa[k]
    return κ * γₖ[k] + inv(k + 1)
end

#This code refers to https://epubs.siam.org/doi/abs/10.1137/S0036144596322507
#Compute all derivatives through k of the polynomials of k+1 points

# In-place version: writes into pre-allocated matrix c
function calc_finite_difference_weights!(c, ts, t, order)
    fill!(c, zero(eltype(c)))
    c1 = one(t)
    c4 = ts[1] - t
    c[1, 1] = one(t)
    for i in 2:(order + 1)
        c2 = one(t)
        c5 = c4
        c4 = ts[i] - t
        @inbounds for j in 1:(i - 1)
            c3 = ts[i] - ts[j]
            c2 *= c3
            if j == i - 1
                for k in i:-1:2
                    c[i, k] = c1 * ((k - 1) * c[i - 1, k - 1] - c5 * c[i - 1, k]) / c2
                end
                c[i, 1] = zero(t)
            end
            for k in i:-1:2
                c[j, k] = (c4 * c[j, k] - (k - 1) * c[j, k - 1]) / c3
            end
            c[j, 1] = zero(t)
        end
        c1 = c2
    end
    return c
end

function calc_finite_difference_weights(ts, t, order, ::Val{N}) where {N}
    max_order = N
    c = zeros(eltype(ts), max_order + 1, max_order + 1)
    calc_finite_difference_weights!(c, ts, t, order)
    return c
end

function reinitFBDF!(integrator, cache)
    # This function is used to initialize arrays that store past history information.
    # It will be used in the first-time step advancing and event handling.
    (;
        consfailcnt, ts, u_history, u_corrector, iters_from_event,
        order,
    ) = cache
    (; t, dt, uprev) = integrator

    if integrator.derivative_discontinuity
        order = cache.order = 1
        consfailcnt = cache.consfailcnt = cache.nconsteps = 0
        cache.qwait = 3 # order + 2, matching nconsteps >= order + 2 for failure-free runs
        iters_from_event = cache.iters_from_event = 0
        if hasproperty(cache, :stald)
            stald_reset!(cache.stald)
        end

        fill!(ts, zero(eltype(ts)))
        for h in u_history
            if h isa AbstractArray && ArrayInterface.ismutable(h)
                fill!(h, zero(eltype(h)))
            end
        end
        for h in u_corrector
            if h isa AbstractArray && ArrayInterface.ismutable(h)
                fill!(h, zero(eltype(h)))
            end
        end
    end

    if uprev isa AbstractArray && ArrayInterface.ismutable(uprev)
        if iters_from_event == 0
            ts[1] = t
            copyto!(u_history[1], uprev)
        elseif iters_from_event == 1 && t != ts[1]
            ts[2] = ts[1]
            ts[1] = t
            copyto!(u_history[2], u_history[1])
            copyto!(u_history[1], uprev)
        elseif consfailcnt == 0
            for i in (order + 2):-1:2
                ts[i] = ts[i - 1]
                copyto!(u_history[i], u_history[i - 1])
            end
            ts[1] = t
            copyto!(u_history[1], uprev)
        end
    else
        if iters_from_event == 0
            ts[1] = t
            u_history[1] = uprev
        elseif iters_from_event == 1 && t != ts[1]
            ts[2] = ts[1]
            ts[1] = t
            u_history[2] = u_history[1]
            u_history[1] = uprev
        elseif consfailcnt == 0
            for i in (order + 2):-1:2
                ts[i] = ts[i - 1]
                u_history[i] = u_history[i - 1]
            end
            ts[1] = t
            u_history[1] = uprev
        end
    end
    return nothing
end

####################################################################
# Chebyshev reference nodes and barycentric weights for FBDF/DFBDF
# dense output. Lagrange interpolation is resampled at these fixed
# nodes during the step (calck block), then evaluated at arbitrary
# Θ during post-solve interpolation using the barycentric formula.
#
# Chebyshev nodes of the first kind on [0,1]:
#   xi_j = (1 + cos((2j-1)π/(2n))) / 2,  j = 1,...,n
#
# Barycentric weights (type 2):
#   w_j = (-1)^(j-1) * sin((2j-1)π/(2n))
#
# These are compile-time constants; no GPU scalar indexing.
####################################################################

const _CHEB_NODES = ntuple(
    n -> ntuple(j -> (1 + cospi((2j - 1) / (2n))) / 2, n), 6
)

const _BARY_WEIGHTS = ntuple(
    n -> ntuple(j -> (-1)^(j - 1) * sinpi((2j - 1) / (2n)), n), 6
)

# Evaluate Lagrange basis L_j(xi) for *actual* (variable) theta nodes.
# All arguments are scalars. Used only in calck/addsteps to resample
# the original interpolant at fixed Chebyshev reference nodes.
@inline function _lagrange_basis_scalar(xi, j, thetas, n)
    theta_j = thetas[j]
    L = one(xi)
    for m in 1:n
        m == j && continue
        L *= (xi - thetas[m]) / (theta_j - thetas[m])
    end
    return L
end

# Evaluate Lagrange interpolant through (thetas[j], values[j]) at scalar xi.
# thetas are scalars, values can be arrays. Out-of-place (allocating).
function _eval_lagrange_oop(xi, thetas, values, n)
    L1 = _lagrange_basis_scalar(xi, 1, thetas, n)
    out = values[1] isa Number ? L1 * values[1] : @.. L1 * values[1]
    for j in 2:n
        Lj = _lagrange_basis_scalar(xi, j, thetas, n)
        if values[1] isa Number
            out += Lj * values[j]
        else
            out = @.. out + Lj * values[j]
        end
    end
    return out
end

# In-place variant: overwrites `out`.
function _eval_lagrange_iip!(out, xi, thetas, values, n)
    L1 = _lagrange_basis_scalar(xi, 1, thetas, n)
    @.. broadcast = false out = L1 * values[1]
    for j in 2:n
        Lj = _lagrange_basis_scalar(xi, j, thetas, n)
        @.. broadcast = false out = out + Lj * values[j]
    end
    return out
end

# Resample: evaluate the Lagrange interpolant defined by
# (thetas[j], values[j]) at fixed Chebyshev reference nodes,
# writing k[i] = p(xi_i). Out-of-place.
function _resample_at_chebyshev!(k, values, thetas, n)
    nodes = _CHEB_NODES[n]
    for i in 1:n
        xi = nodes[i]
        L1 = _lagrange_basis_scalar(xi, 1, thetas, n)
        if values[1] isa Number
            k[i] = L1 * values[1]
            for j in 2:n
                Lj = _lagrange_basis_scalar(xi, j, thetas, n)
                k[i] += Lj * values[j]
            end
        else
            k[i] = @.. L1 * values[1]
            for j in 2:n
                Lj = _lagrange_basis_scalar(xi, j, thetas, n)
                k[i] = @.. k[i] + Lj * values[j]
            end
        end
    end
    return
end

# In-place variant: reads from k[1..n] (the original values),
# copies to scratch[1..n], then writes resampled values back to k[1..n].
function _resample_at_chebyshev_iip!(k, thetas, n, scratch)
    for j in 1:n
        copyto!(scratch[j], k[j])
    end
    nodes = _CHEB_NODES[n]
    for i in 1:n
        xi = nodes[i]
        L1 = _lagrange_basis_scalar(xi, 1, thetas, n)
        @.. broadcast = false k[i] = L1 * scratch[1]
        for j in 2:n
            Lj = _lagrange_basis_scalar(xi, j, thetas, n)
            @.. broadcast = false k[i] = k[i] + Lj * scratch[j]
        end
    end
    return
end

# Direct in-place resampling: writes Chebyshev-resampled values to k[1..n]
# using u (step endpoint) and u_history[1..n-1] (past solution values) as
# sources, without requiring a scratch buffer.  This avoids type mismatches
# during ForwardDiff AD where k has Dual element types but pre-allocated
# scratch buffers are Float64.
function _resample_at_chebyshev_direct_iip!(k, u, u_history, thetas, n)
    nodes = _CHEB_NODES[n]
    for i in 1:n
        xi = nodes[i]
        L1 = _lagrange_basis_scalar(xi, 1, thetas, n)
        @.. broadcast = false k[i] = L1 * u
        for j in 2:n
            Lj = _lagrange_basis_scalar(xi, j, thetas, n)
            @.. broadcast = false k[i] = k[i] + Lj * u_history[j - 1]
        end
    end
    return
end

####################################################################
# MOOSE234 divided-difference coefficient utilities
# DeCaria et al., arXiv:1810.06670v1, Algorithm 7.1 (BACKDIFF)
# and BDFANDFILTCOEFF.
####################################################################

"""
    backdiff(ts::AbstractVector{T}) where {T}

Compute the divided-difference coefficient table for time points `ts` given in
ascending order (`ts[1] < ts[2] < … < ts[n]`).

Implements Algorithm 7.1 (BACKDIFF) from DeCaria et al. (arXiv:1810.06670v1).

Returns an `n × n` matrix `c` where `c[q+1, i]` is the coefficient of the `i`-th
solution value in the `q`-th order divided difference:

    δ^q y = Σᵢ c[q+1, i] · yᵢ

Row `c[1, :]` is the zeroth divided difference (selects the newest value `y[n]`).
"""
function backdiff(ts::AbstractVector{T}) where {T}
    n = length(ts)
    m = n - 1

    # D[j, :] initialised to standard basis vector e_{n+1-j}
    # D[1,:] selects the newest point, D[n,:] selects the oldest.
    D = zeros(T, n, n)
    @inbounds for j in 1:n
        D[j, n + 1 - j] = one(T)
    end

    c = zeros(T, n, n)
    @inbounds for i in 1:n
        c[1, i] = D[1, i]
    end

    @inbounds for q in 1:m
        for j in 1:(n - q)
            denom = ts[n + 1 - j] - ts[n + 1 - j - q]
            for i in 1:n
                D[j, i] = (D[j, i] - D[j + 1, i]) / denom
            end
        end
        for i in 1:n
            c[q + 1, i] = D[1, i]
        end
    end

    return c
end

"""
    bdf_and_filt_coeff(ts::AbstractVector{T}, p::Integer) where {T}

Compute variable-stepsize BDF_p coefficients and FBDF_{p+1} filter weight η^{p+1}.

Implements BDFANDFILTCOEFF from DeCaria et al. (arXiv:1810.06670v1).
Requires `length(ts) ≥ p + 2` (for MOOSE234 with `p = 3`, pass 5 ascending time
points `[tₙ₋₃, tₙ₋₂, tₙ₋₁, tₙ, tₙ₊₁]`).

# Returns
- `α_bar::Vector{T}` — BDF coefficient vector of length `n`. `α_bar[i]` multiplies
  the solution at `ts[i]`. The leading coefficient `α_bar[n]` corresponds to the
  unknown; for the nonlinear solve, `γ = 1 / α_bar[n]`.
- `η::T` — Filter weight η^{p+1} for the FBDF_{p+1} time filter.
- `c::Matrix{T}` — Divided-difference table from [`backdiff`](@ref).
"""
function bdf_and_filt_coeff(ts::AbstractVector{T}, p::Integer) where {T}
    c = backdiff(ts)
    n = length(ts)
    tm = ts[n]   # tₙ₊ₘ (newest time point)

    # η^{p+1} = ∏ᵢ₌₁ᵖ (tm − ts[n−i])  /  Σⱼ₌₁ᵖ⁺¹ (tm − ts[n−j])⁻¹
    num_η = one(T)
    @inbounds for i in 1:p
        num_η *= tm - ts[n - i]
    end
    den_η = zero(T)
    @inbounds for j in 1:(p + 1)
        den_η += one(T) / (tm - ts[n - j])
    end
    η = num_η / den_η

    # ᾱₖ = Σⱼ₌₁ᵖ [∏ᵢ₌₁ʲ⁻¹ (tm − ts[n−i])] · c[j+1, k]
    # for k = (n − p) : n  (only indices that BDF_p touches; earlier ones stay 0)
    α_bar = zeros(T, n)
    @inbounds for k in (n - p):n
        for j in 1:p
            prefactor = one(T)
            for i in 1:(j - 1)
                prefactor *= tm - ts[n - i]
            end
            α_bar[k] += prefactor * c[j + 1, k]
        end
    end

    return α_bar, η, c
end

"""
    bdf3stab_coeff(ts::AbstractVector, µ = 9/125)
    bdf3stab_coeff(c::AbstractMatrix, n::Integer, µ = 9/125)

Compute the BDF3-Stab filter weight for variable stepsizes (equation 3.16 of
DeCaria et al., arXiv:1810.06670v1).

For the filter step

    y² = y³ + weight · δ³y³

returns `weight = µ / c[4, n]`, where `c[4, n]` is the leading coefficient of the
3rd divided difference and `δ³y³ = Σᵢ c[4, i] · yᵢ` (with `yₙ = y³`).

Requires at least 4 time points. The parameter `µ` defaults to `9/125`, the
G-stability-optimal constant from the paper.
"""
function bdf3stab_coeff(ts::AbstractVector, µ = 9 / 125)
    c = backdiff(ts)
    return µ / c[4, length(ts)]
end

function bdf3stab_coeff(c::AbstractMatrix, n::Integer, µ = 9 / 125)
    return µ / c[4, n]
end

function estimate_terk!(integrator, cache, k, ::Val{max_order}) where {max_order}
    #calculate hᵏ⁻¹yᵏ⁻¹
    (; ts_tmp, terk_tmp, u_history, fd_weights) = cache
    (; t, dt, u) = integrator
    calc_finite_difference_weights!(fd_weights, ts_tmp, t + dt, k - 1)
    @.. broadcast = false terk_tmp = fd_weights[1, k] * u
    for i in 2:k
        @.. broadcast = false terk_tmp += fd_weights[i, k] * u_history[i - 1]
    end
    return @.. broadcast = false terk_tmp *= abs(dt^(k - 1))
end

####################################################################
# MOOSE234 history management
# Mirrors reinitFBDF! but with initial order 2 and qwait 4.
####################################################################

function reinitMOOSE!(integrator, cache)
    (;
        consfailcnt, ts, u_history, u_corrector, iters_from_event,
        order,
    ) = cache
    (; t, dt, uprev) = integrator

    if integrator.derivative_discontinuity
        order = cache.order = 2
        consfailcnt = cache.consfailcnt = cache.nconsteps = 0
        cache.qwait = 4 # order + 2
        iters_from_event = cache.iters_from_event = 0

        fill!(ts, zero(eltype(ts)))
        for h in u_history
            if h isa AbstractArray && ArrayInterface.ismutable(h)
                fill!(h, zero(eltype(h)))
            end
        end
        for h in u_corrector
            if h isa AbstractArray && ArrayInterface.ismutable(h)
                fill!(h, zero(eltype(h)))
            end
        end
    end

    if uprev isa AbstractArray && ArrayInterface.ismutable(uprev)
        if iters_from_event == 0
            ts[1] = t
            copyto!(u_history[1], uprev)
        elseif iters_from_event == 1 && t != ts[1]
            ts[2] = ts[1]
            ts[1] = t
            copyto!(u_history[2], u_history[1])
            copyto!(u_history[1], uprev)
        elseif consfailcnt == 0
            for i in (order + 2):-1:2
                ts[i] = ts[i - 1]
                copyto!(u_history[i], u_history[i - 1])
            end
            ts[1] = t
            copyto!(u_history[1], uprev)
        end
    else
        if iters_from_event == 0
            ts[1] = t
            u_history[1] = uprev
        elseif iters_from_event == 1 && t != ts[1]
            ts[2] = ts[1]
            ts[1] = t
            u_history[2] = u_history[1]
            u_history[1] = uprev
        elseif consfailcnt == 0
            for i in (order + 2):-1:2
                ts[i] = ts[i - 1]
                u_history[i] = u_history[i - 1]
            end
            ts[1] = t
            u_history[1] = uprev
        end
    end
    return nothing
end

function estimate_terk(integrator, cache, k, ::Val{max_order}, u) where {max_order}
    (; ts_tmp, u_history, fd_weights) = cache
    (; t, dt) = integrator
    calc_finite_difference_weights!(fd_weights, ts_tmp, t + dt, k - 1)
    terk = @.. broadcast = false fd_weights[1, k] * u
    #@show terk,fd_weights[1,k+1]
    if u isa Number
        for i in 2:k
            terk += fd_weights[i, k] * u_history[i - 1]
        end
        terk *= abs(dt^(k - 1))
    else
        if ArrayInterface.ismutable(terk)
            for i in 2:k
                @.. broadcast = false terk += fd_weights[i, k] * u_history[i - 1]
            end
        else
            for i in 2:k
                terk = @. terk + fd_weights[i, k] * u_history[i - 1]
            end
        end
        terk = @.. broadcast = false terk * abs(dt^(k - 1))
    end
    return terk
end

# Divided-difference filters and residual estimates (DeCaria et al., arXiv:1810.06670).

"""
    backdiff!(c, D, ts, n)

Tabulate divided-difference coefficients for the ascending time points `ts[1:n]`.
On return `c[q + 1, i]` is the coefficient of the solution value at `ts[i]` in the
`q`-th divided difference over the newest `q + 1` points,

    δ^q y = Σᵢ c[q + 1, i] * y(ts[i]),   q = 0, …, n - 1.

`D` is an `n × n` workspace. Both `c` and `D` may be larger than `n × n`; only the
leading block is touched.
"""
function backdiff!(c, D, ts, n)
    T = eltype(c)
    @inbounds for j in 1:n, i in 1:n
        D[j, i] = ifelse(i == n + 1 - j, one(T), zero(T))
    end
    @inbounds for i in 1:n
        c[1, i] = D[1, i]
    end
    @inbounds for q in 1:(n - 1)
        for j in 1:(n - q)
            denom = ts[n + 1 - j] - ts[n + 1 - j - q]
            for i in 1:n
                D[j, i] = (D[j, i] - D[j + 1, i]) / denom
            end
        end
        for i in 1:n
            c[q + 1, i] = D[1, i]
        end
    end
    return c
end

"""
    bdf3stab_coeff(c, n, µ = 9 // 125)

Weight of the BDF3-Stab filter step `y² = y³ + weight * δ³y³` (equation 3.16 of
DeCaria et al., arXiv:1810.06670v1), given a divided-difference table `c` from
[`backdiff!`](@ref) over `n ≥ 4` ascending time points.

`µ = 9/125` is the G-stability-optimal constant of the paper, which makes the
resulting order-2 method A-stable where BDF3 is only α-stable.
"""
bdf3stab_coeff(c, n, µ = 9 // 125) = µ / c[4, n]

# ts_asc[1:n] = [ts[n-1], …, ts[1], tdt], i.e. the cache history in ascending
# order with the point being solved for appended.
function _filt_fill_ts_asc!(ts_asc, ts, tdt, n)
    @inbounds for i in 1:(n - 1)
        ts_asc[i] = ts[n - i]
    end
    @inbounds ts_asc[n] = tdt
    return ts_asc
end

# δ = Σᵢ c[row, i] · y(ts_asc[i]), where `y` is the value at the newest point and
# the older ones come from `u_history` (which is ordered newest-first).
function _filt_divided_diff(c, row, n, y, u_history)
    if y isa Number
        δ = c[row, n] * y
        for i in 1:(n - 1)
            δ += c[row, i] * u_history[n - i]
        end
        return δ
    end
    δ = @.. broadcast = false c[row, n] * y
    for i in 1:(n - 1)
        δ = @.. broadcast = false δ + c[row, i] * u_history[n - i]
    end
    return δ
end

function _filt_divided_diff!(δ, c, row, n, y, u_history)
    @.. broadcast = false δ = c[row, n] * y
    for i in 1:(n - 1)
        @.. broadcast = false δ = δ + c[row, i] * u_history[n - i]
    end
    return δ
end

# Residual of the higher-order BDF formula, used as a defect estimate after
# scaling by its leading coefficient.
function _filt_bdf_residual(α_bar, n, y, u_history, fy)
    if y isa Number
        res = α_bar[n] * y
        for i in 1:(n - 1)
            res += α_bar[i] * u_history[n - i]
        end
        return res - fy
    end
    res = @.. broadcast = false α_bar[n] * y
    for i in 1:(n - 1)
        res = @.. broadcast = false res + α_bar[i] * u_history[n - i]
    end
    return @.. broadcast = false res - fy
end

function _filt_bdf_residual!(res, α_bar, n, y, u_history, fy)
    @.. broadcast = false res = α_bar[n] * y
    for i in 1:(n - 1)
        @.. broadcast = false res = res + α_bar[i] * u_history[n - i]
    end
    @.. broadcast = false res = res - fy
    return res
end

# Cancel the fixed-coefficient solve's error on a monic polynomial of degree
# k+1. Its history interpolant is θ^(k+1) - ∏ⱼ(θ-θⱼ), with θ=(t-tdt)/dt.
function _fbdf_filter_weight(c, ts, n, k, dt, bdf_coeffs)
    defect = zero(eltype(ts))
    for i in 1:k
        θ = -i
        product = one(eltype(ts))
        for j in 1:(n - 1)
            product *= θ - (ts[j] - ts[n]) / dt
        end
        defect -= bdf_coeffs[k, i + 1] * (θ^(k + 1) - product)
    end
    error = defect * dt^(k + 1) / bdf_coeffs[k, 1]
    return error / (1 + c[n, n] * error)
end

function _filt_bdf_coefficients!(α_bar, c, ts, n, p)
    for i in 1:n
        α_bar[i] = zero(eltype(α_bar))
        prefactor = one(eltype(α_bar))
        for j in 1:p
            α_bar[i] += prefactor * c[j + 1, i]
            prefactor *= ts[n] - ts[n - j]
        end
    end
    return α_bar
end

@inline function _filt_step_ratio(est, p)
    iszero(est) && return oftype(est, Inf)
    isfinite(est) && est > 0 || return zero(est)
    return est^(-1 / (p + 1))
end

function _fbdf_select_filter!(integrator, cache, k, max_order, err, err_hi, err_lo)
    order = k
    estimate = err
    ratio = _filt_step_ratio(err, k)
    if k == 3 && _filt_step_ratio(err_lo, 2) > ratio
        order = 2
        estimate = err_lo
        ratio = _filt_step_ratio(err_lo, 2)
    end
    if k < max_order && cache.qwait == 0 && _filt_step_ratio(err_hi, k + 1) > ratio
        order = k + 1
        estimate = err_hi
    end
    cache.filter_order = order
    OrdinaryDiffEqCore.set_EEst!(integrator, estimate)
    return order
end

# NordsieckBDF scales the Newton increment by the test quantity tq[2], which puts
# it in the units of the local error test, so `NLNewton(κ = …)` means CVODE's
# NLSCOEF.
error_constant(integrator, alg::NordsieckBDFAlgs, k) = integrator.cache.tq[2]

function _fbdf_finish_fixed_step!(integrator, cache)
    if cache.time_filter && !integrator.opts.adaptive
        cache.prev_order = cache.order
        cache.order = max(cache.order, cache.filter_order)
        cache.iters_from_event += 1
        cache.nconsteps += 1
        cache.consfailcnt = 0
    end
    return nothing
end
