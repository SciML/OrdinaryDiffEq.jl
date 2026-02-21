### Type unions for dispatch
QNDF_CACHES = Union{QNDFConstantCache, QNDFCache}
FBDF_CACHES = Union{FBDFConstantCache, FBDFCache, FBDFCacheVF64, DFBDFConstantCache, DFBDFCache}
BDF_CACHES_WITH_INTERPOLATIONS = Union{QNDF_CACHES, FBDF_CACHES}

### Fallbacks to capture unsupported derivative orders
function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::BDF_CACHES_WITH_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::BDF_CACHES_WITH_INTERPOLATIONS,
        idxs, T::Type{Val{D}}, differential_vars
    ) where {D}
    throw(DerivativeOrderNotPossibleError())
end

####################################################################
# QNDF: Newton backward difference interpolation
#
# k stores backward differences D[:,j].
# p(Θ) = y₁ + Σ_{j=1}^{order} φ_j(Θ-1) * k[j]
#
# where φ_j(σ) = σ*(σ+1)*...*(σ+j-1)/j! and σ = Θ-1.
# Recurrence: φ₁(σ) = σ, φ_{j+1}(σ) = φ_j(σ) * (σ+j)/(j+1)
#
# At Θ=0: σ=-1, φ₁(-1)=-1, φ_j(-1)=0 for j≥2, so p(0) = y₁ - k[1] = y₀
# At Θ=1: σ=0, all φ_j(0)=0, so p(1) = y₁
####################################################################

## QNDF Val{0}: Function value interpolation

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::QNDF_CACHES,
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    σ = Θ - 1
    num_k = length(k)
    φ = σ
    out = @.. y₁ + φ * k[1]
    for j in 2:num_k
        φ *= (σ + j - 1) / j
        out = @.. out + φ * k[j]
    end
    return out
end

# Out-of-place, with idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::QNDF_CACHES,
        idxs, T::Type{Val{0}}, differential_vars
    )
    σ = Θ - 1
    num_k = length(k)
    φ = σ
    out = @.. y₁[idxs] + φ * k[1][idxs]
    for j in 2:num_k
        φ *= (σ + j - 1) / j
        out = @.. out + φ * k[j][idxs]
    end
    return out
end

# In-place, no idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::QNDF_CACHES,
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    σ = Θ - 1
    num_k = length(k)
    φ = σ
    @.. out = y₁ + φ * k[1]
    for j in 2:num_k
        φ *= (σ + j - 1) / j
        @.. out = out + φ * k[j]
    end
    return out
end

# In-place, with idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::QNDF_CACHES,
        idxs, T::Type{Val{0}}, differential_vars
    )
    σ = Θ - 1
    num_k = length(k)
    φ = σ
    @views @.. out = y₁[idxs] + φ * k[1][idxs]
    for j in 2:num_k
        φ *= (σ + j - 1) / j
        @views @.. out = out + φ * k[j][idxs]
    end
    return out
end

## QNDF Val{1}: First derivative dp/dt = (1/dt) * Σ dφ'_j(σ) * k[j]
#
# Recurrence for φ'_j:
#   φ'₁(σ) = 1
#   φ'_{j+1}(σ) = (φ'_j(σ) * (σ+j) + φ_j(σ)) / (j+1)

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::QNDF_CACHES,
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    σ = Θ - 1
    num_k = length(k)
    invdt = inv(dt)
    φ = σ
    dφ = one(Θ)
    out = @.. dφ * k[1]
    for j in 2:num_k
        dφ_new = (dφ * (σ + j - 1) + φ) / j
        φ *= (σ + j - 1) / j
        dφ = dφ_new
        out = @.. out + dφ * k[j]
    end
    return @.. out * invdt
end

# Out-of-place, with idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::QNDF_CACHES,
        idxs, T::Type{Val{1}}, differential_vars
    )
    σ = Θ - 1
    num_k = length(k)
    invdt = inv(dt)
    φ = σ
    dφ = one(Θ)
    out = @.. dφ * k[1][idxs]
    for j in 2:num_k
        dφ_new = (dφ * (σ + j - 1) + φ) / j
        φ *= (σ + j - 1) / j
        dφ = dφ_new
        out = @.. out + dφ * k[j][idxs]
    end
    return @.. out * invdt
end

# In-place, no idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::QNDF_CACHES,
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    σ = Θ - 1
    num_k = length(k)
    invdt = inv(dt)
    φ = σ
    dφ = one(Θ)
    @.. out = dφ * k[1]
    for j in 2:num_k
        dφ_new = (dφ * (σ + j - 1) + φ) / j
        φ *= (σ + j - 1) / j
        dφ = dφ_new
        @.. out = out + dφ * k[j]
    end
    @.. out = out * invdt
    return out
end

# In-place, with idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::QNDF_CACHES,
        idxs, T::Type{Val{1}}, differential_vars
    )
    σ = Θ - 1
    num_k = length(k)
    invdt = inv(dt)
    φ = σ
    dφ = one(Θ)
    @views @.. out = dφ * k[1][idxs]
    for j in 2:num_k
        dφ_new = (dφ * (σ + j - 1) + φ) / j
        φ *= (σ + j - 1) / j
        dφ = dφ_new
        @views @.. out = out + dφ * k[j][idxs]
    end
    @.. out = out * invdt
    return out
end

####################################################################
# FBDF / DFBDF: Non-equidistant Lagrange interpolation
#
# All data lives in k. y₀ and y₁ are NOT used.
#
# Layout: k has 2*half entries where half = length(k)÷2.
#   k[1..half] = solution values at past/current times
#   k[half+1..2*half] = corresponding Θ positions
#
# For scalar u, k[half+j] is the Θ value directly.
# For vector u, k[half+j] is filled with the same Θ value; use first().
#
# For dense output (stored in calck block at step end):
#   k[1] = u_new (step end solution), Θ₁ = 1
#   k[1+j] = u_history[:,j],          Θ_{1+j} = (ts[j] - t) / dt
#   for j = 1..order. Total: order+1 data points.
#
# For predictor (rebuilt at step start):
#   k[j] = u_history[:,j],            Θ_j = (ts[j] - t) / dt
#   for j = 1..order+1. Total: order+1 data points.
#
# General Lagrange formula:
#   p(Θ) = Σ_{j=1}^n L_j(Θ) * k[j]
#   L_j(Θ) = Π_{m≠j} (Θ - Θ_m) / (Θ_j - Θ_m)
#
# This matches calc_Lagrange_interp: same polynomial through the same
# actual solution values at their actual times.
####################################################################

# Helper: determine number of active data points from k entries.
# Layout: k[1..half] solution values, k[half+1..2*half] Θ positions.
# Active entries are non-zero (counting from the top).
function _bdf_active_order(k)
    half = length(k) ÷ 2
    n = half
    while n > 0 && iszero(k[n])
        n -= 1
    end
    return max(n, 1)
end

# Extract scalar Θ from a k entry (which may be a scalar or vector filled
# with the same value).
_get_theta(x::Number) = x
_get_theta(x) = first(x)

# Compute Lagrange basis value L_i(Θ) for node i among n nodes.
# Node positions: Θ_j = _get_theta(k[half + j]) for j = 1..n.
@inline function _lagrange_basis(Θ, i, n, k)
    half = length(k) ÷ 2
    θ_i = _get_theta(k[half + i])
    Li = one(Θ)
    for m in 1:n
        m == i && continue
        θ_m = _get_theta(k[half + m])
        Li *= (Θ - θ_m) / (θ_i - θ_m)
    end
    return Li
end

## FBDF Val{0}: Lagrange function value interpolation (non-equidistant nodes)

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    n = _bdf_active_order(k)

    L1 = _lagrange_basis(Θ, 1, n, k)
    out = @.. L1 * k[1]

    for j in 2:n
        Lj = _lagrange_basis(Θ, j, n, k)
        out = @.. out + Lj * k[j]
    end

    return out
end

# Out-of-place, with idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs, T::Type{Val{0}}, differential_vars
    )
    n = _bdf_active_order(k)

    L1 = _lagrange_basis(Θ, 1, n, k)
    out = @.. L1 * k[1][idxs]

    for j in 2:n
        Lj = _lagrange_basis(Θ, j, n, k)
        out = @.. out + Lj * k[j][idxs]
    end

    return out
end

# In-place, no idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    n = _bdf_active_order(k)

    L1 = _lagrange_basis(Θ, 1, n, k)
    @.. out = L1 * k[1]

    for j in 2:n
        Lj = _lagrange_basis(Θ, j, n, k)
        @.. out = out + Lj * k[j]
    end

    return out
end

# In-place, with idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs, T::Type{Val{0}}, differential_vars
    )
    n = _bdf_active_order(k)

    L1 = _lagrange_basis(Θ, 1, n, k)
    @views @.. out = L1 * k[1][idxs]

    for j in 2:n
        Lj = _lagrange_basis(Θ, j, n, k)
        @views @.. out = out + Lj * k[j][idxs]
    end

    return out
end

## FBDF Val{1} and Val{2}: Derivatives of Lagrange interpolant
#
# p'(Θ) = (1/dt) * Σ_{i=1}^n k[i] * L'_i(Θ)
# L'_i(Θ) = L_i(Θ) * Σ_{m≠i} 1/(Θ - Θ_m)
#   (using logarithmic derivative of the Lagrange basis)

# Product-rule form of L'_i(Θ):
#   L'_i(Θ) = Σ_{m≠i} [1/(θ_i - θ_m)] * Π_{l≠i,l≠m} (Θ - θ_l)/(θ_i - θ_l)
# This form avoids 0*Inf = NaN at node points (where Θ = θ_j for some j).
@inline function _lagrange_basis_deriv(Θ, i, n, k)
    half = length(k) ÷ 2
    θ_i = _get_theta(k[half + i])
    result = zero(Θ)
    for m in 1:n
        m == i && continue
        θ_m = _get_theta(k[half + m])
        term = inv(θ_i - θ_m)
        for l in 1:n
            (l == i || l == m) && continue
            θ_l = _get_theta(k[half + l])
            term *= (Θ - θ_l) / (θ_i - θ_l)
        end
        result += term
    end
    return result
end

# Product-rule form of L''_i(Θ):
#   L''_i(Θ) = Σ_{m≠i} Σ_{p≠i,p≠m} [1/((θ_i-θ_m)(θ_i-θ_p))]
#              * Π_{l≠i,l≠m,l≠p} (Θ-θ_l)/(θ_i-θ_l)
@inline function _lagrange_basis_deriv2(Θ, i, n, k)
    half = length(k) ÷ 2
    θ_i = _get_theta(k[half + i])
    result = zero(Θ)
    for m in 1:n
        m == i && continue
        θ_m = _get_theta(k[half + m])
        for p in 1:n
            (p == i || p == m) && continue
            θ_p = _get_theta(k[half + p])
            term = inv(θ_i - θ_m) * inv(θ_i - θ_p)
            for l in 1:n
                (l == i || l == m || l == p) && continue
                θ_l = _get_theta(k[half + l])
                term *= (Θ - θ_l) / (θ_i - θ_l)
            end
            result += term
        end
    end
    return result
end

# Product-rule form of L'''_i(Θ):
#   L'''_i(Θ) = Σ_{m≠i} Σ_{p≠i,p≠m} Σ_{q≠i,q≠m,q≠p}
#              [1/((θ_i-θ_m)(θ_i-θ_p)(θ_i-θ_q))]
#              * Π_{l≠i,l≠m,l≠p,l≠q} (Θ-θ_l)/(θ_i-θ_l)
@inline function _lagrange_basis_deriv3(Θ, i, n, k)
    half = length(k) ÷ 2
    θ_i = _get_theta(k[half + i])
    result = zero(Θ)
    for m in 1:n
        m == i && continue
        θ_m = _get_theta(k[half + m])
        for p in 1:n
            (p == i || p == m) && continue
            θ_p = _get_theta(k[half + p])
            for q in 1:n
                (q == i || q == m || q == p) && continue
                θ_q = _get_theta(k[half + q])
                term = inv(θ_i - θ_m) * inv(θ_i - θ_p) * inv(θ_i - θ_q)
                for l in 1:n
                    (l == i || l == m || l == p || l == q) && continue
                    θ_l = _get_theta(k[half + l])
                    term *= (Θ - θ_l) / (θ_i - θ_l)
                end
                result += term
            end
        end
    end
    return result
end

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt = inv(dt)

    dL1 = _lagrange_basis_deriv(Θ, 1, n, k)
    out = @.. (dL1 * invdt) * k[1]

    for j in 2:n
        dLj = _lagrange_basis_deriv(Θ, j, n, k)
        cj = dLj * invdt
        out = @.. out + cj * k[j]
    end
    return out
end

# Out-of-place, with idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs, T::Type{Val{1}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt = inv(dt)

    dL1 = _lagrange_basis_deriv(Θ, 1, n, k)
    out = @.. (dL1 * invdt) * k[1][idxs]

    for j in 2:n
        dLj = _lagrange_basis_deriv(Θ, j, n, k)
        cj = dLj * invdt
        out = @.. out + cj * k[j][idxs]
    end
    return out
end

# In-place, no idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt = inv(dt)

    dL1 = _lagrange_basis_deriv(Θ, 1, n, k)
    @.. out = (dL1 * invdt) * k[1]

    for j in 2:n
        dLj = _lagrange_basis_deriv(Θ, j, n, k)
        cj = dLj * invdt
        @.. out = out + cj * k[j]
    end
    return out
end

# In-place, with idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs, T::Type{Val{1}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt = inv(dt)

    dL1 = _lagrange_basis_deriv(Θ, 1, n, k)
    @views @.. out = (dL1 * invdt) * k[1][idxs]

    for j in 2:n
        dLj = _lagrange_basis_deriv(Θ, j, n, k)
        cj = dLj * invdt
        @views @.. out = out + cj * k[j][idxs]
    end
    return out
end

## FBDF Val{2}: Second derivative of Lagrange interpolant
#
# p''(Θ) = (1/dt²) * Σ_{i=1}^n k[i] * L''_i(Θ)

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{2}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt2 = inv(dt)^2

    d2L1 = _lagrange_basis_deriv2(Θ, 1, n, k)
    out = @.. (d2L1 * invdt2) * k[1]

    for j in 2:n
        d2Lj = _lagrange_basis_deriv2(Θ, j, n, k)
        cj = d2Lj * invdt2
        out = @.. out + cj * k[j]
    end

    if differential_vars !== nothing
        for i in eachindex(differential_vars)
            if !differential_vars[i]
                if out isa Number
                    out = zero(out)
                else
                    out[i] = zero(eltype(out))
                end
            end
        end
    end

    return out
end

# Out-of-place, with idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs, T::Type{Val{2}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt2 = inv(dt)^2

    d2L1 = _lagrange_basis_deriv2(Θ, 1, n, k)
    out = @.. (d2L1 * invdt2) * k[1][idxs]

    for j in 2:n
        d2Lj = _lagrange_basis_deriv2(Θ, j, n, k)
        cj = d2Lj * invdt2
        out = @.. out + cj * k[j][idxs]
    end

    if differential_vars !== nothing
        for (idx_pos, orig_idx) in enumerate(idxs)
            if !differential_vars[orig_idx]
                if out isa Number
                    out = zero(out)
                else
                    out[idx_pos] = zero(eltype(out))
                end
            end
        end
    end

    return out
end

# In-place, no idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{2}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt2 = inv(dt)^2

    d2L1 = _lagrange_basis_deriv2(Θ, 1, n, k)
    @.. out = (d2L1 * invdt2) * k[1]

    for j in 2:n
        d2Lj = _lagrange_basis_deriv2(Θ, j, n, k)
        cj = d2Lj * invdt2
        @.. out = out + cj * k[j]
    end

    if differential_vars !== nothing
        for i in eachindex(differential_vars)
            if !differential_vars[i]
                out[i] = zero(eltype(out))
            end
        end
    end

    return out
end

# In-place, with idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs, T::Type{Val{2}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt2 = inv(dt)^2

    d2L1 = _lagrange_basis_deriv2(Θ, 1, n, k)
    @views @.. out = (d2L1 * invdt2) * k[1][idxs]

    for j in 2:n
        d2Lj = _lagrange_basis_deriv2(Θ, j, n, k)
        cj = d2Lj * invdt2
        @views @.. out = out + cj * k[j][idxs]
    end

    if differential_vars !== nothing
        for (idx_pos, orig_idx) in enumerate(idxs)
            if !differential_vars[orig_idx]
                out[idx_pos] = zero(eltype(out))
            end
        end
    end

    return out
end

## FBDF Val{3}: Third derivative of Lagrange interpolant
#
# p'''(Θ) = (1/dt³) * Σ_{i=1}^n k[i] * L'''_i(Θ)

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{3}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt3 = inv(dt)^3

    d3L1 = _lagrange_basis_deriv3(Θ, 1, n, k)
    out = @.. (d3L1 * invdt3) * k[1]

    for j in 2:n
        d3Lj = _lagrange_basis_deriv3(Θ, j, n, k)
        cj = d3Lj * invdt3
        out = @.. out + cj * k[j]
    end

    if differential_vars !== nothing
        for i in eachindex(differential_vars)
            if !differential_vars[i]
                if out isa Number
                    out = zero(out)
                else
                    out[i] = zero(eltype(out))
                end
            end
        end
    end

    return out
end

# Out-of-place, with idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs, T::Type{Val{3}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt3 = inv(dt)^3

    d3L1 = _lagrange_basis_deriv3(Θ, 1, n, k)
    out = @.. (d3L1 * invdt3) * k[1][idxs]

    for j in 2:n
        d3Lj = _lagrange_basis_deriv3(Θ, j, n, k)
        cj = d3Lj * invdt3
        out = @.. out + cj * k[j][idxs]
    end

    if differential_vars !== nothing
        for (idx_pos, orig_idx) in enumerate(idxs)
            if !differential_vars[orig_idx]
                if out isa Number
                    out = zero(out)
                else
                    out[idx_pos] = zero(eltype(out))
                end
            end
        end
    end

    return out
end

# In-place, no idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{3}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt3 = inv(dt)^3

    d3L1 = _lagrange_basis_deriv3(Θ, 1, n, k)
    @.. out = (d3L1 * invdt3) * k[1]

    for j in 2:n
        d3Lj = _lagrange_basis_deriv3(Θ, j, n, k)
        cj = d3Lj * invdt3
        @.. out = out + cj * k[j]
    end

    if differential_vars !== nothing
        for i in eachindex(differential_vars)
            if !differential_vars[i]
                out[i] = zero(eltype(out))
            end
        end
    end

    return out
end

# In-place, with idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs, T::Type{Val{3}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt3 = inv(dt)^3

    d3L1 = _lagrange_basis_deriv3(Θ, 1, n, k)
    @views @.. out = (d3L1 * invdt3) * k[1][idxs]

    for j in 2:n
        d3Lj = _lagrange_basis_deriv3(Θ, j, n, k)
        cj = d3Lj * invdt3
        @views @.. out = out + cj * k[j][idxs]
    end

    if differential_vars !== nothing
        for (idx_pos, orig_idx) in enumerate(idxs)
            if !differential_vars[orig_idx]
                out[idx_pos] = zero(eltype(out))
            end
        end
    end

    return out
end
