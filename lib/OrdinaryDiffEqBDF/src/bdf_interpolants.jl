### Type unions for dispatch
QNDF_CACHES = Union{QNDFConstantCache, QNDFCache}
FBDF_CACHES = Union{FBDFConstantCache, FBDFCache, DFBDFConstantCache, DFBDFCache}
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
# FBDF / DFBDF: Lagrange interpolation at fixed Chebyshev nodes
#
# All data lives in k. y₀ and y₁ are NOT used.
#
# Layout: k has max_order+1 entries.
#   k[1..n] = solution values resampled at fixed Chebyshev reference
#             nodes xi_1,...,xi_n in [0,1] (see _CHEB_NODES in bdf_utils.jl)
#   k[n+1..end] = zero (unused entries for orders below max_order)
#
# The reference nodes are compile-time constants (no GPU scalar indexing).
# At step time the original Lagrange interpolant through the actual
# (variable) history nodes is resampled at these fixed positions.
# At interpolation time, the polynomial is reconstructed via the
# barycentric formula (Val{0}, O(n)) or product-rule Lagrange basis
# derivatives (Val{1..3}, O(n^2) but n≤6).
####################################################################

# Helper: determine number of active data points from k entries.
# Active entries are non-zero (counting from the top).
function _bdf_active_order(k)
    n = length(k)
    while n > 0 && iszero(k[n])
        n -= 1
    end
    return max(n, 1)
end

####################################################################
# Lagrange basis functions on fixed Chebyshev nodes (and derivatives)
#
# These use the precomputed _CHEB_NODES constants and operate entirely
# on scalars, so no GPU scalar indexing occurs.
####################################################################

# L_i(Θ) for fixed Chebyshev node i among n nodes.
@inline function _fbdf_lagrange_basis(Θ, i, n)
    nodes = _CHEB_NODES[n]
    Li = one(Θ)
    for m in 1:n
        m == i && continue
        Li *= (Θ - nodes[m]) / (nodes[i] - nodes[m])
    end
    return Li
end

# Product-rule form of L'_i(Θ):
#   L'_i(Θ) = Σ_{m≠i} [1/(xi_i - xi_m)] * Π_{l≠i,l≠m} (Θ - xi_l)/(xi_i - xi_l)
@inline function _fbdf_lagrange_basis_deriv(Θ, i, n)
    nodes = _CHEB_NODES[n]
    result = zero(Θ)
    for m in 1:n
        m == i && continue
        term = inv(nodes[i] - nodes[m])
        for l in 1:n
            (l == i || l == m) && continue
            term *= (Θ - nodes[l]) / (nodes[i] - nodes[l])
        end
        result += term
    end
    return result
end

# Product-rule form of L''_i(Θ):
@inline function _fbdf_lagrange_basis_deriv2(Θ, i, n)
    nodes = _CHEB_NODES[n]
    result = zero(Θ)
    for m in 1:n
        m == i && continue
        for p in 1:n
            (p == i || p == m) && continue
            term = inv(nodes[i] - nodes[m]) * inv(nodes[i] - nodes[p])
            for l in 1:n
                (l == i || l == m || l == p) && continue
                term *= (Θ - nodes[l]) / (nodes[i] - nodes[l])
            end
            result += term
        end
    end
    return result
end

# Product-rule form of L'''_i(Θ):
@inline function _fbdf_lagrange_basis_deriv3(Θ, i, n)
    nodes = _CHEB_NODES[n]
    result = zero(Θ)
    for m in 1:n
        m == i && continue
        for p in 1:n
            (p == i || p == m) && continue
            for q in 1:n
                (q == i || q == m || q == p) && continue
                term = inv(nodes[i] - nodes[m]) * inv(nodes[i] - nodes[p]) *
                    inv(nodes[i] - nodes[q])
                for l in 1:n
                    (l == i || l == m || l == p || l == q) && continue
                    term *= (Θ - nodes[l]) / (nodes[i] - nodes[l])
                end
                result += term
            end
        end
    end
    return result
end

####################################################################
# FBDF Val{0}: Barycentric Lagrange interpolation (O(n))
#
# p(Θ) = [Σ_j w_j/(Θ-xi_j) * k[j]] / [Σ_j w_j/(Θ-xi_j)]
#
# where w_j are precomputed barycentric weights and xi_j are
# Chebyshev reference nodes. When Θ equals a node, return k[j].
####################################################################

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    n = _bdf_active_order(k)
    nodes = _CHEB_NODES[n]
    weights = _BARY_WEIGHTS[n]

    # Check for exact node hit
    for j in 1:n
        if Θ == nodes[j]
            return k[j] isa Number ? k[j] : copy(k[j])
        end
    end

    # Barycentric formula: compute scalar coefficients, then linear combination
    den = zero(Θ)
    for j in 1:n
        den += weights[j] / (Θ - nodes[j])
    end
    inv_den = inv(den)

    c1 = weights[1] / (Θ - nodes[1]) * inv_den
    out = @.. c1 * k[1]
    for j in 2:n
        cj = weights[j] / (Θ - nodes[j]) * inv_den
        out = @.. out + cj * k[j]
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
    nodes = _CHEB_NODES[n]
    weights = _BARY_WEIGHTS[n]

    for j in 1:n
        if Θ == nodes[j]
            return k[j][idxs]
        end
    end

    den = zero(Θ)
    for j in 1:n
        den += weights[j] / (Θ - nodes[j])
    end
    inv_den = inv(den)

    c1 = weights[1] / (Θ - nodes[1]) * inv_den
    out = @.. c1 * k[1][idxs]
    for j in 2:n
        cj = weights[j] / (Θ - nodes[j]) * inv_den
        out = @.. out + cj * k[j][idxs]
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
    nodes = _CHEB_NODES[n]
    weights = _BARY_WEIGHTS[n]

    for j in 1:n
        if Θ == nodes[j]
            @.. out = k[j]
            return out
        end
    end

    den = zero(Θ)
    for j in 1:n
        den += weights[j] / (Θ - nodes[j])
    end
    inv_den = inv(den)

    c1 = weights[1] / (Θ - nodes[1]) * inv_den
    @.. out = c1 * k[1]
    for j in 2:n
        cj = weights[j] / (Θ - nodes[j]) * inv_den
        @.. out = out + cj * k[j]
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
    nodes = _CHEB_NODES[n]
    weights = _BARY_WEIGHTS[n]

    for j in 1:n
        if Θ == nodes[j]
            @views @.. out = k[j][idxs]
            return out
        end
    end

    den = zero(Θ)
    for j in 1:n
        den += weights[j] / (Θ - nodes[j])
    end
    inv_den = inv(den)

    c1 = weights[1] / (Θ - nodes[1]) * inv_den
    @views @.. out = c1 * k[1][idxs]
    for j in 2:n
        cj = weights[j] / (Θ - nodes[j]) * inv_den
        @views @.. out = out + cj * k[j][idxs]
    end
    return out
end

####################################################################
# FBDF Val{1}: First derivative via fixed-node Lagrange basis
#
# p'(Θ) = (1/dt) * Σ_{i=1}^n k[i] * L'_i(Θ)
####################################################################

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt = inv(dt)

    dL1 = _fbdf_lagrange_basis_deriv(Θ, 1, n)
    out = @.. (dL1 * invdt) * k[1]

    for j in 2:n
        dLj = _fbdf_lagrange_basis_deriv(Θ, j, n)
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

    dL1 = _fbdf_lagrange_basis_deriv(Θ, 1, n)
    out = @.. (dL1 * invdt) * k[1][idxs]

    for j in 2:n
        dLj = _fbdf_lagrange_basis_deriv(Θ, j, n)
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

    dL1 = _fbdf_lagrange_basis_deriv(Θ, 1, n)
    @.. out = (dL1 * invdt) * k[1]

    for j in 2:n
        dLj = _fbdf_lagrange_basis_deriv(Θ, j, n)
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

    dL1 = _fbdf_lagrange_basis_deriv(Θ, 1, n)
    @views @.. out = (dL1 * invdt) * k[1][idxs]

    for j in 2:n
        dLj = _fbdf_lagrange_basis_deriv(Θ, j, n)
        cj = dLj * invdt
        @views @.. out = out + cj * k[j][idxs]
    end
    return out
end

####################################################################
# FBDF Val{2}: Second derivative
#
# p''(Θ) = (1/dt²) * Σ_{i=1}^n k[i] * L''_i(Θ)
####################################################################

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{2}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt2 = inv(dt)^2

    d2L1 = _fbdf_lagrange_basis_deriv2(Θ, 1, n)
    out = @.. (d2L1 * invdt2) * k[1]

    for j in 2:n
        d2Lj = _fbdf_lagrange_basis_deriv2(Θ, j, n)
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

    d2L1 = _fbdf_lagrange_basis_deriv2(Θ, 1, n)
    out = @.. (d2L1 * invdt2) * k[1][idxs]

    for j in 2:n
        d2Lj = _fbdf_lagrange_basis_deriv2(Θ, j, n)
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

    d2L1 = _fbdf_lagrange_basis_deriv2(Θ, 1, n)
    @.. out = (d2L1 * invdt2) * k[1]

    for j in 2:n
        d2Lj = _fbdf_lagrange_basis_deriv2(Θ, j, n)
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

    d2L1 = _fbdf_lagrange_basis_deriv2(Θ, 1, n)
    @views @.. out = (d2L1 * invdt2) * k[1][idxs]

    for j in 2:n
        d2Lj = _fbdf_lagrange_basis_deriv2(Θ, j, n)
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

####################################################################
# FBDF Val{3}: Third derivative
#
# p'''(Θ) = (1/dt³) * Σ_{i=1}^n k[i] * L'''_i(Θ)
####################################################################

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{3}}, differential_vars
    )
    n = _bdf_active_order(k)
    invdt3 = inv(dt)^3

    d3L1 = _fbdf_lagrange_basis_deriv3(Θ, 1, n)
    out = @.. (d3L1 * invdt3) * k[1]

    for j in 2:n
        d3Lj = _fbdf_lagrange_basis_deriv3(Θ, j, n)
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

    d3L1 = _fbdf_lagrange_basis_deriv3(Θ, 1, n)
    out = @.. (d3L1 * invdt3) * k[1][idxs]

    for j in 2:n
        d3Lj = _fbdf_lagrange_basis_deriv3(Θ, j, n)
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

    d3L1 = _fbdf_lagrange_basis_deriv3(Θ, 1, n)
    @.. out = (d3L1 * invdt3) * k[1]

    for j in 2:n
        d3Lj = _fbdf_lagrange_basis_deriv3(Θ, j, n)
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

    d3L1 = _fbdf_lagrange_basis_deriv3(Θ, 1, n)
    @views @.. out = (d3L1 * invdt3) * k[1][idxs]

    for j in 2:n
        d3Lj = _fbdf_lagrange_basis_deriv3(Θ, j, n)
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
