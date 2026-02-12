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
# FBDF / DFBDF: Lagrange interpolation
#
# k stores past solution values at equidistant points:
#   k[1] = u_{n-1} (= uprev) at normalized position x₁ = 0
#   k[j] = u_{n-j} at normalized position xⱼ = 1-j
# y₁ = u_n at normalized position x₀ = 1
#
# Direct Lagrange formula (no barycentric division singularities):
#   p(Θ) = Σ_{j=0}^n y_j L_j(Θ)
#   L_j(Θ) = Π_{m≠j} (Θ - x_m) / (x_j - x_m)
#
# Since x_m = 1-m, the denominators (x_j - x_m) = (m - j) are constant
# integers, avoiding division by (Θ - x_j). This ensures compatibility
# with ForwardDiff dual numbers.
#
# Same polynomial as calc_Lagrange_interp in bdf_utils.jl.
####################################################################

# Helper: determine active BDF order from k entries
# (unused entries are set to zero; active entries are past solution values)
function _bdf_active_order(k)
    n = length(k)
    while n > 0 && iszero(k[n])
        n -= 1
    end
    return max(n, 1)
end

## FBDF Val{0}: Direct Lagrange function value interpolation

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    num_k = _bdf_active_order(k)
    n = num_k
    σ = Θ - 1

    # L_0(Θ) = Π_{m=1}^n (σ + m) / m
    L0 = one(Θ)
    for m in 1:n
        L0 *= (σ + m) / m
    end
    out = @.. L0 * y₁

    # L_j(Θ) for j = 1..n
    for j in 1:n
        Lj = one(Θ)
        for m in 0:n
            if m != j
                Lj *= (σ + m) / (m - j)
            end
        end
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
    num_k = _bdf_active_order(k)
    n = num_k
    σ = Θ - 1

    L0 = one(Θ)
    for m in 1:n
        L0 *= (σ + m) / m
    end
    out = @.. L0 * y₁[idxs]

    for j in 1:n
        Lj = one(Θ)
        for m in 0:n
            if m != j
                Lj *= (σ + m) / (m - j)
            end
        end
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
    num_k = _bdf_active_order(k)
    n = num_k
    σ = Θ - 1

    L0 = one(Θ)
    for m in 1:n
        L0 *= (σ + m) / m
    end
    @.. out = L0 * y₁

    for j in 1:n
        Lj = one(Θ)
        for m in 0:n
            if m != j
                Lj *= (σ + m) / (m - j)
            end
        end
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
    num_k = _bdf_active_order(k)
    n = num_k
    σ = Θ - 1

    L0 = one(Θ)
    for m in 1:n
        L0 *= (σ + m) / m
    end
    @views @.. out = L0 * y₁[idxs]

    for j in 1:n
        Lj = one(Θ)
        for m in 0:n
            if m != j
                Lj *= (σ + m) / (m - j)
            end
        end
        @views @.. out = out + Lj * k[j][idxs]
    end

    return out
end

## FBDF Val{1}: First derivative via Newton backward difference conversion
#
# Solution values in k are converted to backward differences as scalar
# coefficients, then the Newton derivative formula is applied.
#
# p'(Θ) = (1/dt) * Σ dφ_j * D_j
# where D_j = y₁ + Σ_{i=1}^j (-1)^i C(j,i) k[i]
#
# Expanding: p'(Θ) = (α/dt)*y₁ + Σ (β_i/dt)*k[i]
# where α = Σ dφ_j, β_i = Σ_{j≥i} dφ_j * (-1)^i * C(j,i)

# Out-of-place, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    num_k = _bdf_active_order(k)
    σ = Θ - 1
    invdt = inv(dt)

    # Pre-compute scalar coefficients
    φ = σ
    dφ = one(Θ)
    α = dφ
    β = zeros(typeof(Θ), num_k)
    β[1] += dφ * (-one(Θ))

    for j in 2:num_k
        dφ_new = (dφ * (σ + j - 1) + φ) / j
        φ *= (σ + j - 1) / j
        dφ = dφ_new
        α += dφ
        for i in 1:j
            β[i] += dφ * ((-1)^i * binomial(j, i))
        end
    end

    c_y1 = α * invdt
    out = @.. c_y1 * y₁
    for i in 1:num_k
        ci = β[i] * invdt
        out = @.. out + ci * k[i]
    end
    return out
end

# Out-of-place, with idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs, T::Type{Val{1}}, differential_vars
    )
    num_k = _bdf_active_order(k)
    σ = Θ - 1
    invdt = inv(dt)

    φ = σ
    dφ = one(Θ)
    α = dφ
    β = zeros(typeof(Θ), num_k)
    β[1] += dφ * (-one(Θ))

    for j in 2:num_k
        dφ_new = (dφ * (σ + j - 1) + φ) / j
        φ *= (σ + j - 1) / j
        dφ = dφ_new
        α += dφ
        for i in 1:j
            β[i] += dφ * ((-1)^i * binomial(j, i))
        end
    end

    c_y1 = α * invdt
    out = @.. c_y1 * y₁[idxs]
    for i in 1:num_k
        ci = β[i] * invdt
        out = @.. out + ci * k[i][idxs]
    end
    return out
end

# In-place, no idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    num_k = _bdf_active_order(k)
    σ = Θ - 1
    invdt = inv(dt)

    φ = σ
    dφ = one(Θ)
    α = dφ
    β = zeros(typeof(Θ), num_k)
    β[1] += dφ * (-one(Θ))

    for j in 2:num_k
        dφ_new = (dφ * (σ + j - 1) + φ) / j
        φ *= (σ + j - 1) / j
        dφ = dφ_new
        α += dφ
        for i in 1:j
            β[i] += dφ * ((-1)^i * binomial(j, i))
        end
    end

    c_y1 = α * invdt
    @.. out = c_y1 * y₁
    for i in 1:num_k
        ci = β[i] * invdt
        @.. out = out + ci * k[i]
    end
    return out
end

# In-place, with idxs
@muladd function _ode_interpolant!(
        out, Θ, dt, y₀, y₁, k,
        cache::FBDF_CACHES,
        idxs, T::Type{Val{1}}, differential_vars
    )
    num_k = _bdf_active_order(k)
    σ = Θ - 1
    invdt = inv(dt)

    φ = σ
    dφ = one(Θ)
    α = dφ
    β = zeros(typeof(Θ), num_k)
    β[1] += dφ * (-one(Θ))

    for j in 2:num_k
        dφ_new = (dφ * (σ + j - 1) + φ) / j
        φ *= (σ + j - 1) / j
        dφ = dφ_new
        α += dφ
        for i in 1:j
            β[i] += dφ * ((-1)^i * binomial(j, i))
        end
    end

    c_y1 = α * invdt
    @views @.. out = c_y1 * y₁[idxs]
    for i in 1:num_k
        ci = β[i] * invdt
        @views @.. out = out + ci * k[i][idxs]
    end
    return out
end
