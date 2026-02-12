### Fallbacks to capture unsupported derivative orders
BDF_CACHES_WITH_INTERPOLATIONS = Union{
    QNDFConstantCache, QNDFCache,
    FBDFConstantCache, FBDFCache,
    DFBDFConstantCache, DFBDFCache,
}

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

####
# Newton backward difference interpolation
#
# p(Θ) = y₁ + Σ_{j=1}^{order} φ_j(Θ-1) * k[j]
#
# where φ_j(σ) = σ*(σ+1)*...*(σ+j-1)/j! and σ = Θ-1.
# Recurrence: φ₁(σ) = σ, φ_{j+1}(σ) = φ_j(σ) * (σ+j)/(j+1)
#
# At Θ=0: σ=-1, φ₁(-1)=-1, φ_j(-1)=0 for j≥2, so p(0) = y₁ - k[1] = y₀
# At Θ=1: σ=0, all φ_j(0)=0, so p(1) = y₁
####

## Val{0}: Function value interpolation

# Out-of-place, ConstantCache, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{
            QNDFConstantCache, FBDFConstantCache, DFBDFConstantCache,
        },
        idxs::Nothing, T::Type{Val{0}}, differential_vars
    )
    σ = Θ - 1
    num_k = length(k)
    # Compute via recurrence: φ₁ = σ, φ_{j+1} = φ_j * (σ+j)/(j+1)
    φ = σ
    out = @.. y₁ + φ * k[1]
    for j in 2:num_k
        φ *= (σ + j - 1) / j
        out = @.. out + φ * k[j]
    end
    return out
end

# Out-of-place, MutableCache, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{
            QNDFCache, FBDFCache, DFBDFCache,
        },
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
        cache::BDF_CACHES_WITH_INTERPOLATIONS,
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
        cache::BDF_CACHES_WITH_INTERPOLATIONS,
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
        cache::BDF_CACHES_WITH_INTERPOLATIONS,
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

## Val{1}: First derivative dp/dt = (1/dt) * Σ_{j=1}^{order} φ'_j(σ) * k[j]
#
# Recurrence for φ'_j:
#   φ'₁(σ) = 1
#   φ'_{j+1}(σ) = (φ'_j(σ) * (σ+j) + φ_j(σ)) / (j+1)

# Out-of-place, ConstantCache, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{
            QNDFConstantCache, FBDFConstantCache, DFBDFConstantCache,
        },
        idxs::Nothing, T::Type{Val{1}}, differential_vars
    )
    σ = Θ - 1
    num_k = length(k)
    invdt = inv(dt)
    # j=1: φ₁ = σ, φ'₁ = 1
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

# Out-of-place, MutableCache, no idxs
@muladd function _ode_interpolant(
        Θ, dt, y₀, y₁, k,
        cache::Union{
            QNDFCache, FBDFCache, DFBDFCache,
        },
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
        cache::BDF_CACHES_WITH_INTERPOLATIONS,
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
        cache::BDF_CACHES_WITH_INTERPOLATIONS,
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
        cache::BDF_CACHES_WITH_INTERPOLATIONS,
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
