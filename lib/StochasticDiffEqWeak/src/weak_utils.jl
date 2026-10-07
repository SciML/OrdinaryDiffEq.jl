################################################################################
# Ihat2 methods for weak-order cache types
# These extend StochasticDiffEqCore.Ihat2 with cache-type-specific dispatch.

function Ihat2(cache::Union{DRI1ConstantCache, DRI1Cache}, _dW, _dZ, sqdt, k, l)
    # compute elements of I^_(k,l) which is a mxm matrix
    if k < l
        return (_dW[k] * _dW[l] - sqdt * _dZ[k]) / 2
    elseif l < k
        return (_dW[k] * _dW[l] + sqdt * _dZ[l]) / 2
    else
        # l == k
        return (_dW[k]^2 - sqdt^2) / 2
    end
end

function Ihat2(cache::Union{RSConstantCache, RSCache}, _dW, _dZ, sqdt, k, l)
    # compute elements of I^_(k,l) which is a mxm matrix
    if k < l
        return -_dW[k] * _dZ[l]
    elseif l < k
        return _dW[k] * _dZ[l]
    else
        # l == k
        return zero(_dW[k])
    end
end

function Ihat2(cache::Union{PL1WMConstantCache, PL1WMCache}, _dW, _dZ, sqdt, k, l)
    # compute elements of I^_(k,l) which is a mxm matrix
    if k < l
        return -_dZ[Int(1 + 1 // 2 * (l - 3) * l + k)]
    elseif l < k
        return _dZ[Int(1 + 1 // 2 * (k - 3) * k + l)]
    else
        # l == k
        return -sqdt^2
    end
end

function Ihat2(cache::Union{NONConstantCache, NONCache}, _dW, _dZ, sqdt, k, l)
    # compute elements of I^_(k,l) which is a mxm matrix
    if k < l
        return _dZ[k]
    elseif l < k
        return _dW[l]
    else
        # l == k
        return _dW[k]
    end
end

################################################################################
# resize!/deleteat!/addat! support for the per-noise-dimension stage buffers.
# With diagonal noise the noise dimension equals the state length, so the outer
# vectors grow and shrink with the state. Every element is scratch space that
# `perform_step!` overwrites before reading, so only the lengths need fixing.

noise_stage_vectors(c::Union{DRI1Cache, IRI1Cache}) =
    (c.g2, c.g3, c.H12, c.H13, c.H22, c.H23)
noise_stage_vectors(c::RSCache) = (c.g2, c.g3, c.g4, c.H12, c.H13, c.H14, c.H22, c.H23)
noise_stage_vectors(c::W2Ito1Cache) = (c.H12, c.H13)
noise_stage_vectors(c::PL1WMCache) = (c.Yp, c.Ym)

const NoiseStageCache = Union{DRI1Cache, IRI1Cache, RSCache, W2Ito1Cache, PL1WMCache}

function fit_stage_vector!(v, n)
    isempty(v) && return v
    while length(v) > n
        pop!(v)
    end
    while length(v) < n
        push!(v, similar(v[1]))
    end
    for x in v
        length(x) == n || resize!(x, n)
        fill!(x, zero(eltype(x)))
    end
    return v
end

function fit_noise_stages!(integrator, cache::NoiseStageCache)
    is_diagonal_noise(integrator.sol.prob) || return nothing
    n = length(integrator.u)
    foreach(v -> fit_stage_vector!(v, n), noise_stage_vectors(cache))
    return nothing
end

for f in (:resize_non_user_cache!, :deleteat_non_user_cache!, :addat_non_user_cache!)
    @eval function SciMLBase.$f(integrator::SDEIntegrator, cache::NoiseStageCache, i)
        invoke(SciMLBase.$f, Tuple{SDEIntegrator, Any, Any}, integrator, cache, i)
        return fit_noise_stages!(integrator, cache)
    end
end
