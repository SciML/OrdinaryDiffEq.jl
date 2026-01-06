# This function calculates the largest eigenvalue
# (absolute value wise) by power iteration.
const RKCAlgs = Union{RKC, ESERK4, ESERK5, SERK2, ROCK2, ROCK4}

function maxeig!(integrator, cache::OrdinaryDiffEqConstantCache)
    isfirst = integrator.iter == 1 || integrator.u_modified
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    maxiter = (integrator.alg isa Union{ESERK4, ESERK5, SERK2}) ? 100 : 50

    safe = (integrator.alg isa RKCAlgs) ? 1.0 : 1.2
    # Initial guess for eigenvector `z`
    if isfirst
        if integrator.alg isa RKCAlgs
            z = fsalfirst
        else
            fz = fsalfirst
            z = f(fz, p, t)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
    else
        z = cache.zprev
    end
    # Perturbation
    u_norm = integrator.opts.internalnorm(uprev, t)
    z_norm = integrator.opts.internalnorm(z, t)
    pert = eps(u_norm)
    sqrt_pert = sqrt(pert)
    is_u_zero = u_norm == zero(u_norm)
    is_z_zero = z_norm == zero(z_norm)
    # Normalize `z` such that z-u lie in a circle
    if (!is_u_zero && !is_z_zero)
        dz_u = u_norm * sqrt_pert
        quot = dz_u / z_norm
        z = uprev + quot * z
    elseif !is_u_zero
        dz_u = u_norm * sqrt_pert
        z = uprev + uprev * dz_u
    elseif !is_z_zero
        dz_u = pert
        quot = dz_u / z_norm
        z *= quot
    else
        dz_u = pert
        z = dz_u * ones(z)
    end # endif
    # Start power iteration
    integrator.eigen_est = 0
    for iter in 1:maxiter
        fz = f(z, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        tmp = fz - fsalfirst
        Δ = integrator.opts.internalnorm(tmp, t)
        eig_prev = integrator.eigen_est
        integrator.eigen_est = Δ / dz_u * safe
        # Convergence
        if integrator.alg isa RKCAlgs # To match the constants given in the paper
            if iter >= 2 &&
                    abs(eig_prev - integrator.eigen_est) <
                    max(integrator.eigen_est, 1.0 / integrator.opts.dtmax) * 0.01
                integrator.eigen_est *= 1.2
                # Store the eigenvector
                cache.zprev = z - uprev
                return true
            end
        else
            if iter >= 2 &&
                    abs(eig_prev - integrator.eigen_est) < integrator.eigen_est * 0.05
                # Store the eigenvector
                cache.zprev = z - uprev
                return true
            end
        end

        # Next `z`
        if Δ != zero(Δ)
            quot = dz_u / Δ
            z = uprev + quot * tmp
        else
            # An arbitrary change on `z`
            nind = length(z)
            if (nind != 1)
                ind = 1 + iter % nind
                # val = (uprev[ind] - (z[ind] - uprev[ind]))*one(eltype(z))*2
                _vec(z) .= _vec(z) .* (1 .- 2 .* ((1:length(z)) .== ind))
            else
                z = -z
            end
        end
    end
    return false
end

function maxeig!(integrator, cache::OrdinaryDiffEqMutableCache)
    isfirst = integrator.iter == 1 || integrator.u_modified
    (; t, dt, uprev, u, f, p, fsalfirst) = integrator
    fz, z, atmp = cache.k, cache.tmp, cache.atmp
    ccache = cache.constantcache
    maxiter = (integrator.alg isa Union{ESERK4, ESERK5, SERK2}) ? 100 : 50
    safe = (integrator.alg isa RKCAlgs) ? 1.0 : 1.2
    # Initial guess for eigenvector `z`
    if isfirst
        if integrator.alg isa RKCAlgs
            @.. broadcast = false z = fsalfirst
        else
            @.. broadcast = false fz = fsalfirst
            f(z, fz, p, t)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
    else
        @.. broadcast = false z = ccache.zprev
    end
    # Perturbation
    u_norm = integrator.opts.internalnorm(uprev, t)
    z_norm = integrator.opts.internalnorm(z, t)
    pert = eps(u_norm)
    sqrt_pert = sqrt(pert)
    is_u_zero = u_norm == zero(u_norm)
    is_z_zero = z_norm == zero(z_norm)
    # Normalize `z` such that z-u lie in a circle
    if (!is_u_zero && !is_z_zero)
        dz_u = u_norm * sqrt_pert
        quot = dz_u / z_norm
        @.. broadcast = false z = uprev + quot * z
    elseif !is_u_zero
        dz_u = u_norm * sqrt_pert
        @.. broadcast = false z = uprev + uprev * dz_u
    elseif !is_z_zero
        dz_u = pert
        quot = dz_u / z_norm
        @.. broadcast = false z *= quot
    else
        dz_u = pert
        @.. broadcast = false z = dz_u * one(eltype(z))
    end # endif
    # Start power iteration
    integrator.eigen_est = 0
    for iter in 1:maxiter
        f(fz, z, p, t)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        @.. broadcast = false atmp = fz - fsalfirst
        Δ = integrator.opts.internalnorm(atmp, t)
        eig_prev = integrator.eigen_est
        integrator.eigen_est = Δ / dz_u * safe
        # Convergence
        if integrator.alg isa RKCAlgs # To match the constants given in the paper
            if iter >= 2 &&
                    abs(eig_prev - integrator.eigen_est) <
                    max(integrator.eigen_est, 1.0 / integrator.opts.dtmax) * 0.01
                integrator.eigen_est *= 1.2
                # Store the eigenvector
                @.. broadcast = false ccache.zprev = z - uprev
                return true
            end
        else
            if iter >= 2 &&
                    abs(eig_prev - integrator.eigen_est) < integrator.eigen_est * 0.05
                # Store the eigenvector
                @.. broadcast = false ccache.zprev = z - uprev
                return true
            end
        end
        # Next `z`
        if Δ != zero(Δ)
            quot = dz_u / Δ
            @.. broadcast = false z = uprev + quot * atmp
        else
            # An arbitrary change on `z`
            nind = length(uprev)
            if (nind != 1)
                ind = 1 + iter % nind
                # val = (uprev[ind] - (z[ind] - uprev[ind]))*one(eltype(z))
                _vec(z) .= _vec(z) .* (1 .- 2 .* ((1:length(z)) .== ind))
            else
                z = -z
            end
        end
    end
    return false
end
"""
    choosedeg!(cache) -> nothing

Calculate `mdeg = ms[deg_index]` (the degree of the Chebyshev polynomial)
and `cache.start` (the start index of recurrence parameters for that
degree), where `recf` are the `μ,κ` pairs
for the `mdeg` degree method. The `κ` for `stage-1` for every degree
is 0 therefore it's not included in `recf`
"""
function choosedeg!(cache::T) where {T}
    isconst = T <: OrdinaryDiffEqConstantCache
    isconst || (cache = cache.constantcache)
    start = 1
    @inbounds for i in 1:size(cache.ms, 1)
        if cache.ms[i] >= cache.mdeg
            cache.deg_index = i
            cache.mdeg = cache.ms[i]
            cache.start = start
            break
        end
        start += cache.ms[i] * 2 - 1
    end
    return nothing
end

function choosedeg_SERK!(integrator, cache::T) where {T}
    isconst = T <: OrdinaryDiffEqConstantCache
    isconst || (cache = cache.constantcache)
    (; ms) = cache
    start = 1
    @inbounds for i in 1:size(ms, 1)
        if ms[i] < cache.mdeg
            start += ms[i] + 1
        else
            cache.start = start
            cache.mdeg = ms[i]
            break
        end
    end
    if integrator.alg isa ESERK5
        if cache.mdeg <= 20
            cache.internal_deg = 2
        elseif cache.mdeg <= 50
            cache.internal_deg = 5
        elseif cache.mdeg <= 100
            cache.internal_deg = 10
        elseif cache.mdeg <= 500
            cache.internal_deg = 50
        elseif cache.mdeg <= 1000
            cache.internal_deg = 100
        elseif cache.mdeg <= 2000
            cache.internal_deg = 200
        end
    end

    if integrator.alg isa ESERK4
        if cache.mdeg <= 20
            cache.internal_deg = 2
        elseif cache.mdeg <= 100
            cache.internal_deg = 10
        elseif cache.mdeg <= 500
            cache.internal_deg = 25
        elseif cache.mdeg <= 1000
            cache.internal_deg = 100
        elseif cache.mdeg <= 4000
            cache.internal_deg = 200
        end
    end

    if integrator.alg isa SERK2
        cache.internal_deg = cache.mdeg / 10
    end
    return nothing
end
