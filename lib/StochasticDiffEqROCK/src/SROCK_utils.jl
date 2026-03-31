# This function calculates the largest eigenvalue
# (absolute value wise) by power iteration.
function maxeig!(integrator, cache::StochasticDiffEqConstantCache)
    isfirst = integrator.iter == 1 || integrator.u_modified
    (; t, dt, uprev, u, p) = integrator
    maxiter = 50
    safe = 1.2
    fsalfirst = integrator.f(uprev, p, t)
    # Initial guess for eigenvector `z`
    if isfirst
        fz = fsalfirst
        z = integrator.f(fz, p, t)
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
        z = dz_u .* (false .* z .+ one(eltype(z)))
    end # endif
    # Start power iteration
    integrator.eigen_est = 0
    for iter in 1:maxiter
        fz = integrator.f(z, p, t)
        fz = fz - fsalfirst

        Δ = integrator.opts.internalnorm(fz, t)
        eig_prev = integrator.eigen_est
        integrator.eigen_est = Δ / dz_u * safe
        # Convergence
        if iter >= 2 && abs(eig_prev - integrator.eigen_est) < integrator.eigen_est * 0.05
            # Store the eigenvector
            cache.zprev = z - uprev
            return true
        end

        # Next `z`
        if Δ != zero(Δ)
            quot = dz_u / Δ
            z = uprev + quot * fz
        else
            # An arbitrary change on `z`
            nind = length(z)
            if (nind != 1)
                ind = 1 + iter % nind
                # val = (uprev[ind] - (z[ind] - uprev[ind]))*one(eltype(z))*2
                vec(z) = vec(z) .* (1 .- 2 .* ((1:length(z)) .== ind))
            else
                z = -z
            end
        end
    end
    return false
end

function maxeig!(integrator, cache::StochasticDiffEqMutableCache)
    isfirst = integrator.iter == 1 || integrator.u_modified
    (; t, dt, uprev, u, p) = integrator
    fz, z, fsalfirst = cache.atmp, cache.tmp, cache.fsalfirst
    integrator.f(fsalfirst, uprev, p, t)
    ccache = cache.constantcache
    maxiter = 50
    safe = 1.2

    # Initial guess for eigenvector `z`
    if isfirst
        @.. fz = fsalfirst
        integrator.f(z, fz, p, t)
        # integrator.stats.nf += 1
    else
        @.. z = ccache.zprev
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
        @.. z = uprev + quot * z
    elseif !is_u_zero
        dz_u = u_norm * sqrt_pert
        @.. z = uprev + uprev * dz_u
    elseif !is_z_zero
        dz_u = pert
        quot = dz_u / z_norm
        @.. z *= quot
    else
        dz_u = pert
        @.. z = dz_u * (false * z + one(eltype(z)))
    end # endif
    # Start power iteration
    integrator.eigen_est = 0
    for iter in 1:maxiter
        integrator.f(fz, z, p, t)
        # integrator.stats.nf += 1
        @.. fz = fz - fsalfirst

        Δ = integrator.opts.internalnorm(fz, t)
        eig_prev = integrator.eigen_est
        integrator.eigen_est = Δ / dz_u * safe
        # Convergence
        if iter >= 2 && abs(eig_prev - integrator.eigen_est) < integrator.eigen_est * 0.05
            # Store the eigenvector
            @.. ccache.zprev = z - uprev
            return true
        end
        # Next `z`
        if Δ != zero(Δ)
            quot = dz_u / Δ
            @.. z = uprev + quot * fz
        else
            # An arbitrary change on `z`
            nind = length(uprev)
            if (nind != 1)
                ind = 1 + iter % nind
                # val = (uprev[ind] - (z[ind] - uprev[ind]))*one(eltype(z))
                vec(z) .= vec(z) .* (1 .- 2 .* ((1:length(z)) .== ind))
            else
                z = -z
            end
        end
    end
    return false
end

function choose_deg!(integrator, cache::T) where {T}
    isconst = T <: StochasticDiffEqConstantCache
    isconst || (cache = cache.constantcache)

    if integrator.alg isa SROCK1
        # binary search as stability domain is monotonically incrasing with number of stages
        mn_st, mx_st, mid_st = 3, 200, 3
        while (mx_st - mn_st > 1)
            mid_st = Int(floor((mn_st + mx_st) * 0.5))
            if (mid_st^2 * (131.97 / 197 - 0.45 * mid_st / 197) > cache.mdeg)
                mx_st = mid_st
            else
                mn_st = mid_st
            end
        end
        cache.mdeg = (cache.mdeg > mn_st^2 * (131.97 / 197 - 0.45 * mn_st / 197)) ? mx_st : mn_st

        @inbounds for i in 1:size(cache.ms, 1)
            if cache.ms[i] <= cache.mdeg
                cache.optimal_η = cache.mη[i]
            else
                break
            end
        end
    end

    if integrator.alg isa SROCK2
        start = 1
        @inbounds for i in 1:size(cache.ms, 1)
            if cache.ms[i] >= cache.mdeg
                cache.deg_index = i
                cache.mdeg = cache.ms[i]
                cache.start = start
                break
            else
                start += cache.ms[i] * 2 - 1
            end
        end
    end

    if integrator.alg isa TangXiaoSROCK2
        start = 1
        @inbounds for i in 1:size(cache.ms, 1)
            if cache.ms[i] >= cache.mdeg
                cache.deg_index = i
                cache.mdeg = cache.ms[i]
                cache.start = start
                break
            else
                start += cache.ms[i] * 2 - 1
            end
        end

        start = 1
        @inbounds for i in 1:5
            if integrator.alg.version_num == i
                cache.start_mcs = start
                break
            else
                start += cache.mn̂[i]
            end
        end
    end

    if integrator.alg isa SROCKEM
        @inbounds for i in 1:size(cache.ms, 1)
            if cache.ms[i] <= cache.mdeg
                cache.optimal_η = cache.mη[i]
            else
                break
            end
        end
    end

    if integrator.alg isa KomBurSROCK2
        start = 1
        @inbounds for i in 1:size(cache.ms, 1)
            if cache.ms[i] >= cache.mdeg
                cache.deg_index = i
                cache.mdeg = cache.ms[i]
                cache.start = start
                break
            else
                start += cache.ms[i] * 2 - 1
            end
        end
    end

    if integrator.alg isa SROCKC2
        start = 1
        @inbounds for i in 1:size(cache.ms, 1)
            if cache.ms[i] >= cache.mdeg
                cache.deg_index = i
                cache.mdeg = cache.ms[i]
                cache.start = start
                break
            else
                start += cache.ms[i] * 2 - 1
            end
        end
    end

    return nothing
end
