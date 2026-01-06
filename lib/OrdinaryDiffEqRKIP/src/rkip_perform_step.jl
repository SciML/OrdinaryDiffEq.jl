"""
Helper function for the (nonlinear) part of the problem maybe in place.
Overwrite compute f(u, p, t) and if possible, overwrite u with the result.
Same principle of operation as `matvec_prod_mip` for mutability/in-place handling.
"""
function nl_part_mip(
        tmp::uType, f!, u::uType, p, t::tType, ::Val{true}
    ) where {tType, uType}
    f!(tmp, u, p, t)
    copyto!(u, tmp)
    return u
end
function nl_part_mip(_::uType, f, u::uType, p, t::tType, ::Val{false}) where {tType, uType}
    return f(u, p, t)
end

"""
Helper function to compute Au + f(u, p, t) and if in place, store the result in res.
Return res if in place, otherwise return Au + f(u, p, t)
"""
function f_mip!(
        res::uType, tmp::uType, A::opType, f, u::uType, p,
        t::tType, ::Val{true}
    ) where {tType, uType, opType}
    res .= u
    res = matvec_prod_mip(tmp, A, res, Val(true), p, t)
    f(tmp, u, p, t)
    res .+= tmp
    return res
end

function f_mip!(
        _::uType, tmp::uType, A::opType, f, u::uType, p,
        t::tType, ::Val{false}
    ) where {tType, uType, opType}
    return matvec_prod_mip(tmp, A, u, Val(false), p, t) + f(u, p, t)
end

"""
Helper function for the residual maybe in place.
Same principle of operation as `_safe_matvec_prod` for mutability/in-place handling.
"""
function calculate_residuals_mip(
        tmp::uType, utilde::uType, uprev::uType, u::uType, abstol,
        reltol, internalnorm, t, ::Val{true}
    ) where {uType}
    calculate_residuals!(tmp, utilde, uprev, u, abstol, reltol, internalnorm, t)
    return tmp
end
function calculate_residuals_mip(
        _::uType, utilde::uType, uprev::uType, u::uType, abstol,
        reltol, internalnorm, t, ::Val{false}
    ) where {uType}
    return calculate_residuals(utilde, uprev, u, abstol, reltol, internalnorm, t)
end

@fastmath function perform_step!(
        integrator,
        cache::RKIPCache{expOpType, cacheType, tType, opType, uType, iip}
    ) where {
        expOpType, cacheType, tType, opType, uType, iip,
    }
    (; t, dt, uprev, u, f, p, fsalfirst, fsallast, alg) = integrator
    (; c, α, αEEst, stages, A) = alg.tableau
    (; kk, utilde, tmp) = cache
    (; adaptive, abstol, reltol, internalnorm) = integrator.opts

    Â::opType = f.f1.f # Linear Operator

    cache_exp_op_for_this_step!(cache, Â, dt, alg) # compute/reuse cached exp(± Â * [dt * cᵢ]) for this dt and cache them if possible

    @bb u .= uprev

    for i in 1:(stages)
        @bb kk[i] .= uprev

        for j in 1:(i - 1)
            g = A[i, j] * dt
            kk[i] = axpy_mip(g, kk[j], kk[i], iip) # kk_i += dt*A[i, j] * kk_j
        end

        t_ = t + c[i] * dt

        # if mutable/heaps type, assignment does nothing as the function is in-place,
        kk[i] = expmv_rkip_mip(cache, kk[i], dt, i, p, t_) # kᵢ = exp(Â * [dt * cᵢ])*kᵢ ➡ Change from interaction picture to "true" coordinate

        kk[i] = nl_part_mip(cache.tmp, f.f2, kk[i], p, t_, iip) # kᵢ = f(u + Σⱼ dt*(Aᵢⱼ kⱼ), t + dt*cᵢ)

        kk[i] = expmv_rkip_mip(cache, kk[i], -dt, i, p, t_) # kᵢ = exp(-Â * [dt * cᵢ])*kᵢ ➡ Going back in interaction picture

        integrator.stats.nf += 2 # two exp vec product
        integrator.stats.nf2 += 1 # one function evaluation

        u = axpy_mip(α[i] * dt, kk[i], u, iip) # uₙ = uₙ₋₁ + dt Σᵢ * (αᵢ - α*ᵢ) kᵢ

        if adaptive # error estimation ũ = Σᵢ dt * (αᵢ - α*ᵢ) kᵢ
            if i == 1
                @bb @. utilde = dt * (α[i] - αEEst[i]) * kk[i] # to avoid filling with zero, we set it to i == 1
            else
                utilde = axpy_mip(dt * (α[i] - αEEst[i]), kk[i], utilde, iip) # otherwise we add it
            end
        end
    end

    u = expmv_rkip_mip(cache, u, dt, p, t) # stepping forward into the interaction u = exp(Â dt)*u

    if adaptive
        utilde = expmv_rkip_mip(cache, utilde, dt, p, t) # stepping forward into the interaction ũ = exp(Â dt)*ũ
        tmp = calculate_residuals_mip(
            tmp, utilde, uprev, u, abstol, reltol, internalnorm, t, iip
        ) # error computation maybe in place
        integrator.EEst = internalnorm(tmp, t)
    end

    fsallast = f_mip!(fsallast, cache.tmp, Â, f.f2, u, p, t + dt, iip) # derivative estimation for interpolation
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = fsallast
    integrator.u = u

    @bb copyto!(integrator.k[1], fsalfirst)
    @bb copyto!(integrator.k[2], fsallast)
end

function initialize!(
        integrator,
        cache::RKIPCache{expOpType, cacheType, tType, opType, uType, iip}
    ) where {
        expOpType, cacheType, tType, opType, uType, iip,
    }
    (; f, u, p, t, fsalfirst, fsallast) = integrator

    kshortsize = 2
    k = [zero(u) for _ in 1:kshortsize]

    fsalfirst = f_mip!(fsalfirst, cache.tmp, f.f1.f, f.f2, u, p, t, iip) # first derivative for interpolation computation, maybe in place
    integrator.stats.nf += 1
    integrator.stats.nf2 += 1

    integrator.kshortsize = kshortsize
    integrator.k = k
    integrator.fsalfirst = fsalfirst
    integrator.fsallast = fsallast

    @bb copyto!(integrator.k[1], fsalfirst)
    return @bb copyto!(integrator.k[2], fsallast)
end
