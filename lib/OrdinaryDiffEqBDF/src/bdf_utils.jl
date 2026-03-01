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

function calc_R(ρ, k, ::Val{N}) where {N}
    R = zero(MMatrix{N, N, typeof(ρ)})
    @inbounds for r in 1:k
        R[1, r] = -r * ρ
        for j in 2:k
            R[j, r] = R[j - 1, r] * ((j - 1) - r * ρ) / j
        end
    end
    return SArray(R)
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

const γₖ = @SVector[sum(Int64(1) // j for j in 1:k) for k in 1:6]

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
    c = zero(MMatrix{max_order + 1, max_order + 1, eltype(ts)})
    calc_finite_difference_weights!(c, ts, t, order)
    return SArray(c)
end

function reinitFBDF!(integrator, cache)
    # This function is used to initialize arrays that store past history information.
    # It will be used in the first-time step advancing and event handling.
    (;
        consfailcnt, ts, u_history, u_corrector, iters_from_event,
        order,
    ) = cache
    (; t, dt, uprev) = integrator

    if integrator.u_modified
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

function estimate_terk(integrator, cache, k, ::Val{max_order}, u) where {max_order}
    (; ts_tmp, u_history) = cache
    (; t, dt) = integrator
    fd_weights = calc_finite_difference_weights(ts_tmp, t + dt, k - 1, Val(max_order))
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
