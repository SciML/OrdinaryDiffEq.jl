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
function calc_finite_difference_weights(ts, t, order, ::Val{N}) where {N}
    max_order = N
    c = zero(MMatrix{max_order + 1, max_order + 1, eltype(ts)})
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

function estimate_terk!(integrator, cache, k, ::Val{max_order}) where {max_order}
    #calculate hᵏ⁻¹yᵏ⁻¹
    (; ts_tmp, terk_tmp, u_history) = cache
    (; t, dt, u) = integrator
    fd_weights = calc_finite_difference_weights(ts_tmp, t + dt, k - 1, Val(max_order))
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
