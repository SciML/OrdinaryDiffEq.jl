@inline function U!(k, U)
    @inbounds for r in 1:k
        U[1, r] = -r
        for j in 2:k
            U[j, r] = U[j - 1, r] * ((j - 1) - r) / j
        end
    end
    nothing
end

function R!(k, ρ, cache)
    @unpack R = cache
    @inbounds for r in 1:k
        R[1, r] = -r * ρ
        for j in 2:k
            R[j, r] = R[j - 1, r] * ((j - 1) - r * ρ) / j
        end
    end
    nothing
end

# This functions takes help of D2 array to create backward differences array D
# Ith row of D2 keeps Ith order backward differences (∇ⁱyₙ)

function backward_diff!(cache::OrdinaryDiffEqMutableCache, D, D2, k, flag = true)
    flag && copyto!(D[1], D2[1, 1])
    for i in 2:k
        for j in 1:(k - i + 1)
            @.. broadcast=false D2[i, j]=D2[i - 1, j] - D2[i - 1, j + 1]
        end
        flag && copyto!(D[i], D2[i, 1])
    end
end

function backward_diff!(cache::OrdinaryDiffEqConstantCache, D, D2, k, flag = true)
    flag && (D[1] = D2[1, 1])
    for i in 2:k
        for j in 1:(k - i + 1)
            D2[i, j] = D2[i - 1, j] - D2[i - 1, j + 1]
        end
        flag && (D[i] = D2[i, 1])
    end
end

# this function updates backward difference array D when stepsize gets change
# Formula -> D = D * (R * U)
# and it is taken from the paper -
# Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s paper
# E. Alberdi Celaya, J. J. Anza Aguirrezabala, and P. Chatzipantelidis
function reinterpolate_history!(cache::OrdinaryDiffEqMutableCache, D, R, k)
    @unpack tmp = cache.nlsolver
    fill!(tmp, zero(eltype(D[1])))
    for j in 1:k
        for k in 1:k
            @. tmp += D[k] * R[k, j]
        end
        D[j] .= tmp
        fill!(tmp, zero(eltype(tmp)))
    end
end

function reinterpolate_history!(cache::OrdinaryDiffEqConstantCache, D, R, k)
    tmp = zero(D[1])
    for j in 1:k
        for k in 1:k
            tmp += D[k] * R[k, j]
        end
        D[j] = tmp
    end
end

function calc_R(ρ, k, ::Val{N}) where {N}
    R = zero(MMatrix{N, N, typeof(ρ)})
    @inbounds for r in 1:k
        R[1, r] = -r * ρ
        for j in 2:k
            R[j, r] = R[j - 1, r] * ((j - 1) - r * ρ) / j
        end
    end
    SArray(R)
end

function update_D!(D, dd, k)
    dd = _vec(dd)
    @views @.. broadcast=false D[:, k + 2]=dd - D[:, k + 1]
    @views @.. broadcast=false D[:, k + 1]=dd
    for i in k:-1:1
        @views @.. broadcast=false D[:, i]=D[:, i] + D[:, i + 1]
    end
    return nothing
end

const γₖ = @SVector[sum(1 // j for j in 1:k) for k in 1:6]

error_constant(integrator, order) = error_constant(integrator, integrator.alg, order)
function error_constant(integrator, alg::QNDF, k)
    @unpack γₖ = integrator.cache
    κ = alg.kappa[k]
    κ * γₖ[k] + inv(k + 1)
end

function compute_weights!(ts, k, weights)
    for i in 1:(k + 1)
        weights[i] = one(eltype(weights))
        for j in 1:(i - 1)
            weights[i] *= ts[i] - ts[j]
        end
        for j in (i + 1):(k + 1)
            weights[i] *= ts[i] - ts[j]
        end
        weights[i] = 1 / weights[i]
    end
end

#This uses lagrangian interpolation to calculate the predictor
function calc_Lagrange_interp(k, weights, t, ts, u_history, u::Number)
    if t in ts
        i = searchsortedfirst(ts, t, rev = true)
        return u_history[i]
    else
        for i in 1:(k + 1)
            u += weights[i] / (t - ts[i]) * u_history[i]
        end
        for i in 1:(k + 1)
            u *= t - ts[i]
        end
    end
    u
end

function calc_Lagrange_interp(k, weights, t, ts, u_history, u)
    vc = _vec(u)
    if t in ts
        i = searchsortedfirst(ts, t, rev = true)
        return reshape(u_history[:, i], size(u))
    else
        for i in 1:(k + 1)
            vc += @.. broadcast=false weights[i] / (t - ts[i])*u_history[:, i]
            u = reshape(vc, size(u))
        end
        for i in 1:(k + 1)
            u *= @.. broadcast=false t-ts[i]
        end
    end
    u
end

function calc_Lagrange_interp!(k, weights, t, ts, u_history, u)
    vc = _vec(u)
    # ts is short enough that linear search will be faster
    i = findfirst(isequal(t), ts)
    if i !== nothing
        @.. broadcast=false @views vc = u_history[:, i]
    else
        for i in 1:(k + 1)
            wdt = weights[i] / (t - ts[i])
            @.. broadcast=false @views vc = muladd(wdt, u_history[:, i], vc)
        end
        for i in 1:(k + 1)
            dt = t - ts[i]
            @.. broadcast=false u*=dt
        end
    end
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
    # This function is used for initialize weights and arrays that store past history information. It will be used in the first-time step advancing and event handling.
    @unpack weights, consfailcnt, ts, u_history, u_corrector, iters_from_event, order = cache
    @unpack t, dt, uprev = integrator

    if integrator.u_modified
        order = cache.order = 1
        consfailcnt = cache.consfailcnt = cache.nconsteps = 0
        iters_from_event = cache.iters_from_event = 0

        fill!(weights, zero(eltype(weights)))
        fill!(ts, zero(eltype(ts)))
        fill!(u_history, zero(eltype(u_history)))
        fill!(u_corrector, zero(eltype(u_corrector)))
    end

    vuprev = _vec(uprev)
    @views if iters_from_event == 0
        weights[1] = 1 / dt
        ts[1] = t
        @.. broadcast=false u_history[:, 1]=vuprev
    elseif iters_from_event == 1 && t != ts[1]
        ts[2] = ts[1]
        ts[1] = t
        @.. broadcast=false u_history[:, 2]=u_history[:, 1]
        @.. broadcast=false u_history[:, 1]=vuprev
    elseif consfailcnt == 0
        for i in (order + 2):-1:2
            ts[i] = ts[i - 1]
            @.. broadcast=false u_history[:, i]=u_history[:, i - 1]
        end
        ts[1] = t
        @.. broadcast=false u_history[:, 1]=vuprev
    end
    if iters_from_event >= 1
        compute_weights!(ts, order, weights)
    end
end

function estimate_terk!(integrator, cache, k, ::Val{max_order}) where {max_order}
    #calculate hᵏ⁻¹yᵏ⁻¹
    @unpack ts_tmp, terk_tmp, u_history = cache
    @unpack t, dt, u = integrator
    fd_weights = calc_finite_difference_weights(ts_tmp, t + dt, k - 1, Val(max_order))
    @.. broadcast=false terk_tmp=fd_weights[1, k] * u
    vc = _vec(terk_tmp)
    for i in 2:k
        @.. broadcast=false @views vc += fd_weights[i, k] * u_history[:, i - 1]
    end
    @.. broadcast=false terk_tmp*=abs(dt^(k - 1))
end

function estimate_terk(integrator, cache, k, ::Val{max_order}, u) where {max_order}
    @unpack ts_tmp, u_history = cache
    @unpack t, dt = integrator
    fd_weights = calc_finite_difference_weights(ts_tmp, t + dt, k - 1, Val(max_order))
    terk = @.. broadcast=false fd_weights[1, k]*u
    #@show terk,fd_weights[1,k+1]
    if u isa Number
        for i in 2:k
            terk += fd_weights[i, k] * u_history[i - 1]
        end
        terk *= abs(dt^(k - 1))
    else
        vc = _vec(terk)
        if ArrayInterface.ismutable(vc)
            for i in 2:k
                @.. broadcast=false vc+=fd_weights[i, k] * @view(u_history[:, i - 1])
            end
        else
            for i in 2:k
                vc = @. vc + fd_weights[i, k] * @view(u_history[:, i - 1])
            end
        end
        terk = reshape(vc, size(terk))
        terk *= abs(dt^(k - 1))
    end
    return terk
end
