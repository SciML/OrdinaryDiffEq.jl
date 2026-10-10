function nordsieck_adjust!(integrator, cache::T) where {T}
    (; nextorder, order) = cache
    if nextorder != order
        if is_bdf(integrator.alg)
            nordsieck_adjust_order_bdf!(cache, nextorder - order)
        else
            nordsieck_adjust_order!(cache, nextorder - order)
        end
        cache.order = cache.nextorder
        cache.L = cache.order + 1
        cache.n_wait = cache.L
    end
    nordsieck_rescale!(cache)
    return nothing
end

# TODO: SUNDIALS NLsolve handling

function nordsieck_finalize!(integrator, cache::T) where {T}
    isconst = T <: OrdinaryDiffEqConstantCache
    (; order, dts) = cache
    update_nordsieck_vector!(cache)
    cache.n_wait -= 1
    return if is_nordsieck_change_order(cache, 1) && cache.order != nordsieck_qmax(integrator.alg)
        if isconst
            cache.z[end] = cache.Δ
        else
            @.. broadcast = false cache.z[end] = cache.Δ
        end
        cache.prev_𝒟 = cache.c_𝒟
    end
end

function nordsieck_prepare_next!(integrator, cache::T) where {T}
    isconst = T <: OrdinaryDiffEqConstantCache
    (; maxη, order, L) = cache
    # TODO: further clean up
    (; bias1, bias2, bias3, addon) = integrator.alg
    if get_EEst(integrator) > one(get_EEst(integrator))
        nordsieck_rewind!(cache)
        cache.n_wait = max(2, cache.n_wait)
        cache.nextorder = order
        cache.η = inv((bias2 * get_EEst(integrator))^inv(L) + addon)
        return nothing
    end
    cache.ηq = inv((bias2 * get_EEst(integrator))^inv(L) + addon)
    stepsize_η!(integrator, cache, cache.order)
    if !is_nordsieck_change_order(cache)
        cache.η = cache.ηq
        cache.nextorder = order
        setη!(integrator, cache)
        return nothing
    end
    # On an order change (cache.n_wait == 0), we are going to compute the η for
    # order q+1 and q-1, where η = dt_next/dt
    cache.n_wait = 2
    stepsize_η₊₁!(integrator, cache, order)
    stepsize_η₋₁!(integrator, cache, order)
    chooseη!(integrator, cache)
    setη!(integrator, cache)
    # TODO: Maybe not here
    if isconst
        cache.Δ = cache.c_LTE * cache.Δ
    else
        @.. broadcast = false cache.Δ = cache.c_LTE * cache.Δ
    end
    return nothing
end

function nordsieck_complete_step!(integrator, cache, oldest_dt)
    if get_EEst(integrator) > one(get_EEst(integrator))
        nordsieck_prepare_next!(integrator, cache)
        nordsieck_restore_dts!(cache, oldest_dt)
    else
        nordsieck_finalize!(integrator, cache)
        nordsieck_prepare_next!(integrator, cache)
    end
    return nothing
end

##############################################################
# Lower level functions
##############################################################

# This function computes the integral, from -1 to 0, of a polynomial
# `P(x)` from the coefficients of `P` with an offset `k`.
function ∫₋₁⁰dx(a, deg, k)
    @inbounds begin
        int = zero(eltype(a))
        sign = 1
        for i in 0:deg
            int += flipsign(a[i + 1] / (i + k), sign)
            sign = -sign
        end
        return int
    end
end

# `l` is the coefficients of the polynomial `Λ` that satisfies conditions
# Λ(0) = 1, Λ(-1) = 0, and Λ̇(-ξᵢ) = 0, where ξᵢ = (tₙ-tₙ₋₁)/dt.
# It is described in the paper "A Polyalgorithm for the Numerical Solution
# of Ordinary Differential Equations" by G. D. Byrne and A. C. Hindmarsh in
# the page 86.
# https://dl.acm.org/citation.cfm?id=355636

# More implementation details are in the
# https://github.com/JuliaDiffEq/DiffEqDevMaterials repository
function calc_coeff!(cache::T) where {T}
    isvode = (T <: JVODECache || T <: JVODEConstantCache)
    @inbounds begin
        isconst = T <: OrdinaryDiffEqConstantCache
        isvarorder = is_nordsieck_change_order(cache, 1)
        (; m, l, dts, order) = cache
        dtsum = dt = dts[1]
        if order == 1
            l[1] = l[2] = cache.c_LTE₋₁ = cache.c_𝒟 = 1
            cache.c_LTE = 1 // 2
            cache.c_LTE₊₁ = 1 // 12
            cache.c_conv = 1 // 10 / cache.c_LTE
            return nothing
        end
        m[1] = 1
        for i in 2:(order + 1)
            m[i] = 0
        end
        # initialize ξ_inv
        ξ_inv = dt / dtsum
        # compute coefficients from the Newton polynomial
        # check the `JuliaDiffEq/DiffEqDevMaterials` repository for more details
        for j in 1:(order - 1)
            if isvarorder && j == order - 1
                M₋₁ = ∫₋₁⁰dx(m, order - 2, 2)
                # It is the same with `tq[1]` in SUNDIALS cvode.c
                cache.c_LTE₋₁ = order * M₋₁ / m[order - 1]
            end
            ξ_inv = dt / dtsum
            for i in j:-1:1
                m[i + 1] = muladd(m[i], ξ_inv, m[i + 1])
            end
            dtsum += dts[j + 1]
        end
        ξ_inv = dt / dtsum

        M0 = ∫₋₁⁰dx(m, order - 1, 1)
        M1 = ∫₋₁⁰dx(m, order - 1, 2)
        M0_inv = inv(M0)
        l[1] = 1
        for i in 1:order
            l[i + 1] = M0_inv * m[i] / i
        end
        # TODO: simplify LTE calculation
        # This is the error estimation coefficient for the current order `q`
        # ||Δ||⋅c_LTE yields the difference between a `q` degree interpolating
        # polynomial and a `q+1` degree interpolating polynomial at time `t`.
        # It is the same with `tq[2]` in SUNDIALS cvode.c
        cache.c_LTE = M1 * M0_inv * ξ_inv
        # It is the same with `tq[5]` in SUNDIALS cvode.c
        isvode && (cache.c_𝒟 = inv(ξ_inv) / l[order + 1])
        if isvarorder
            for i in order:-1:1
                m[i + 1] = muladd(ξ_inv, m[i], m[i + 1])
            end
            M2 = ∫₋₁⁰dx(m, order, 2)
            # It is the same with `tq[3]` in SUNDIALS cvode.c
            cache.c_LTE₊₁ = M2 * M0_inv / (order + 1)
        end # endif isvarorder
        # It is the same with `tq[4]` in SUNDIALS cvode.c
        cache.c_conv = 1 // 10 / cache.c_LTE
        return nothing
    end # end @inbounds
end

# Fixed-leading-coefficient BDF coefficients `l` and error constants, following
# `cvSetBDF` and `cvSetTqBDF` in SUNDIALS CVODE. CVODE's `l[i]` is `l[i + 1]` here
# and its `tau[j]` (the j-th previous step size) is `dts[j + 1]`.
function calc_coeff_bdf!(cache)
    @inbounds begin
        (; l, dts, order) = cache
        q = order
        dt = dts[1]
        T = eltype(l)
        for i in 1:(q + 1)
            l[i] = zero(T)
        end
        l[1] = l[2] = one(T)
        ξ_inv = ξstar_inv = one(T)
        α0 = α0_hat = -one(T)
        dtsum = dt
        if q > 1
            for j in 2:(q - 1)
                dtsum += dts[j]
                ξ_inv = dt / dtsum
                α0 -= one(T) / j
                for i in j:-1:1
                    l[i + 1] = muladd(l[i], ξ_inv, l[i + 1])
                end
            end
            α0 -= one(T) / q
            ξstar_inv = -l[2] - α0
            dtsum += dts[q]
            ξ_inv = dt / dtsum
            α0_hat = -l[2] - ξ_inv
            for i in q:-1:1
                l[i + 1] = muladd(l[i], ξstar_inv, l[i + 1])
            end
        end
        A1 = 1 - α0_hat + α0
        A2 = 1 + q * A1
        cache.c_LTE = abs(A1 / (α0 * A2))
        cache.c_𝒟 = abs(A2 * ξstar_inv / (l[q + 1] * ξ_inv))
        if is_nordsieck_change_order(cache, 1)
            if q > 1
                A3 = α0 + one(T) / q
                A4 = α0_hat + ξ_inv
                cache.c_LTE₋₁ = abs(ξstar_inv / l[q + 1] * (1 - A4 + A3) / A3)
            else
                cache.c_LTE₋₁ = one(T)
            end
            dtsum += dts[q + 1]
            ξ_inv = dt / dtsum
            A5 = α0 - one(T) / (q + 1)
            A6 = α0_hat - ξ_inv
            cache.c_LTE₊₁ = abs((1 - A6 + A5) / A2 / (ξ_inv * (q + 2) * A5))
        end
        cache.c_conv = 1 // 10 / cache.c_LTE
        return nothing
    end
end

# Apply the Pascal linear operator
function perform_predict!(cache::T, rewind = false) where {T}
    return @inbounds begin
        isconst = T <: OrdinaryDiffEqConstantCache
        (; z, order) = cache
        # This can be parallelized
        if !rewind
            if isconst
                for i in 1:order, j in order:-1:i

                    z[j] = z[j] + z[j + 1]
                end
            else
                for i in 1:order, j in order:-1:i

                    @.. broadcast = false z[j] = z[j] + z[j + 1]
                end
            end # endif const cache
        else
            if isconst
                for i in 1:order, j in order:-1:i

                    z[j] = z[j] - z[j + 1]
                end
            else
                for i in 1:order, j in order:-1:i

                    @.. broadcast = false z[j] = z[j] - z[j + 1]
                end
            end # endif const cache
        end # endif !rewind
    end # end @inbounds
end

# Apply corrections on the Nordsieck vector
function update_nordsieck_vector!(cache::T) where {T}
    isvode = (T <: JVODECache || T <: JVODEConstantCache)
    return @inbounds begin
        isconst = T <: OrdinaryDiffEqConstantCache
        (; z, Δ, l, order) = cache
        if isconst
            for i in 1:(order + 1)
                z[i] = muladd.(l[i], Δ, z[i])
            end
        else
            for i in 1:(order + 1)
                @.. broadcast = false z[i] = muladd(l[i], Δ, z[i])
            end
        end # endif not const cache
    end # end @inbounds
end

function nlsolve_functional!(integrator, cache::T) where {T}
    (; f, dt, t, p) = integrator
    isconstcache = T <: OrdinaryDiffEqConstantCache
    (; z, l, c_conv, Δ) = cache
    if isconstcache
        ratetmp = integrator.f(z[1], p, dt + t)
    else
        (; ratetmp) = cache
        integrator.f(ratetmp, z[1], p, dt + t)
    end
    increment_nf!(integrator.stats, 1)
    maxiters = 3
    div_rate = 2
    # Zero out the difference vector
    isconstcache ? (cache.Δ = zero(cache.Δ)) : (Δ .= zero(eltype(Δ)))
    # `k` is a counter for convergence test
    k = 0
    # `conv_rate` is used in convergence rate estimation
    conv_rate = 1.0
    # initialize `δ_prev`
    δ_prev = 0
    # Start the functional iteration & store the difference into `Δ`
    for k in 1:maxiters
        if isconstcache
            ratetmp = inv(l[2]) * muladd.(dt, ratetmp, -z[2])
            integrator.u = ratetmp + z[1]
            cache.Δ = ratetmp - cache.Δ
        else
            @.. broadcast = false integrator.u = -z[2]
            @.. broadcast = false ratetmp = inv(l[2]) * muladd(dt, ratetmp, integrator.u)
            @.. broadcast = false integrator.u = ratetmp + z[1]
            @.. broadcast = false cache.Δ = ratetmp - cache.Δ
        end
        # @show norm(dt*ratetmp - ( z[2] + (integrator.u - z[1])*l[2] ))
        # @show norm(cache.Δ - (integrator.u - z[1]))
        # It only makes sense to calculate convergence rate in the second iteration
        δ = integrator.opts.internalnorm(cache.Δ, t)
        isconstcache ? (cache.Δ = copy(ratetmp)) : copyto!(cache.Δ, ratetmp)
        if k >= 1
            conv_rate = max(1 // 10 * conv_rate, δ / δ_prev)
        end
        test_rate = δ * min(one(conv_rate), conv_rate) / c_conv
        if test_rate <= one(test_rate)
            return true
        end
        # Divergence criteria
        if ((k == maxiters) || (k >= 2 && δ > div_rate * δ_prev))
            return false
        end
        δ_prev = δ
        increment_nf!(integrator.stats, 1)
        isconstcache ? (ratetmp = integrator.f(integrator.u, p, dt + t)) :
            integrator.f(ratetmp, integrator.u, p, dt + t)
    end
    return
end

function nordsieck_rescale!(cache::T, rewind = false) where {T}
    isconstcache = T <: OrdinaryDiffEqConstantCache
    (; z, dts, order) = cache
    eta = rewind ? dts[2] / dts[1] : dts[1] / dts[2]
    factor = eta
    for i in 2:(order + 1)
        if isconstcache
            z[i] = z[i] * factor
        else
            rmul!(z[i], factor)
        end
        factor *= eta
    end
    return nothing
end

function nordsieck_rewind!(cache)
    perform_predict!(cache, true)
    return nordsieck_rescale!(cache, true)
end

function nordsieck_restore_failed_step!(cache, oldest_dt)
    nordsieck_rewind!(cache)
    return nordsieck_restore_dts!(cache, oldest_dt)
end

function nordsieck_restore_dts!(cache, oldest_dt)
    dts = cache.dts
    for i in 1:(length(dts) - 1)
        dts[i] = dts[i + 1]
    end
    dts[end] = oldest_dt
    return nothing
end

function is_nordsieck_change_order(cache::T, n = 0) where {T}
    isconstcache = T <: OrdinaryDiffEqConstantCache
    isvode = (T <: JVODECache || T <: JVODEConstantCache)
    isvode || return false
    return cache.n_wait == 0 + n
end

function nordsieck_decrement_wait!(cache::T) where {T}
    isvode = (T <: JVODECache || T <: JVODEConstantCache)
    isvode || return nothing
    isconstcache = T <: OrdinaryDiffEqConstantCache
    cache.n_wait = max(0, cache.n_wait - 1)
    return nothing
end

function nordsieck_adjust_order!(cache::T, dorder) where {T}
    isconstcache = T <: OrdinaryDiffEqConstantCache
    (; order, dts) = cache
    # WIP: uncomment when finished
    #@inbound begin
    return begin
        # Adams order increase
        if dorder == 1
            if isconstcache
                cache.z[order + 2] = zero(cache.z[order + 2])
            else
                cache.z[order + 2] .= 0
            end
        else
            # Adams order decrease
            # One needs to rescale the Nordsieck vector on an order decrease
            cache.l .= 0
            cache.l[2] = 1
            dt = dts[1]
            hsum = zero(eltype(cache.dts))
            for j in 1:(order - 2)
                hsum += cache.dts[j + 1]
                # TODO: `hscale`?
                ξ = hsum / dt
                for i in (j + 1):-1:1
                    cache.l[i + 1] = cache.l[i + 1] * ξ + cache.l[i]
                end # for i
            end # for j

            for j in 2:(order - 1)
                cache.l[j + 1] = order * cache.l[j] / j
            end
            for j in 3:order
                if isconstcache
                    cache.z[j] = muladd.(-cache.l[j], cache.z[order + 1], cache.z[j])
                else
                    @.. broadcast = false cache.z[j] = muladd(
                        -cache.l[j], cache.z[order + 1],
                        cache.z[j]
                    )
                end
            end # for j
        end # else
    end # @inbound
end

# BDF order change of the Nordsieck vector, following `cvIncreaseBDF` and
# `cvDecreaseBDF` in SUNDIALS CVODE. It runs before the rescale, so `z` is still
# scaled to the previous step `dts[2]`, and CVODE's `tau[j]` is `dts[j + 1]`.
function nordsieck_adjust_order_bdf!(cache::T, dorder) where {T}
    isconstcache = T <: OrdinaryDiffEqConstantCache
    (; z, l, order, dts) = cache
    q = order
    dtscale = dts[2]
    fill!(l, zero(eltype(l)))
    l[3] = one(eltype(l))
    if dorder == 1
        α0 = -one(eltype(l))
        α1 = prod = ξold = one(eltype(l))
        dtsum = dtscale
        for j in 1:(q - 1)
            dtsum += dts[j + 2]
            ξ = dtsum / dtscale
            prod *= ξ
            α0 -= one(eltype(l)) / (j + 1)
            α1 += inv(ξ)
            for i in (j + 2):-1:2
                l[i + 1] = muladd(l[i + 1], ξold, l[i])
            end
            ξold = ξ
        end
        # `z[end]` holds `Δ` of the last step at order `q` (stored by `chooseη!`)
        A1 = (-α0 - α1) / prod
        if isconstcache
            z[q + 2] = A1 * z[end]
            for j in 2:q
                z[j + 1] = muladd.(l[j + 1], z[q + 2], z[j + 1])
            end
        else
            @.. broadcast = false z[q + 2] = A1 * z[end]
            for j in 2:q
                @.. broadcast = false z[j + 1] = muladd(l[j + 1], z[q + 2], z[j + 1])
            end
        end
    elseif q > 2
        dtsum = zero(dtscale)
        for j in 1:(q - 2)
            dtsum += dts[j + 1]
            ξ = dtsum / dtscale
            for i in (j + 2):-1:2
                l[i + 1] = muladd(l[i + 1], ξ, l[i])
            end
        end
        for j in 2:(q - 1)
            if isconstcache
                z[j + 1] = muladd.(-l[j + 1], z[q + 1], z[j + 1])
            else
                @.. broadcast = false z[j + 1] = muladd(-l[j + 1], z[q + 1], z[j + 1])
            end
        end
    end
    return nothing
end

is_bdf(alg) = alg isa JVODE && alg.algorithm === :BDF

# An Adams order change is only applied together with a step size change
function nordsieck_needs_adjust(integrator, cache)
    integrator.dt != cache.dts[2] && return true
    return is_bdf(integrator.alg) && cache.nextorder != cache.order
end

# Maximum order: 5 for BDF (as in CVODE), the length of the Nordsieck array for Adams
nordsieck_qmax(alg) = is_bdf(alg) ? 5 : 12

# `η` is `dtₙ₊₁/dtₙ`
function setη!(integrator, cache::T) where {T}
    if cache.η < get_qsteady_max(integrator)
        cache.η = 1
    else
        # TODO: Not the same with SUNDIALS
        (integrator.iter == 1 || integrator.derivative_discontinuity) &&
            (cache.η = min(1.0e5, cache.η); return nothing)
        cache.η = min(get_qmax(integrator), max(get_qmin(integrator), cache.η))
    end
    return nothing
end

function chooseη!(integrator, cache::T) where {T}
    isconst = T <: OrdinaryDiffEqConstantCache
    (; ηq, η₋₁, η₊₁, order, z, Δ) = cache
    η = max(ηq, η₋₁, η₊₁)
    if η < get_qsteady_max(integrator)
        cache.η = 1
        cache.nextorder = order
        is_bdf(integrator.alg) && return nothing
    end

    if η == ηq
        cache.η = cache.ηq
        cache.nextorder = order
    elseif η == η₋₁
        cache.η = cache.η₋₁
        cache.nextorder = order - 1
    else
        cache.η = cache.η₊₁
        cache.nextorder = order + 1
        if is_bdf(integrator.alg)
            if isconst
                z[end] = Δ
            else
                @.. broadcast = false z[end] = Δ
            end
        end #endif BDF
    end # endif η == ηq
    return nothing
end

function stepsize_η!(integrator, cache, order)
    bias2 = integrator.alg.bias2
    addon = integrator.alg.addon
    L = order + 1
    cache.ηq = inv((bias2 * get_EEst(integrator))^inv(L) + addon)
    return cache.ηq
end

# TODO: Check them
function stepsize_η₊₁!(integrator, cache::T, order) where {T}
    isconstcache = T <: OrdinaryDiffEqConstantCache
    atmp = ratetmp = integrator.uprev  # Initialize for JET
    isconstcache || ((; atmp, ratetmp) = cache)
    (; uprev, t, u) = integrator
    (; z, c_LTE₊₁, dts, c_𝒟) = cache
    bias3 = integrator.alg.bias3
    addon = integrator.alg.addon
    q = order
    cache.η₊₁ = 0
    qmax = nordsieck_qmax(integrator.alg)
    L = q + 1
    if q != qmax
        cache.prev_𝒟 == 0 && return cache.η₊₁
        cquot = (c_𝒟 / cache.prev_𝒟) * (dts[1] / dts[2])^L
        if isconstcache
            atmp = muladd.(-cquot, z[end], cache.Δ)
            atmp = calculate_residuals(
                atmp, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t
            )
        else
            @.. broadcast = false ratetmp = muladd(-cquot, z[end], cache.Δ)
            calculate_residuals!(
                atmp, ratetmp, uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
        end
        dup = integrator.opts.internalnorm(atmp, t) * c_LTE₊₁
        cache.η₊₁ = inv((bias3 * dup)^inv(L + 1) + addon)
    end
    return cache.η₊₁
end

function stepsize_η₋₁!(integrator, cache::T, order) where {T}
    isconstcache = T <: OrdinaryDiffEqConstantCache
    atmp = integrator.uprev  # Initialize for JET
    isconstcache || (atmp = cache.atmp)
    (; uprev, t, u) = integrator
    (; z, c_LTE₋₁) = cache
    bias1 = integrator.alg.bias1
    addon = integrator.alg.addon
    cache.η₋₁ = 0
    if order > 1
        if isconstcache
            atmp = calculate_residuals(
                z[order + 1], uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm,
                t
            )
        else
            calculate_residuals!(
                atmp, z[order + 1], uprev, u, integrator.opts.abstol,
                integrator.opts.reltol, integrator.opts.internalnorm, t
            )
        end
        approx = integrator.opts.internalnorm(atmp, t) * c_LTE₋₁
        cache.η₋₁ = inv((bias1 * approx)^inv(order) + addon)
    end
    return cache.η₋₁
end
