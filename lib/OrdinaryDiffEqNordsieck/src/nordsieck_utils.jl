function nordsieck_adjust!(integrator, cache::T) where {T}
    (; nextorder, order) = cache
    if nextorder != order
        # TODO: optimize?
        nordsieck_adjust_order!(cache, nextorder - order)
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
    return if is_nordsieck_change_order(cache, 1) && cache.order != 12
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
    if OrdinaryDiffEqCore.get_EEst(integrator) > one(OrdinaryDiffEqCore.get_EEst(integrator))
        nordsieck_rewind!(cache)
        cache.n_wait = max(2, cache.n_wait)
        cache.nextorder = order
        cache.η = inv((bias2 * OrdinaryDiffEqCore.get_EEst(integrator))^inv(L) + addon)
        return nothing
    end
    cache.ηq = inv((bias2 * OrdinaryDiffEqCore.get_EEst(integrator))^inv(L) + addon)
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
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
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
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
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

# `η` is `dtₙ₊₁/dtₙ`
function setη!(integrator, cache::T) where {T}
    if cache.η < integrator.alg.qsteady_max
        cache.η = 1
    else
        # TODO: Not the same with SUNDIALS
        (integrator.iter == 1 || integrator.derivative_discontinuity) &&
            (cache.η = min(1.0e5, cache.η); return nothing)
        cache.η = min(integrator.alg.qmax, max(integrator.alg.qmin, cache.η))
    end
    return nothing
end

function chooseη!(integrator, cache::T) where {T}
    isconst = T <: OrdinaryDiffEqConstantCache
    (; ηq, η₋₁, η₊₁, order, z, Δ) = cache
    η = max(ηq, η₋₁, η₊₁)
    if η < integrator.alg.qsteady_max
        cache.η = 1
        cache.nextorder = order
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
        # TODO: BDF
        if integrator.alg.algorithm == :BDF
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
    cache.ηq = inv((bias2 * OrdinaryDiffEqCore.get_EEst(integrator))^inv(L) + addon)
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
    qmax = length(z) - 1
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
