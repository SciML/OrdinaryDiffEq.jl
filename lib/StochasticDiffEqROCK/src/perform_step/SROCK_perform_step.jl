@muladd function perform_step!(integrator, cache::SROCK1ConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    cache.mdeg = Int(floor(sqrt(2 * abs(dt) * integrator.opts.internalnorm(integrator.eigen_est, t)) + 1)) # this is the spectral radius estimate to choose optimal stage
    choose_deg!(integrator, cache)

    mdeg = cache.mdeg
    η = cache.optimal_η
    ω₀ = 1 + (η / (mdeg^2))
    ωSq = (ω₀^2) - 1
    Sqrt_ω = sqrt(ωSq)
    cosh_inv = log(ω₀ + Sqrt_ω)             # arcosh(ω₀)
    ω₁ = (Sqrt_ω * cosh(mdeg * cosh_inv)) / (mdeg * sinh(mdeg * cosh_inv))

    # only used on the Stratonovich paths, but computed unconditionally (cheap scalars)
    α = cosh(mdeg * cosh_inv) / (2 * ω₀ * cosh((mdeg - 1) * cosh_inv))
    γ = 1 / (2 * α)
    β = -γ

    uᵢ₋₂ = copy(uprev)
    k = integrator.f(uprev, p, t)
    Tᵢ₋₂ = oneunit(t)
    Tᵢ₋₁ = oftype(t, ω₀)
    Tᵢ = Tᵢ₋₁
    tᵢ₋₁ = t + dt * (ω₁ / ω₀)
    tᵢ = tᵢ₋₁
    tᵢ₋₂ = t
    gₘ₋₁ = zero(k)
    gₘ₋₂ = zero(k)

    #stage 1
    uᵢ₋₁ = uprev + (dt * ω₁ / ω₀) * k

    for i in 2:mdeg
        Tᵢ = 2 * ω₀ * Tᵢ₋₁ - Tᵢ₋₂
        μ = 2 * ω₁ * (Tᵢ₋₁ / Tᵢ)
        ν = 2 * ω₀ * (Tᵢ₋₁ / Tᵢ)
        κ = (-Tᵢ₋₂ / Tᵢ)
        k = integrator.f(uᵢ₋₁, p, tᵢ₋₁)

        u = dt * μ * k + ν * uᵢ₋₁ + κ * uᵢ₋₂
        if (i > mdeg - 2) &&
                SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Stratonovich
            if i == mdeg - 1
                gₘ₋₂ = integrator.f.g(uᵢ₋₁, p, tᵢ₋₁)
                if W.dW isa Number || !is_diagonal_noise(integrator.sol.prob)
                    u += α * (gₘ₋₂ * W.dW)
                else
                    u .+= α .* gₘ₋₂ .* W.dW
                end
            else
                gₘ₋₁ = integrator.f.g(uᵢ₋₁, p, tᵢ₋₁)
                if W.dW isa Number || !is_diagonal_noise(integrator.sol.prob)
                    u += β * (gₘ₋₂ * W.dW) + γ * (gₘ₋₁ * W.dW)
                else
                    u .+= (β .* gₘ₋₂ .+ γ .* gₘ₋₁) .* W.dW
                end
            end
        elseif (i == mdeg) &&
                SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Ito
            if W.dW isa Number
                gₘ₋₂ = integrator.f.g(uᵢ₋₁, p, tᵢ₋₁)
                uᵢ₋₂ = uᵢ₋₁ + sqrt(abs(dt)) * gₘ₋₂
                gₘ₋₁ = integrator.f.g(uᵢ₋₂, p, tᵢ₋₁)
                u += gₘ₋₂ * W.dW + 1 / (2 * sqrt(abs(dt))) * (gₘ₋₁ - gₘ₋₂) * (W.dW^2 - abs(dt))
            elseif is_diagonal_noise(integrator.sol.prob)
                gₘ₋₂ = integrator.f.g(uᵢ₋₁, p, tᵢ₋₁)
                uᵢ₋₂ .= uᵢ₋₁ .+ sqrt(abs(dt)) .* gₘ₋₂
                gₘ₋₁ = integrator.f.g(uᵢ₋₂, p, tᵢ₋₁)
                u .+= gₘ₋₂ .* W.dW .+
                    (1 / (2 * sqrt(abs(dt)))) .* (gₘ₋₁ .- gₘ₋₂) .* (W.dW .^ 2 .- abs(dt))
            else
                gₘ₋₂ = integrator.f.g(uᵢ₋₁, p, tᵢ₋₁)
                u += gₘ₋₂ * W.dW
            end
        end

        if i < mdeg
            tᵢ = μ * dt + ν * tᵢ₋₁ + κ * tᵢ₋₂
            uᵢ₋₂ = uᵢ₋₁
            uᵢ₋₁ = u
            tᵢ₋₂ = tᵢ₋₁
            tᵢ₋₁ = tᵢ
            Tᵢ₋₂ = Tᵢ₋₁
            Tᵢ₋₁ = Tᵢ
        end
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SROCK1Cache)
    (; uᵢ₋₁, uᵢ₋₂, k, gₘ₋₁, gₘ₋₂) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    ccache.mdeg = Int(floor(sqrt(2 * abs(dt) * integrator.opts.internalnorm(integrator.eigen_est, t)) + 1))   # this is the spectral radius estimate to choose optimal stage
    choose_deg!(integrator, cache)

    mdeg = ccache.mdeg
    η = ccache.optimal_η
    ω₀ = 1 + η / (mdeg^2)
    ωSq = ω₀^2 - 1
    Sqrt_ω = sqrt(ωSq)
    cosh_inv = log(ω₀ + Sqrt_ω)             # arcosh(ω₀)
    ω₁ = (Sqrt_ω * cosh(mdeg * cosh_inv)) / (mdeg * sinh(mdeg * cosh_inv))

    # only used on the Stratonovich paths, but computed unconditionally (cheap scalars)
    α = cosh(mdeg * cosh_inv) / (2 * ω₀ * cosh((mdeg - 1) * cosh_inv))
    γ = 1 / (2 * α)
    β = -γ

    @.. uᵢ₋₂ = uprev
    Tᵢ₋₂ = oneunit(t)
    Tᵢ₋₁ = oftype(t, ω₀)
    Tᵢ = Tᵢ₋₁
    tᵢ₋₁ = t + dt * (ω₁ / ω₀)
    tᵢ = tᵢ₋₁
    tᵢ₋₂ = t

    #stage 1
    #this take advantage of the fact that cache.k === cache.fsalfirst
    #and this has already been done i maxeig!  i.e. integrator.f(fsalfirst, uprev, p, t)
    # integrator.f(k,uprev,p,t)
    @.. uᵢ₋₁ = uprev + (dt * ω₁ / ω₀) * k

    for i in 2:mdeg
        Tᵢ = 2 * ω₀ * Tᵢ₋₁ - Tᵢ₋₂
        μ = 2 * ω₁ * Tᵢ₋₁ / Tᵢ
        ν = 2 * ω₀ * Tᵢ₋₁ / Tᵢ
        κ = - Tᵢ₋₂ / Tᵢ
        integrator.f(k, uᵢ₋₁, p, tᵢ₋₁)
        @.. u = dt * μ * k + ν * uᵢ₋₁ + κ * uᵢ₋₂
        if (i > mdeg - 2) &&
                SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Stratonovich
            if i == mdeg - 1
                integrator.f.g(gₘ₋₂, uᵢ₋₁, p, tᵢ₋₁)
                if W.dW isa Number || is_diagonal_noise(integrator.sol.prob)
                    @.. u += α * gₘ₋₂ * W.dW
                else
                    mul!(k, gₘ₋₂, W.dW)
                    @.. u += α * k
                end
            else
                integrator.f.g(gₘ₋₁, uᵢ₋₁, p, tᵢ₋₁)
                if W.dW isa Number || is_diagonal_noise(integrator.sol.prob)
                    @.. u += (β * gₘ₋₂ + γ * gₘ₋₁) * W.dW
                else
                    mul!(k, gₘ₋₂, W.dW)
                    @.. u += β * k
                    mul!(k, gₘ₋₁, W.dW)
                    @.. u += γ * k
                end
            end
        elseif (i == mdeg) &&
                SciMLBase.alg_interpretation(integrator.alg) ==
                SciMLBase.AlgorithmInterpretation.Ito
            if W.dW isa Number || is_diagonal_noise(integrator.sol.prob)
                integrator.f.g(gₘ₋₂, uᵢ₋₁, p, tᵢ₋₁)
                @.. uᵢ₋₂ = uᵢ₋₁ + sqrt(abs(dt)) * gₘ₋₂
                integrator.f.g(gₘ₋₁, uᵢ₋₂, p, tᵢ₋₁)
                @.. u += gₘ₋₂ * W.dW + 1 / (2 * sqrt(abs(dt))) * (gₘ₋₁ - gₘ₋₂) * (W.dW^2 - abs(dt))
            else
                integrator.f.g(gₘ₋₂, uᵢ₋₁, p, tᵢ₋₁)
                mul!(uᵢ₋₁, gₘ₋₂, W.dW)
                @.. u += uᵢ₋₁
            end
        end

        if i < mdeg
            tᵢ = dt * μ + ν * tᵢ₋₁ + κ * tᵢ₋₂
            @.. uᵢ₋₂ = uᵢ₋₁
            @.. uᵢ₋₁ = u
            tᵢ₋₂ = tᵢ₋₁
            tᵢ₋₁ = tᵢ
            Tᵢ₋₂ = Tᵢ₋₁
            Tᵢ₋₁ = Tᵢ
        end
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SROCK2ConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    (; recf, recf2, mα, mσ, mτ) = cache

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    cache.mdeg = Int(floor(sqrt((2 * abs(dt) * integrator.opts.internalnorm(integrator.eigen_est, t) + 1.5) / 0.811) + 1))
    cache.mdeg = max(3, min(cache.mdeg, 200)) - 2
    choose_deg!(integrator, cache)

    mdeg = cache.mdeg
    start = cache.start
    deg_index = cache.deg_index
    α = mα[deg_index]
    σ = (1 - α) * 1 // 2 + α * mσ[deg_index]
    τ = 1 // 2 * ((1 - α)^2) + 2 * α * (1 - α) * mσ[deg_index] +
        (α^2) * (mσ[deg_index] * (mσ[deg_index] + mτ[deg_index]))

    sqrt_dt = sqrt(abs(dt))

    μ = recf[start]  # here κ = 0
    tᵢ = t + α * dt * μ
    tᵢ₋₁ = tᵢ
    tᵢ₋₂ = t

    # stage 1
    uᵢ₋₂ = uprev
    uᵢ = integrator.f(uprev, p, t)
    uᵢ₋₁ = uprev + α * dt * μ * uᵢ

    # stages 2 upto s-2
    for i in 2:mdeg
        μ, κ = recf[start + 2 * (i - 2) + 1], recf[start + 2 * (i - 2) + 2]
        ν = 1 + κ
        uᵢ = integrator.f(uᵢ₋₁, p, tᵢ₋₁)
        uᵢ = α * dt * μ * uᵢ + ν * uᵢ₋₁ - κ * uᵢ₋₂
        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = uᵢ
        tᵢ = α * dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂
        tᵢ₋₂ = tᵢ₋₁
        tᵢ₋₁ = tᵢ
    end

    #2 stage-finishing procedure
    #stage s-1
    μ, κ = recf2[(deg_index - 1) * 4 + 1], recf2[(deg_index - 1) * 4 + 2]
    ν = 1 + κ
    uᵢ = integrator.f(uᵢ₋₁, p, tᵢ₋₁)

    tₓ = tᵢ₋₁ + 2 * dt * τ                    # So that we don't have to calculate f(uₛ₋₂) again
    uₓ = uᵢ₋₁ + 2 * dt * τ * uᵢ                 # uₓ and tₓ represent u_star
    u = uᵢ₋₁ + (2 * σ - 1 // 2) * dt * uᵢ

    uᵢ = α * dt * μ * uᵢ + ν * uᵢ₋₁ - κ * uᵢ₋₂
    tᵢ = α * dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂
    tᵢ₋₂ = tᵢ₋₁
    tᵢ₋₁ = tᵢ
    uᵢ₋₂ = uᵢ₋₁
    uᵢ₋₁ = uᵢ

    #stage s
    μ, κ = recf2[(deg_index - 1) * 4 + 3], recf2[(deg_index - 1) * 4 + 4]
    ν = 1 + κ
    uᵢ = integrator.f(uᵢ₋₁, p, tᵢ₋₁)
    uᵢ = α * dt * μ * uᵢ + ν * uᵢ₋₁ - κ * uᵢ₋₂
    tᵢ = α * dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂

    # Now uᵢ₋₂ = uₛ₋₂, uᵢ₋₁ = uₛ₋₁, uᵢ = uₛ
    # Similarly tᵢ₋₂ = tₛ₋₂, tᵢ₋₁ = tₛ₋₁, tᵢ = tₛ

    if (W.dW isa Number) || (length(W.dW) == 1) || is_diagonal_noise(integrator.sol.prob)
        Gₛ = integrator.f.g(uᵢ₋₁, p, tᵢ₋₁)
        u += Gₛ .* W.dW
        Gₛ = integrator.f.g(uᵢ, p, tᵢ)
        uₓ += Gₛ .* W.dW

        uₓ = integrator.f(uₓ, p, tₓ)
        u += (1 // 2 * dt) * uₓ
        uₓ = 1 // 2 .* Gₛ .* (W.dW .^ 2 .- abs(dt))
        uᵢ₋₂ = uᵢ + uₓ
        Gₛ₁ = integrator.f.g(uᵢ₋₂, p, tᵢ)
        u += (1 // 2) * Gₛ₁
        uᵢ₋₂ = uᵢ - uₓ
        Gₛ₁ = integrator.f.g(uᵢ₋₂, p, tᵢ)
        u -= (1 // 2) * Gₛ₁

        uₓ = Gₛ * sqrt_dt
        uᵢ₋₂ = uᵢ + uₓ
        Gₛ₁ = integrator.f.g(uᵢ₋₂, p, tᵢ)
        u += 1 // 4 .* W.dW .* (Gₛ₁ .- 2 .* Gₛ)
        uᵢ₋₂ = uᵢ - uₓ
        Gₛ₁ = integrator.f.g(uᵢ₋₂, p, tᵢ)
        u += 1 // 4 .* W.dW .* Gₛ₁
    else
        vec_χ = 2 .* floor.(false .* W.dW .+ 1 // 2 .+ oftype(W.dW, rand(W.rng, length(W.dW)))) .- true

        Gₛ = integrator.f.g(uᵢ₋₁, p, tᵢ₋₁)
        u += Gₛ * W.dW

        Gₛ = integrator.f.g(uᵢ, p, tᵢ)
        uₓ += Gₛ * W.dW

        uₓ = integrator.f(uₓ, p, tₓ)
        u += (1 // 2) * dt * uₓ
        for i in 1:length(W.dW)
            WikJ = W.dW[i]
            WikJ2 = vec_χ[i]
            WikRange = 1 // 2 .* (W.dW .* WikJ .- (1:length(W.dW) .== i) .* abs(dt)) #.- (1:length(W.dW) .> i) .* dt .* vec_χ .+ (1:length(W.dW) .< i) .* dt .* WikJ2)
            uₓ = Gₛ * WikRange
            WikRange = 1 // 2 .* (1:length(W.dW) .== i)
            uᵢ₋₂ = uᵢ + uₓ
            Gₛ₁ = integrator.f.g(uᵢ₋₂, p, tᵢ)
            u += (Gₛ₁ * WikRange)
            uᵢ₋₂ = uᵢ - uₓ
            Gₛ₁ = integrator.f.g(uᵢ₋₂, p, tᵢ)
            u -= (Gₛ₁ * WikRange)
        end

        uₓ = sqrt_dt * (Gₛ * vec_χ)
        uᵢ₋₂ = uᵢ + uₓ
        Gₛ₁ = integrator.f.g(uᵢ₋₂, p, tᵢ)
        u += 1 // 4 * (Gₛ₁ * W.dW) - 1 // 2 * (Gₛ * W.dW)

        uᵢ₋₂ = uᵢ - uₓ
        Gₛ₁ = integrator.f.g(uᵢ₋₂, p, tᵢ)
        u += 1 // 4 * (Gₛ₁ * W.dW)
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SROCK2Cache)
    (; uᵢ, uₓ, uᵢ₋₁, uᵢ₋₂, k, Gₛ, Gₛ₁, vec_χ, WikRange) = cache

    (; t, dt, uprev, u, W, p, f) = integrator

    (; recf, recf2, mα, mσ, mτ) = cache.constantcache
    ccache = cache.constantcache
    gen_prob = !(
        (is_diagonal_noise(integrator.sol.prob)) || (W.dW isa Number) ||
            (length(W.dW) == 1)
    )

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    ccache.mdeg = Int(
        floor(
            sqrt(
                (
                    2 * abs(dt) * integrator.opts.internalnorm(
                        integrator.eigen_est, t
                    ) + 1.5
                ) / 0.811
            ) + 1
        )
    )
    ccache.mdeg = max(3, min(ccache.mdeg, 200)) - 2
    choose_deg!(integrator, cache)

    mdeg = ccache.mdeg
    start = ccache.start
    deg_index = ccache.deg_index
    α = mα[deg_index]
    σ = (1 - α) * 1 // 2 + α * mσ[deg_index]
    τ = 1 // 2 * ((1 - α)^2) + 2 * α * (1 - α) * mσ[deg_index] +
        (α^2) * (mσ[deg_index] * (mσ[deg_index] + mτ[deg_index]))

    sqrt_dt = sqrt(abs(dt))
    if gen_prob
        vec_χ .= 1 // 2 .+ oftype(W.dW, rand(W.rng, length(W.dW)))
        @.. vec_χ = 2 * floor(vec_χ) - 1
    end

    μ = recf[start]  # here κ = 0
    tᵢ = t + α * dt * μ
    tᵢ₋₁ = tᵢ
    tᵢ₋₂ = t

    # stage 1
    @.. uᵢ₋₂ = uprev
    #this take advantage of the fact that cache.k === cache.fsalfirst
    #and this has already been done i maxeig!  i.e. integrator.f(fsalfirst, uprev, p, t)
    # integrator.f(k,uprev,p,t)
    @.. uᵢ₋₁ = uprev + α * dt * μ * k

    # stages 2 upto s-2
    for i in 2:mdeg
        μ, κ = recf[start + 2 * (i - 2) + 1], recf[start + 2 * (i - 2) + 2]
        ν = 1 + κ
        integrator.f(k, uᵢ₋₁, p, t)
        @.. uᵢ = α * dt * μ * k + ν * uᵢ₋₁ - κ * uᵢ₋₂
        @.. uᵢ₋₂ = uᵢ₋₁
        @.. uᵢ₋₁ = uᵢ
        tᵢ = α * dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂
        tᵢ₋₂ = tᵢ₋₁
        tᵢ₋₁ = tᵢ
    end

    #2 stage-finishing procedure
    #stage s-1
    μ, κ = recf2[(deg_index - 1) * 4 + 1], recf2[(deg_index - 1) * 4 + 2]
    ν = 1 + κ
    integrator.f(k, uᵢ₋₁, p, tᵢ₋₁)

    tₓ = tᵢ₋₁ + 2 * dt * τ                    # So that we don't have to calculate f(uₛ₋₂) again
    @.. uₓ = uᵢ₋₁ + 2 * dt * τ * k                 # uₓ and tₓ represent u_star
    @.. u = uᵢ₋₁ + (2 * σ - 1 // 2) * dt * k

    @.. uᵢ = α * dt * μ * k + ν * uᵢ₋₁ - κ * uᵢ₋₂
    tᵢ = α * dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂
    tᵢ₋₂ = tᵢ₋₁
    tᵢ₋₁ = tᵢ
    @.. uᵢ₋₂ = uᵢ₋₁
    @.. uᵢ₋₁ = uᵢ

    #stage s
    μ, κ = recf2[(deg_index - 1) * 4 + 3], recf2[(deg_index - 1) * 4 + 4]
    ν = 1 + κ
    integrator.f(k, uᵢ₋₁, p, tᵢ₋₁)
    @.. uᵢ = α * dt * μ * k + ν * uᵢ₋₁ - κ * uᵢ₋₂
    tᵢ = α * dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂

    # Now uᵢ₋₂ = uₛ₋₂, uᵢ₋₁ = uₛ₋₁, uᵢ = uₛ
    # Similarly tᵢ₋₂ = tₛ₋₂, tᵢ₋₁ = tₛ₋₁, tᵢ = tₛ

    if (W.dW isa Number) || (length(W.dW) == 1) || is_diagonal_noise(integrator.sol.prob)
        integrator.f.g(Gₛ, uᵢ₋₁, p, tᵢ₋₁)
        @.. u += Gₛ * W.dW
        integrator.f.g(Gₛ, uᵢ, p, tᵢ)
        @.. uₓ += Gₛ * W.dW
        integrator.f(k, uₓ, p, tₓ)
        @.. u += (1 // 2) * dt * k

        @.. uₓ = Gₛ * ((W.dW^2 - abs(dt)) / 2)
        @.. uᵢ₋₂ = uᵢ + uₓ
        integrator.f.g(Gₛ₁, uᵢ₋₂, p, tᵢ)
        @.. u += (1 // 2) * Gₛ₁
        @.. uᵢ₋₂ = uᵢ - uₓ
        integrator.f.g(Gₛ₁, uᵢ₋₂, p, tᵢ)
        @.. u -= (1 // 2) * Gₛ₁

        @.. uₓ = sqrt_dt * Gₛ
        @.. uᵢ₋₂ = uᵢ + uₓ
        integrator.f.g(Gₛ₁, uᵢ₋₂, p, tᵢ)
        @.. u += (1 // 4 * W.dW) * (Gₛ₁ - 2 * Gₛ)
        @.. uᵢ₋₂ = uᵢ - uₓ
        integrator.f.g(Gₛ₁, uᵢ₋₂, p, tᵢ)
        @.. u += (1 // 4 * W.dW) * Gₛ₁
    else
        integrator.f.g(Gₛ, uᵢ₋₁, p, tᵢ₋₁)
        mul!(uᵢ₋₂, Gₛ, W.dW)
        @.. u += uᵢ₋₂

        integrator.f.g(Gₛ, uᵢ, p, tᵢ)
        mul!(uᵢ₋₁, Gₛ, W.dW)
        @.. uₓ += uᵢ₋₁

        integrator.f(k, uₓ, p, tₓ)
        @.. u += (1 // 2) * dt * k

        for i in 1:length(W.dW)
            WikJ = W.dW[i]
            WikJ2 = vec_χ[i]
            dwrange = 1:length(W.dW)
            abs_dt = abs(dt)
            @.. WikRange = 1 // 2 * (W.dW * WikJ - (dwrange == i) * abs_dt) #+ (dwrange < i) * dt * WikJ2 - (dwrange > i) * dt * vec_χ)
            mul!(uₓ, Gₛ, WikRange)
            @.. uᵢ₋₂ = uᵢ + uₓ
            @.. WikRange = 1 // 2 * (dwrange == i)
            integrator.f.g(Gₛ₁, uᵢ₋₂, p, tᵢ)
            mul!(uᵢ₋₂, Gₛ₁, WikRange)
            @.. u += uᵢ₋₂
            @.. uᵢ₋₂ = uᵢ - uₓ
            integrator.f.g(Gₛ₁, uᵢ₋₂, p, tᵢ)
            mul!(uᵢ₋₂, Gₛ₁, WikRange)
            @.. u -= uᵢ₋₂
        end
        @.. WikRange = vec_χ * sqrt_dt
        mul!(uₓ, Gₛ, WikRange)

        @.. uᵢ₋₂ = uᵢ + uₓ
        integrator.f.g(Gₛ₁, uᵢ₋₂, p, tᵢ)
        mul!(uᵢ₋₂, Gₛ₁, W.dW)
        # mul!(uᵢ₋₁,Gₛ,W.dW) This is already been calculated
        @.. u += 1 // 4 * (uᵢ₋₂ - 2 * uᵢ₋₁)

        @.. uᵢ₋₂ = uᵢ - uₓ
        integrator.f.g(Gₛ₁, uᵢ₋₂, p, tᵢ)
        mul!(uᵢ₋₂, Gₛ₁, W.dW)
        @.. u += 1 // 4 * uᵢ₋₂
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SROCKEMConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    if integrator.alg.strong_order_1
        cache.mdeg = Int(floor(sqrt(abs(dt) * integrator.opts.internalnorm(integrator.eigen_est, t) / 0.19) + 1))
    else
        cache.mdeg = Int(floor(sqrt(abs(dt) * integrator.opts.internalnorm(integrator.eigen_est, t) / 0.33) + 1))
    end
    cache.mdeg = max(3, min(cache.mdeg, 200))
    choose_deg!(integrator, cache)

    mdeg = cache.mdeg
    η = cache.optimal_η
    ω₀ = oneunit(t) + (η / (mdeg^2))
    ωSq = ω₀^2 - oneunit(t)
    Sqrt_ω = sqrt(ωSq)
    cosh_inv = log(ω₀ + Sqrt_ω)             # arcosh(ω₀)
    ω₁ = (Sqrt_ω * cosh(mdeg * cosh_inv)) / (mdeg * sinh(mdeg * cosh_inv))

    uᵢ₋₂ = copy(uprev)
    k = integrator.f(uprev, p, t)
    Tᵢ₋₂ = oneunit(t)
    Tᵢ₋₁ = oftype(t, ω₀)
    Tᵢ = Tᵢ₋₁
    tᵢ₋₁ = t + dt * (ω₁ / ω₀)
    tᵢ = tᵢ₋₁
    tᵢ₋₂ = t

    #stage 1
    uᵢ₋₁ = uprev + dt * (ω₁ / ω₀) * k

    for i in 2:mdeg
        Tᵢ = 2 * ω₀ * Tᵢ₋₁ - Tᵢ₋₂
        μ = 2 * ω₁ * (Tᵢ₋₁ / Tᵢ)
        ν = 2 * ω₀ * (Tᵢ₋₁ / Tᵢ)
        κ = -Tᵢ₋₂ / Tᵢ
        k = integrator.f(uᵢ₋₁, p, tᵢ₋₁)

        u = @. dt * μ * k + ν * uᵢ₋₁ + κ * uᵢ₋₂
        tᵢ = @. μ * dt + ν * tᵢ₋₁ + κ * tᵢ₋₂

        if i < mdeg
            uᵢ₋₂ = uᵢ₋₁
            uᵢ₋₁ = u
            tᵢ₋₂ = tᵢ₋₁
            tᵢ₋₁ = tᵢ
            Tᵢ₋₂ = Tᵢ₋₁
            Tᵢ₋₁ = Tᵢ
        end
    end

    Gₛ = integrator.f.g(u, p, tᵢ)
    if (W.dW isa Number) || (length(W.dW) == 1) || is_diagonal_noise(integrator.sol.prob)
        u += Gₛ .* W.dW
    else
        u += Gₛ * W.dW
    end

    if integrator.alg.strong_order_1
        if (W.dW isa Number) || (length(W.dW) == 1) ||
                (is_diagonal_noise(integrator.sol.prob))
            uᵢ₋₂ = @. 1 // 2 * Gₛ * (W.dW^2 - abs(dt))
            tmp = @. u + uᵢ₋₂
            Gₛ = integrator.f.g(tmp, p, tᵢ)
            uᵢ₋₁ = @. 1 // 2 * Gₛ
            tmp = @. u - uᵢ₋₂
            Gₛ = integrator.f.g(tmp, p, tᵢ)
            uᵢ₋₁ = @. uᵢ₋₁ - 1 // 2 * Gₛ
            u = @. u + uᵢ₋₁
        else
            for i in 1:length(W.dW)
                WikJ = W.dW[i]
                WikRange = 1 // 2 .* (W.dW .* WikJ - (1:length(W.dW) .== i) .* abs(dt))
                uᵢ₋₂ = Gₛ * WikRange
                WikRange = 1 // 2 .* (1:length(W.dW) .== i)
                tmp = u + uᵢ₋₂
                Gₛ₁ = integrator.f.g(tmp, p, tᵢ)
                uᵢ₋₁ = Gₛ₁ * WikRange
                tmp = u - uᵢ₋₂
                Gₛ₁ = integrator.f.g(tmp, p, tᵢ)
                uᵢ₋₁ -= Gₛ₁ * WikRange
                u += uᵢ₋₁
            end
        end
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SROCKEMCache)
    (; uᵢ₋₁, uᵢ₋₂, tmp, k, Gₛ, Gₛ₁, WikRange) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    if integrator.alg.strong_order_1
        ccache.mdeg = Int(floor(sqrt(abs(dt) * integrator.opts.internalnorm(integrator.eigen_est, t) / 0.19) + 1))
    else
        ccache.mdeg = Int(floor(sqrt(abs(dt) * integrator.opts.internalnorm(integrator.eigen_est, t) / 0.33) + 1))
    end
    ccache.mdeg = max(3, min(ccache.mdeg, 200))
    choose_deg!(integrator, cache)

    mdeg = ccache.mdeg
    η = ccache.optimal_η
    ω₀ = oneunit(t) + η / (mdeg^2)
    ωSq = ω₀^2 - oneunit(t)
    Sqrt_ω = sqrt(ωSq)
    cosh_inv = log(ω₀ + Sqrt_ω)             # arcosh(ω₀)
    ω₁ = (Sqrt_ω * cosh(mdeg * cosh_inv)) / (mdeg * sinh(mdeg * cosh_inv))

    @.. uᵢ₋₂ = uprev
    integrator.f(k, uprev, p, t)
    Tᵢ₋₂ = oneunit(t)
    Tᵢ₋₁ = convert(eltype(u), ω₀)
    Tᵢ = Tᵢ₋₁
    tᵢ₋₁ = t + dt * (ω₁ / ω₀)
    tᵢ = tᵢ₋₁
    tᵢ₋₂ = t

    #stage 1
    @.. uᵢ₋₁ = uprev + (dt * ω₁ / ω₀) * k

    for i in 2:mdeg
        Tᵢ = 2 * ω₀ * Tᵢ₋₁ - Tᵢ₋₂
        μ = 2 * ω₁ * Tᵢ₋₁ / Tᵢ
        ν = 2 * ω₀ * Tᵢ₋₁ / Tᵢ
        κ = - Tᵢ₋₂ / Tᵢ
        integrator.f(k, uᵢ₋₁, p, tᵢ₋₁)
        @.. u = dt * μ * k + ν * uᵢ₋₁ + κ * uᵢ₋₂
        tᵢ = dt * μ + ν * tᵢ₋₁ + κ * tᵢ₋₂

        if i < mdeg
            @.. uᵢ₋₂ = uᵢ₋₁
            @.. uᵢ₋₁ = u
            tᵢ₋₂ = tᵢ₋₁
            tᵢ₋₁ = tᵢ
            Tᵢ₋₂ = Tᵢ₋₁
            Tᵢ₋₁ = Tᵢ
        end
    end

    integrator.f.g(Gₛ, u, p, tᵢ)
    if (W.dW isa Number) || (length(W.dW) == 1) || is_diagonal_noise(integrator.sol.prob)
        @.. u += Gₛ * W.dW
    else
        mul!(uᵢ₋₁, Gₛ, W.dW)
        u += uᵢ₋₁
    end

    if integrator.alg.strong_order_1
        if (W.dW isa Number) || (length(W.dW) == 1) ||
                (is_diagonal_noise(integrator.sol.prob))
            @.. uᵢ₋₂ = 1 // 2 * Gₛ * (W.dW^2 - abs(dt))
            @.. tmp = u + uᵢ₋₂
            integrator.f.g(Gₛ, tmp, p, tᵢ)
            @.. uᵢ₋₁ = 1 // 2 * Gₛ
            @.. tmp = u - uᵢ₋₂
            integrator.f.g(Gₛ, tmp, p, tᵢ)
            @.. uᵢ₋₁ -= 1 // 2 * Gₛ
            u += uᵢ₋₁
        else
            for i in 1:length(W.dW)
                WikJ = W.dW[i]
                dwrange = 1:length(W.dW)
                abs_dt = abs(dt)
                @.. WikRange = 1 // 2 * (W.dW * WikJ - (dwrange == i) * abs_dt)
                mul!(uᵢ₋₂, Gₛ, WikRange)
                @.. WikRange = 1 // 2 * (dwrange == i)
                @.. tmp = u + uᵢ₋₂
                integrator.f.g(Gₛ₁, tmp, p, tᵢ)
                mul!(uᵢ₋₁, Gₛ₁, WikRange)
                @.. tmp = u - uᵢ₋₂
                integrator.f.g(Gₛ₁, tmp, p, tᵢ)
                @.. u += uᵢ₋₁
                mul!(uᵢ₋₁, Gₛ₁, WikRange)
                @.. u -= uᵢ₋₁
            end
        end
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SKSROCKConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    η = oftype(t, 0.05)
    mdeg = Int(
        floor(
            sqrt(
                (
                    abs(dt) * integrator.opts.internalnorm(integrator.eigen_est, t) +
                        1.5
                ) / (2 - η * 4 / 3)
            ) + 1
        )
    )
    mdeg = max(3, min(mdeg, 200))

    ω₀ = 1 + (η / (mdeg^2))
    ωSq = (ω₀^2) - 1
    Sqrt_ω = sqrt(ωSq)
    cosh_inv = log(ω₀ + Sqrt_ω)             # arcosh(ω₀)
    ω₁ = (Sqrt_ω * cosh(mdeg * cosh_inv)) / (mdeg * sinh(mdeg * cosh_inv))

    μ, ν, κ = ω₁ / ω₀, mdeg * ω₁ / 2, mdeg * ω₁ / ω₀
    Tᵢ₋₂ = oneunit(t)
    Tᵢ₋₁ = convert(eltype(u), ω₀)
    Tᵢ = Tᵢ₋₁
    tᵢ₋₁ = t + dt * (ω₁ / ω₀)
    tᵢ = tᵢ₋₁
    tᵢ₋₂ = t

    #stage 1
    Gₛ = integrator.f.g(uprev, p, t)
    if (W.dW isa Number) || !is_diagonal_noise(integrator.sol.prob)
        u = Gₛ * W.dW
    else
        u = Gₛ .* W.dW
    end

    if integrator.alg.post_processing
        uᵢ₋₁ = uprev + ν * u
        uᵢ₋₂ = integrator.f(uprev, p, t)
        uᵢ₋₁ = integrator.f(uᵢ₋₁, p, t)
        uᵢ₋₁ = uprev + (μ * dt) * uᵢ₋₁ + κ * u + cache.mα[mdeg - 1] * dt * (uᵢ₋₁ - 2 * uᵢ₋₂)
        uᵢ₋₂ = uprev - ν * u
        uᵢ₋₂ = integrator.f(uᵢ₋₂, p, t)
        uᵢ₋₁ += (cache.mα[mdeg - 1] * dt) * uᵢ₋₂
    else
        uᵢ₋₁ = uprev + ν * u
        uᵢ₋₂ = integrator.f(uᵢ₋₁, p, t)
        uᵢ₋₁ = uprev + (μ * dt) * uᵢ₋₂ + κ * u
    end

    uᵢ₋₂ = uprev

    for i in 2:mdeg
        Tᵢ = 2 * ω₀ * Tᵢ₋₁ - Tᵢ₋₂
        μ = 2 * ω₁ * (Tᵢ₋₁ / Tᵢ)
        ν = 2 * ω₀ * (Tᵢ₋₁ / Tᵢ)
        κ = (-Tᵢ₋₂ / Tᵢ)
        u = integrator.f(uᵢ₋₁, p, tᵢ₋₁)

        u = dt * μ * u + ν * uᵢ₋₁ + κ * uᵢ₋₂
        tᵢ = μ * dt + ν * tᵢ₋₁ + κ * tᵢ₋₂

        if i < mdeg
            uᵢ₋₂ = uᵢ₋₁
            uᵢ₋₁ = u
            tᵢ₋₂ = tᵢ₋₁
            tᵢ₋₁ = tᵢ
            Tᵢ₋₂ = Tᵢ₋₁
            Tᵢ₋₁ = Tᵢ
        end
    end
    if integrator.alg.post_processing && (t + dt >= integrator.sol.prob.tspan[2])
        Gₛ = integrator.f.g(u, p, tᵢ)
        if (W.dW isa Number) || is_diagonal_noise(integrator.sol.prob)
            uᵢ₋₁ = Gₛ
        else
            WikRange = 1 .* (1:length(W.dW) .== 1)
            uᵢ₋₁ = Gₛ * WikRange
        end
        winc = rand() * 6
        if winc < 1
            # cache IS the constant cache here (cf. the SKSROCKCache method, which
            # accesses this through ccache = cache.constantcache); `ccache` was an
            # undefined-variable bug on this path
            u -= (sqrt(3 * dt) * cache.mc[mdeg - 1]) * uᵢ₋₁
        elseif winc < 2
            u += (sqrt(3 * dt) * cache.mc[mdeg - 1]) * uᵢ₋₁
        end
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SKSROCKCache)
    (; uᵢ₋₁, uᵢ₋₂, k, Gₛ, WikRange) = cache
    (; t, dt, uprev, u, W, p, f) = integrator

    ccache = cache.constantcache
    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    η = oftype(t, 0.05)
    mdeg = Int(
        floor(
            sqrt(
                (
                    abs(dt) * integrator.opts.internalnorm(integrator.eigen_est, t) +
                        1.5
                ) / (2 - η * 4 / 3)
            ) + 1
        )
    )
    mdeg = max(3, min(mdeg, 200))

    ω₀ = 1 + (η / (mdeg^2))
    ωSq = (ω₀^2) - 1
    Sqrt_ω = sqrt(ωSq)
    cosh_inv = log(ω₀ + Sqrt_ω)             # arcosh(ω₀)
    ω₁ = (Sqrt_ω * cosh(mdeg * cosh_inv)) / (mdeg * sinh(mdeg * cosh_inv))

    μ, ν, κ = ω₁ / ω₀, mdeg * ω₁ / 2, mdeg * ω₁ / ω₀
    Tᵢ₋₂ = oneunit(t)
    Tᵢ₋₁ = convert(eltype(u), ω₀)
    Tᵢ = Tᵢ₋₁
    tᵢ₋₁ = t + dt * (ω₁ / ω₀)
    tᵢ = tᵢ₋₁
    tᵢ₋₂ = t

    #stage 1
    integrator.f.g(Gₛ, uprev, p, t)
    if (W.dW isa Number) || is_diagonal_noise(integrator.sol.prob)
        @.. u = Gₛ * W.dW
    else
        mul!(u, Gₛ, W.dW)
    end

    if integrator.alg.post_processing
        @.. uᵢ₋₂ = uprev + ν * u
        integrator.f(k, uᵢ₋₂, p, t)
        @.. uᵢ₋₁ = uprev + (μ * dt) * k + κ * u + (ccache.mα[mdeg - 1] * dt) * k
        integrator.f(k, uprev, p, t)
        @.. uᵢ₋₁ -= (ccache.mα[mdeg - 1] * dt * 2) * k
        @.. uᵢ₋₂ = uprev - ν * u
        integrator.f(k, uᵢ₋₂, p, t)
        @.. uᵢ₋₁ += (ccache.mα[mdeg - 1] * dt) * k
    else
        @.. uᵢ₋₁ = uprev + ν * u
        integrator.f(k, uᵢ₋₁, p, t)
        @.. uᵢ₋₁ = uprev + (μ * dt) * k + κ * u
    end

    @.. uᵢ₋₂ = uprev

    for i in 2:mdeg
        Tᵢ = 2 * ω₀ * Tᵢ₋₁ - Tᵢ₋₂
        μ = 2 * ω₁ * (Tᵢ₋₁ / Tᵢ)
        ν = 2 * ω₀ * (Tᵢ₋₁ / Tᵢ)
        κ = (-Tᵢ₋₂ / Tᵢ)
        integrator.f(k, uᵢ₋₁, p, tᵢ₋₁)

        @.. u = dt * μ * k + ν * uᵢ₋₁ + κ * uᵢ₋₂
        tᵢ = μ * dt + ν * tᵢ₋₁ + κ * tᵢ₋₂

        if i < mdeg
            @.. uᵢ₋₂ = uᵢ₋₁
            @.. uᵢ₋₁ = u
            tᵢ₋₂ = tᵢ₋₁
            tᵢ₋₁ = tᵢ
            Tᵢ₋₂ = Tᵢ₋₁
            Tᵢ₋₁ = Tᵢ
        end
    end

    if integrator.alg.post_processing && (t + dt >= integrator.sol.prob.tspan[2])
        integrator.f.g(Gₛ, u, p, tᵢ)
        # println(Gₛ/length(W.dW))
        if (W.dW isa Number) || is_diagonal_noise(integrator.sol.prob)
            @.. uᵢ₋₁ = Gₛ
        else
            WikRange .= 1 .* (1:length(W.dW) .== 1)
            mul!(uᵢ₋₁, Gₛ, WikRange)
        end
        winc = rand() * 6
        if winc < 1
            @.. u -= (sqrt(3 * dt) * ccache.mc[mdeg - 1]) * uᵢ₋₁
        elseif winc < 2
            @.. u += (sqrt(3 * dt) * ccache.mc[mdeg - 1]) * uᵢ₋₁
        end
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::TangXiaoSROCK2ConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    (; recf, recf2, mα, mσ, mτ, mn̂, c1, c2) = cache

    n̂ = mn̂[integrator.alg.version_num]

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    (integrator.alg.version_num <= 2) &&
        (
        cache.mdeg = Int(
            floor(
                sqrt(
                    (
                        abs(dt) * integrator.opts.internalnorm(
                            integrator.eigen_est, t
                        ) + 1.5
                    ) / 0.811
                ) + 1
            )
        )
    )
    (integrator.alg.version_num > 2) &&
        (
        cache.mdeg = Int(
            floor(
                sqrt(
                    (
                        abs(dt) * integrator.opts.internalnorm(
                            integrator.eigen_est, t
                        ) + 1.5
                    ) / 0.611
                ) + 1
            )
        )
    )

    cache.mdeg = max(4, min(cache.mdeg, 200)) - 2
    choose_deg!(integrator, cache)

    mdeg = cache.mdeg
    start = cache.start
    deg_index = cache.deg_index
    α = mα[integrator.alg.version_num]
    σ = (1 - α) * 1 // 2 + α * mσ[deg_index]
    τ = 1 // 2 * ((1 - α)^2) + 2 * α * (1 - α) * mσ[deg_index] +
        (α^2) * (mσ[deg_index] * (mσ[deg_index] + mτ[deg_index]))

    η₁ = (rand() < 1 // 2) ? -1 : 1
    η₂ = (rand() < 1 // 2) ? -1 : 1
    sqrt_dt = sqrt(abs(dt))

    Û₁ = zero(u)
    Û₂ = zero(u)
    t̂₁ = t̂₂ = zero(t)
    tᵢ = tᵢ₋₁ = tᵢ₋₂ = tₓ = t
    uᵢ = uᵢ₋₁ = uₓ = uᵢ₋₂ = uprev

    for i in 0:(mdeg + 1)
        if i == 1
            μ = recf[start]
            tᵢ = tᵢ₋₁ = t + α * dt * μ

            uᵢ = integrator.f(uprev, p, t)
            uᵢ₋₁ = uprev + α * dt * μ * uᵢ
        elseif i > 1 && i <= mdeg
            μ, ν,
                κ = recf[start + 2 * (i - 2) + 1], 1 + recf[start + 2 * (i - 2) + 2],
                recf[start + 2 * (i - 2) + 2]

            uᵢ = integrator.f(uᵢ₋₁, p, tᵢ₋₁)
            uᵢ = α * dt * μ * uᵢ + ν * uᵢ₋₁ - κ * uᵢ₋₂
            tᵢ = α * dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂
        elseif i == mdeg + 1
            μ, ν,
                κ = recf2[(deg_index - 1) * 4 + 1], 1 + recf2[(deg_index - 1) * 4 + 2],
                recf2[(deg_index - 1) * 4 + 2]

            uᵢ = integrator.f(uᵢ₋₁, p, tᵢ₋₁)

            tₓ = tᵢ₋₁ + 2 * τ * dt
            uₓ = uᵢ₋₁ + (2 * τ * dt) * uᵢ
            u = uᵢ₋₁ + (2 * σ - 1 // 2) * dt * uᵢ

            uᵢ = α * dt * μ * uᵢ + ν * uᵢ₋₁ - κ * uᵢ₋₂
            tᵢ = α * dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂
        end

        j = i - mdeg - 1 + n̂
        if j > 0
            j += cache.start_mcs - 1
            if i == 0
                Û₁ += c1[j] * uprev
                t̂₁ += c1[j] * t
                Û₂ += c2[j] * uprev
                t̂₂ += c2[j] * t
            elseif i == 1
                Û₁ += c1[j] * uᵢ₋₁
                t̂₁ += c1[j] * tᵢ₋₁
                Û₂ += c2[j] * uᵢ₋₁
                t̂₂ += c2[j] * tᵢ₋₁
            else
                Û₁ += c1[j] * uᵢ
                t̂₁ += c1[j] * tᵢ
                Û₂ += c2[j] * uᵢ
                t̂₂ += c2[j] * tᵢ
            end
        end
        if i > 1 && i < mdeg + 1
            uᵢ₋₂ = uᵢ₋₁
            uᵢ₋₁ = uᵢ
            tᵢ₋₂ = tᵢ₋₁
            tᵢ₋₁ = tᵢ
        end
    end

    if (W.dW isa Number) || (length(W.dW) == 1)
        Gₛ = integrator.f.g(Û₁, p, t̂₁)
        uₓ += Gₛ * W.dW

        uₓ = integrator.f(uₓ, p, tₓ)
        u += (1 // 2 * dt) * uₓ + Gₛ * ((W.dW^2 - abs(dt)) / (η₁ * sqrt_dt) - W.dW)
        Û₁ -= (η₁ * sqrt_dt / 2) * Gₛ
        Û₂ += (η₁ * sqrt_dt / 2) * Gₛ

        Gₛ = integrator.f.g(Û₂, p, t̂₂)
        u += Gₛ * W.dW

        Gₛ = integrator.f.g(Û₁, p, t̂₁)
        u += Gₛ * (W.dW - (W.dW^2 - abs(dt)) / (η₁ * sqrt_dt))
    elseif is_diagonal_noise(integrator.sol.prob)
        Gₛ = integrator.f.g(Û₁, p, t̂₁)
        uᵢ₋₁ = Gₛ .* W.dW

        Û₁ -= (1 // 2 * η₁ * sqrt_dt) * Gₛ
        Û₂ += (1 // 2 * η₁ * sqrt_dt) * Gₛ

        uₓ += uᵢ₋₁
        uₓ = integrator.f(uₓ, p, tₓ)
        u += (1 // 2) * dt * uₓ

        u .+= Gₛ .* ((W.dW .^ 2 .- abs(dt)) ./ (η₁ * sqrt_dt) .- W.dW)

        Gₛ = integrator.f.g(Û₂, p, t̂₂)
        u .+= Gₛ .* W.dW

        Gₛ = integrator.f.g(Û₁, p, t̂₁)
        u .-= Gₛ .* ((W.dW .^ 2 .- abs(dt)) ./ (η₁ * sqrt_dt) .- W.dW)
    else
        Gₛ = integrator.f.g(Û₁, p, t̂₁)

        for i in 1:length(W.dW)
            (i == 1) && (uᵢ₋₁ = @view(Gₛ[:, i]) * W.dW[i])
            (i != 1) && (uᵢ₋₁ += @view(Gₛ[:, i]) * W.dW[i])
        end

        uₓ += uᵢ₋₁
        uₓ = integrator.f(uₓ, p, tₓ)

        u += (1 // 2 * dt) * uₓ - uᵢ₋₁

        for i in 1:length(W.dW)
            uᵢ₋₁ = Û₁ - (1 // 2 * η₁ * sqrt_dt) * @view(Gₛ[:, i])
            Gₛ₁ = integrator.f.g(uᵢ₋₁, p, t̂₁)
            u += @view(Gₛ₁[:, i]) * W.dW[i] +
                (@view(Gₛ[:, i]) - @view(Gₛ₁[:, i])) * ((W.dW[i]^2 - abs(dt)) / (η₁ * sqrt_dt))
        end

        for i in 1:length(W.dW)
            for j in 1:length(W.dW)
                WikJ = i > j ? (1 // 2) * (1 + η₂) * W.dW[j] :
                    i < j ? (1 // 2) * (1 - η₂) * W.dW[j] : (1 // 2) * (η₁ * sqrt_dt)

                uᵢ₋₁ += @view(Gₛ[:, j]) * WikJ
            end
            Gₛ₁ = integrator.f.g(uᵢ₋₁, p, t̂₂)
            u += @view(Gₛ₁[:, i]) * W.dW[i]
        end
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::TangXiaoSROCK2Cache)
    (; uᵢ, uₓ, uᵢ₋₁, uᵢ₋₂, Û₁, Û₂, k, Gₛ, Gₛ₁) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    (; recf, recf2, mα, mσ, mτ, mn̂, c1, c2) = cache.constantcache

    n̂ = mn̂[integrator.alg.version_num]
    ccache = cache.constantcache

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    (integrator.alg.version_num <= 2) &&
        (
        ccache.mdeg = Int(
            floor(
                sqrt(
                    (
                        abs(dt) * integrator.opts.internalnorm(
                            integrator.eigen_est, t
                        ) + 1.5
                    ) / 0.811
                ) + 1
            )
        )
    )
    (integrator.alg.version_num > 2) &&
        (
        ccache.mdeg = Int(
            floor(
                sqrt(
                    (
                        abs(dt) * integrator.opts.internalnorm(
                            integrator.eigen_est, t
                        ) + 1.5
                    ) / 0.611
                ) + 1
            )
        )
    )
    ccache.mdeg = max(4, min(ccache.mdeg, 200)) - 2
    choose_deg!(integrator, cache)

    mdeg = ccache.mdeg
    start = ccache.start
    deg_index = ccache.deg_index
    α = oftype(t, 1.33)
    σ = (1 - α) * 1 // 2 + α * mσ[deg_index]
    τ = 1 // 2 * ((1 - α)^2) + 2 * α * (1 - α) * mσ[deg_index] +
        (α^2) * (mσ[deg_index] * (mσ[deg_index] + mτ[deg_index]))

    η₁ = (rand() < 1 // 2) ? -1 : 1
    η₂ = (rand() < 1 // 2) ? -1 : 1
    sqrt_dt = sqrt(abs(dt))

    @.. Û₁ = zero(eltype(u))
    @.. Û₂ = zero(eltype(u))
    t̂₁ = t̂₂ = tₓ = zero(t)
    tᵢ = tᵢ₋₁ = tᵢ₋₂ = t

    for i in 0:(mdeg + 1)
        if i == 1
            μ = recf[start]
            tᵢ = tᵢ₋₁ = t + α * dt * μ

            @.. uᵢ₋₂ = uprev
            integrator.f(k, uprev, p, t)
            @.. uᵢ₋₁ = uprev + α * dt * μ * k
        elseif i > 1 && i <= mdeg
            μ, ν,
                κ = recf[start + 2 * (i - 2) + 1], 1 + recf[start + 2 * (i - 2) + 2],
                recf[start + 2 * (i - 2) + 2]

            integrator.f(k, uᵢ₋₁, p, tᵢ₋₁)
            @.. uᵢ = α * dt * μ * k + ν * uᵢ₋₁ - κ * uᵢ₋₂
            tᵢ = α * dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂
        elseif i == mdeg + 1
            μ, ν,
                κ = recf2[(deg_index - 1) * 4 + 1], 1 + recf2[(deg_index - 1) * 4 + 2],
                recf2[(deg_index - 1) * 4 + 2]

            integrator.f(k, uᵢ₋₁, p, tᵢ₋₁)

            tₓ = tᵢ₋₁ + 2 * τ * dt
            @.. uₓ = uᵢ₋₁ + (2 * τ * dt) * k
            @.. u = uᵢ₋₁ + (2 * σ - 1 // 2) * dt * k

            @.. uᵢ = α * dt * μ * k + ν * uᵢ₋₁ - κ * uᵢ₋₂
            tᵢ = α * dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂
        end

        j = i - mdeg - 1 + n̂
        if j > 0
            j += ccache.start_mcs - 1
            if i == 0
                @.. Û₁ += c1[j] * uprev
                t̂₁ += c1[j] * t
                @.. Û₂ += c2[j] * uprev
                t̂₂ += c2[j] * t
            elseif i == 1
                @.. Û₁ += c1[j] * uᵢ₋₁
                t̂₁ += c1[j] * tᵢ₋₁
                @.. Û₂ += c2[j] * uᵢ₋₁
                t̂₂ += c2[j] * tᵢ₋₁
            else
                @.. Û₁ += c1[j] * uᵢ
                t̂₁ += c1[j] * tᵢ
                @.. Û₂ += c2[j] * uᵢ
                t̂₂ += c2[j] * tᵢ
            end
        end

        if i > 1 && i < mdeg + 1
            @.. uᵢ₋₂ = uᵢ₋₁
            @.. uᵢ₋₁ = uᵢ
            tᵢ₋₂ = tᵢ₋₁
            tᵢ₋₁ = tᵢ
        end
    end

    if (W.dW isa Number) || (length(W.dW) == 1) || is_diagonal_noise(integrator.sol.prob)
        integrator.f.g(Gₛ, Û₁, p, t̂₁)
        @.. uₓ += Gₛ * W.dW

        integrator.f(k, uₓ, p, tₓ)
        @.. u += (1 // 2 * dt) * k + Gₛ * ((W.dW^2 - abs(dt)) / (η₁ * sqrt_dt) - W.dW)
        @.. Û₁ -= (η₁ * sqrt_dt / 2) * Gₛ
        @.. Û₂ += (η₁ * sqrt_dt / 2) * Gₛ

        integrator.f.g(Gₛ, Û₂, p, t̂₂)
        @.. u += Gₛ * W.dW

        integrator.f.g(Gₛ, Û₁, p, t̂₁)
        @.. u += Gₛ * (W.dW - (W.dW^2 - abs(dt)) / (η₁ * sqrt_dt))
    else
        integrator.f.g(Gₛ, Û₁, p, t̂₁)

        for i in 1:length(W.dW)
            (i == 1) && (@.. uᵢ₋₁ = @view(Gₛ[:, i]) * W.dW[i])
            (i > 1) && (@.. uᵢ₋₁ += @view(Gₛ[:, i]) * W.dW[i])
        end

        @.. uₓ += uᵢ₋₁
        integrator.f(k, uₓ, p, tₓ)

        @.. u += (1 // 2 * dt) * k - uᵢ₋₁

        for i in 1:length(W.dW)
            @.. uᵢ₋₁ = Û₁ - (1 // 2 * η₁ * sqrt_dt) * @view(Gₛ[:, i])
            integrator.f.g(Gₛ₁, uᵢ₋₁, p, t̂₁)
            @.. u += @view(Gₛ₁[:, i]) * W.dW[i] +
                (
                @view(Gₛ[:, i]) -
                    @view(Gₛ₁[:, i])
            ) * ((W.dW[i]^2 - abs(dt)) / (η₁ * sqrt_dt))
        end

        for i in 1:length(W.dW)
            for j in 1:length(W.dW)
                WikJ = i > j ? (1 // 2) * (1 + η₂) * W.dW[j] :
                    i < j ? (1 // 2) * (1 - η₂) * W.dW[j] : (1 // 2) * (η₁ * sqrt_dt)

                @.. uᵢ₋₁ += @view(Gₛ[:, j]) * WikJ
            end
            @.. uᵢ₋₁ += Û₂
            integrator.f.g(Gₛ₁, uᵢ₋₁, p, t̂₂)
            @.. u += @view(Gₛ₁[:, i]) * W.dW[i]
        end
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::KomBurSROCK2ConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    (; recf, mσ, mτ, mδ) = cache

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    cache.mdeg = Int(floor(sqrt((2 * abs(dt) * integrator.opts.internalnorm(integrator.eigen_est, t) + 1.5) / 0.811) + 1))
    cache.mdeg = max(6, min(cache.mdeg, 200)) - 2
    choose_deg!(integrator, cache)

    # here mdeg == s in the paper
    mdeg = cache.mdeg + 2
    start = cache.start
    deg_index = cache.deg_index
    σ = mσ[deg_index]
    τ = mτ[deg_index]

    sqrt_dt = sqrt(abs(dt))

    tᵢ₋₂ = t
    uᵢ₋₂ = uprev
    μ = recf[start]
    tᵢ = tᵢ₋₁ = t + dt * μ
    u = integrator.f(uprev, p, t)
    uᵢ₋₁ = uprev + dt * μ * u

    for i in 2:(mdeg - 4)
        μ, θ = recf[start + 2 * (i - 2) + 1], recf[start + 2 * (i - 2) + 2]
        u = integrator.f(uᵢ₋₁, p, tᵢ₋₁)
        u = dt * μ * u + (1 + θ) * uᵢ₋₁ - θ * uᵢ₋₂
        tᵢ = dt * μ + (1 + θ) * tᵢ₋₁ - θ * tᵢ₋₂

        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = u
        tᵢ₋₂ = tᵢ₋₁
        tᵢ₋₁ = tᵢ
    end

    μₛ₋₃ = dt * recf[start + 2 * (mdeg - 5) + 1]
    θₛ₋₃ = recf[start + 2 * (mdeg - 5) + 2]
    μₛ₋₂ = dt * recf[start + 2 * (mdeg - 4) + 1]
    θₛ₋₂ = recf[start + 2 * (mdeg - 4) + 2]

    δ₁ = dt * mδ[(deg_index - 1) * 8 + 1]
    δ₂ = dt * mδ[(deg_index - 1) * 8 + 2]
    δ₃ = dt * mδ[(deg_index - 1) * 8 + 3]
    δ₄ = mδ[(deg_index - 1) * 8 + 4]
    δ₅ = mδ[(deg_index - 1) * 8 + 5]
    δ₆ = mδ[(deg_index - 1) * 8 + 6]
    δ₇ = mδ[(deg_index - 1) * 8 + 7]
    δ₈ = mδ[(deg_index - 1) * 8 + 8]

    C₁ = (1 + θₛ₋₂) * μₛ₋₃

    u = uᵢ₋₁ + θₛ₋₃ * (1 + θₛ₋₂) * (uᵢ₋₁ - uᵢ₋₂)
    uᵢ₋₁ += θₛ₋₃ * (uᵢ₋₁ - uᵢ₋₂)
    uᵢ₋₂ = u
    ttmp = tᵢ₋₁ + θₛ₋₃ * (1 + θₛ₋₂) * (tᵢ₋₁ - tᵢ₋₂)
    tᵢ₋₁ += θₛ₋₃ * (tᵢ₋₁ - tᵢ₋₂)
    tᵢ₋₂ = ttmp

    if W.dW isa Number || length(W.dW) == 1 || is_diagonal_noise(integrator.sol.prob)
        # stage s-3
        yₛ₋₃ = integrator.f(uᵢ₋₁, p, tᵢ₋₁)
        utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃
        ttmp = tᵢ₋₁ + μₛ₋₃
        Xₛ₋₃ = integrator.f.g(utmp, p, ttmp)
        u += C₁ * yₛ₋₃ + 1 // 8 .* W.dW .* Xₛ₋₃

        #stage s-2
        utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + δ₄ .* W.dW .* Xₛ₋₃
        ttmp = tᵢ₋₁ + μₛ₋₃
        yₛ₋₂ = integrator.f(utmp, p, ttmp)

        utmp = uᵢ₋₂ + C₁ * yₛ₋₃ + 2 // 3 .* W.dW .* Xₛ₋₃
        ttmp = tᵢ₋₂ + C₁
        Xₛ₋₂ = integrator.f.g(utmp, p, ttmp)
        u += μₛ₋₂ * yₛ₋₂ + 3 // 8 .* W.dW .* Xₛ₋₂

        #stage s-1
        utmp = uᵢ₋₂ + μₛ₋₂ * yₛ₋₂ + C₁ * yₛ₋₃ + δ₅ .* W.dW .* Xₛ₋₃ + δ₄ .* W.dW .* Xₛ₋₂
        ttmp = tᵢ₋₂ + μₛ₋₂ + C₁
        yₛ₋₁ = integrator.f(utmp, p, ttmp)

        utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + δ₁ * yₛ₋₂ + 1 // 12 .* W.dW .* Xₛ₋₃ + 1 // 4 .* W.dW .* Xₛ₋₂
        ttmp = tᵢ₋₁ + μₛ₋₃ + δ₁
        Xₛ₋₁ = integrator.f.g(utmp, p, ttmp)
        u += (σ - τ) * dt * yₛ₋₁ + 3 // 8 .* W.dW .* Xₛ₋₁

        #stage s
        utmp = uᵢ₋₂ + C₁ * yₛ₋₃ + μₛ₋₂ * yₛ₋₂ + σ * dt * yₛ₋₁ + δ₆ .* W.dW .* Xₛ₋₃ +
            δ₇ .* W.dW .* Xₛ₋₂ + δ₈ .* W.dW .* Xₛ₋₁
        ttmp = tᵢ₋₂ + C₁ + μₛ₋₂ + σ * dt
        utmp = integrator.f(utmp, p, ttmp)
        u += (σ + τ) * dt * utmp

        utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + δ₂ * yₛ₋₂ + δ₃ * yₛ₋₁ - 5 // 4 .* W.dW .* Xₛ₋₃ +
            1 // 4 .* W.dW .* Xₛ₋₂ + 2 .* W.dW .* Xₛ₋₁
        ttmp = tᵢ₋₁ + μₛ₋₃ + δ₂ + δ₃
        Xₛ₋₁ = integrator.f.g(utmp, p, ttmp)
        u += 1 // 8 .* W.dW .* Xₛ₋₁
    else
        vec_χ = 2 .* floor.(1 // 2 .+ false .* W.dW .+ rand(length(W.dW))) .- 1

        # stage s-3
        yₛ₋₃ = integrator.f(uᵢ₋₁, p, tᵢ₋₁)
        utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃
        ttmp = tᵢ₋₁ + μₛ₋₃
        Xₛ₋₃ = integrator.f.g(utmp, p, ttmp)
        SXₛ₋₃ = Xₛ₋₃ * W.dW

        u += C₁ * yₛ₋₃ + 1 // 8 * SXₛ₋₃

        #stage s-2
        utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + δ₄ * SXₛ₋₃
        ttmp = tᵢ₋₁ + μₛ₋₃
        yₛ₋₂ = integrator.f(utmp, p, ttmp)
        Xₛ₋₂ = zero(Xₛ₋₃)

        for i in 1:length(W.dW)
            WikRange = W.dW .* (1:length(W.dW) .== i)
            # utmp = uᵢ₋₂ + C₁*yₛ₋₃ + 2//3*@view(Xₛ₋₃[:,i])*W.dW[i]
            utmp = uᵢ₋₂ + C₁ * yₛ₋₃ + 2 // 3 * (Xₛ₋₃ * WikRange)
            ttmp = tᵢ₋₂ + C₁
            Gₛ = integrator.f.g(utmp, p, ttmp)
            WikRange = 1 .* (1:length(W.dW) .== i)
            # @view(Xₛ₋₂[:,i]) .=  @view(Gₛ[:,i])
            Xₛ₋₂ .+= Gₛ .* WikRange
        end
        SXₛ₋₂ = Xₛ₋₂ * W.dW
        u += μₛ₋₂ * yₛ₋₂ + 3 // 8 * SXₛ₋₂

        #stage s-1
        utmp = uᵢ₋₂ + μₛ₋₂ * yₛ₋₂ + C₁ * yₛ₋₃ + δ₅ * SXₛ₋₃ + δ₄ * SXₛ₋₂
        ttmp = tᵢ₋₂ + μₛ₋₂ + C₁
        yₛ₋₁ = integrator.f(utmp, p, ttmp)

        Xₛ₋₁ = zero(Xₛ₋₂)
        for i in 1:length(W.dW)
            WikRange = W.dW .* (1:length(W.dW) .== i)
            # utmp = uᵢ₋₁ + μₛ₋₃*yₛ₋₃ + δ₁*yₛ₋₂ - 1//6*W.dW[i]*@view(Xₛ₋₃[:,i]) - 1//2*W.dW[i]*@view(Xₛ₋₂[:,i]) + 1//4*SXₛ₋₃ + 3//4*SXₛ₋₂
            utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + δ₁ * yₛ₋₂ + 1 // 4 * SXₛ₋₃ + 3 // 4 * SXₛ₋₂ -
                1 // 6 * (Xₛ₋₃ * WikRange) - 1 // 2 * (Xₛ₋₂ * WikRange)
            ttmp = tᵢ₋₁ + μₛ₋₃ + δ₁
            Gₛ = integrator.f.g(utmp, p, ttmp)
            WikRange = 1 .* (1:length(W.dW) .== i)
            #@view(Xₛ₋₁[:,i]) .= @view(Gₛ[:,i])
            Xₛ₋₁ .+= Gₛ .* WikRange
        end
        SXₛ₋₁ = Xₛ₋₁ * W.dW
        u += (σ - τ) * dt * yₛ₋₁ + 3 // 8 * SXₛ₋₁

        #stage s
        utmp = uᵢ₋₂ + C₁ * yₛ₋₃ + μₛ₋₂ * yₛ₋₂ + σ * dt * yₛ₋₁ + δ₆ * SXₛ₋₃ + δ₇ * SXₛ₋₂ + δ₈ * SXₛ₋₁
        ttmp = tᵢ₋₂ + C₁ + μₛ₋₂ + σ * dt
        utmp = integrator.f(utmp, p, ttmp)
        u += (σ + τ) * dt * utmp

        SXₛ₋₁ = zero(uprev)
        for i in 1:length(W.dW)
            WikRange = W.dW .* (1:length(W.dW) .== i)
            # utmp = uᵢ₋₁ + μₛ₋₃*yₛ₋₃ + δ₂*yₛ₋₂ + δ₃*yₛ₋₁ - 3//2*W.dW[i]*@view(Xₛ₋₃[:,i]) - 1//2*W.dW[i]*@view(Xₛ₋₂[:,i]) + 2*W.dW[i]*@view(Xₛ₋₁[:,i]) + 1//4*SXₛ₋₃ + 3//4*SXₛ₋₂
            utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + δ₂ * yₛ₋₂ + δ₃ * yₛ₋₁ + 1 // 4 * SXₛ₋₃ + 3 // 4 * SXₛ₋₂ -
                3 // 2 * (Xₛ₋₃ * WikRange) - 1 // 2 * (Xₛ₋₂ * WikRange) + 2 * (Xₛ₋₁ * WikRange)
            ttmp = tᵢ₋₁ + μₛ₋₃ + δ₂ + δ₃
            Gₛ = integrator.f.g(utmp, p, ttmp)
            # SXₛ₋₁ += W.dW[i]*@view(Gₛ[:,i])
            SXₛ₋₁ += Gₛ * WikRange
        end
        u += 1 // 8 * SXₛ₋₁

        for i in 1:length(W.dW)
            SXₛ₋₁ = zero(uprev)
            for j in 1:length(W.dW)
                WikRange = 1 .* (1:length(W.dW) .== j)
                if j != i
                    # SXₛ₋₁ += (i > j ? -vec_χ[i]*W.dW[j] : vec_χ[j]*W.dW[i] )*@view(Xₛ₋₃[:,j])
                    SXₛ₋₁ += (i > j ? -vec_χ[i] * W.dW[j] : vec_χ[j] * W.dW[i]) * (Xₛ₋₃ * WikRange)
                end
            end
            ttmp = tᵢ₋₁ + μₛ₋₃
            utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ - 1 // 4 * SXₛ₋₁
            Gₛ = integrator.f.g(utmp, p, ttmp)
            utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + 1 // 4 * SXₛ₋₁
            Xₛ₋₁ = integrator.f.g(utmp, p, ttmp)
            sqrt_dt *= (length(W.dW) - 1)
            WikRange = 1 .* (1:length(W.dW) .== i)
            # u += sqrt_dt*(@view(Xₛ₋₁[:,i]) - @view(Gₛ[:,i]))
            u += sqrt_dt * (Xₛ₋₁ * WikRange - Gₛ * WikRange)
        end
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::KomBurSROCK2Cache)
    (;
        utmp, uᵢ₋₁, uᵢ₋₂, k, yₛ₋₁, yₛ₋₂, yₛ₋₃, SXₛ₋₁, SXₛ₋₂,
        SXₛ₋₃, Gₛ, Xₛ₋₁, Xₛ₋₂, Xₛ₋₃, vec_χ, WikRange,
    ) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    (; recf, mσ, mτ, mδ) = cache.constantcache

    ccache = cache.constantcache
    gen_prob = !(
        (is_diagonal_noise(integrator.sol.prob)) || (W.dW isa Number) ||
            (length(W.dW) == 1)
    )

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    ccache.mdeg = Int(
        floor(
            sqrt(
                (
                    2 * abs(dt) * integrator.opts.internalnorm(
                        integrator.eigen_est, t
                    ) + 1.5
                ) / 0.811
            ) + 1
        )
    )
    ccache.mdeg = max(6, min(ccache.mdeg, 200)) - 2
    choose_deg!(integrator, cache)

    # here mdeg == s in the paper
    mdeg = ccache.mdeg + 2
    start = ccache.start
    deg_index = ccache.deg_index
    σ = mσ[deg_index]
    τ = mτ[deg_index]

    sqrt_dt = sqrt(abs(dt))
    (gen_prob) && (vec_χ .= 2 .* floor.(1 // 2 .+ false .* vec_χ .+ rand(length(vec_χ))) .- 1)

    tᵢ₋₂ = t
    @.. uᵢ₋₂ = uprev
    μ = recf[start]
    tᵢ = tᵢ₋₁ = t + dt * μ
    integrator.f(k, uprev, p, t)
    @.. uᵢ₋₁ = uprev + dt * μ * k

    for i in 2:(mdeg - 4)
        μ, θ = recf[start + 2 * (i - 2) + 1], recf[start + 2 * (i - 2) + 2]
        integrator.f(k, uᵢ₋₁, p, tᵢ₋₁)
        @.. u = dt * μ * k + (1 + θ) * uᵢ₋₁ - θ * uᵢ₋₂
        tᵢ = dt * μ + (1 + θ) * tᵢ₋₁ - θ * tᵢ₋₂

        @.. uᵢ₋₂ = uᵢ₋₁
        @.. uᵢ₋₁ = u
        tᵢ₋₂ = tᵢ₋₁
        tᵢ₋₁ = tᵢ
    end
    μₛ₋₃ = dt * recf[start + 2 * (mdeg - 5) + 1]
    θₛ₋₃ = recf[start + 2 * (mdeg - 5) + 2]
    μₛ₋₂ = dt * recf[start + 2 * (mdeg - 4) + 1]
    θₛ₋₂ = recf[start + 2 * (mdeg - 4) + 2]

    δ₁ = dt * mδ[(deg_index - 1) * 8 + 1]
    δ₂ = dt * mδ[(deg_index - 1) * 8 + 2]
    δ₃ = dt * mδ[(deg_index - 1) * 8 + 3]
    δ₄ = mδ[(deg_index - 1) * 8 + 4]
    δ₅ = mδ[(deg_index - 1) * 8 + 5]
    δ₆ = mδ[(deg_index - 1) * 8 + 6]
    δ₇ = mδ[(deg_index - 1) * 8 + 7]
    δ₈ = mδ[(deg_index - 1) * 8 + 8]

    C₁ = (1 + θₛ₋₂) * μₛ₋₃

    @.. u = uᵢ₋₁ + θₛ₋₃ * (1 + θₛ₋₂) * (uᵢ₋₁ - uᵢ₋₂)
    @.. uᵢ₋₁ += θₛ₋₃ * (uᵢ₋₁ - uᵢ₋₂)
    @.. uᵢ₋₂ = u
    ttmp = tᵢ₋₁ + θₛ₋₃ * (1 + θₛ₋₂) * (tᵢ₋₁ - tᵢ₋₂)
    tᵢ₋₁ += θₛ₋₃ * (tᵢ₋₁ - tᵢ₋₂)
    tᵢ₋₂ = ttmp

    if W.dW isa Number || length(W.dW) == 1 || is_diagonal_noise(integrator.sol.prob)
        # stage s-3
        integrator.f(yₛ₋₃, uᵢ₋₁, p, tᵢ₋₁)
        @.. utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃
        ttmp = tᵢ₋₁ + μₛ₋₃
        integrator.f.g(Xₛ₋₃, utmp, p, ttmp)
        @.. u += C₁ * yₛ₋₃ + 1 // 8 * W.dW * Xₛ₋₃

        #stage s-2
        @.. utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + δ₄ * W.dW * Xₛ₋₃
        ttmp = tᵢ₋₁ + μₛ₋₃
        integrator.f(yₛ₋₂, utmp, p, ttmp)

        @.. utmp = uᵢ₋₂ + C₁ * yₛ₋₃ + 2 // 3 * W.dW * Xₛ₋₃
        ttmp = tᵢ₋₂ + C₁
        integrator.f.g(Xₛ₋₂, utmp, p, ttmp)
        @.. u += μₛ₋₂ * yₛ₋₂ + 3 // 8 * W.dW * Xₛ₋₂

        #stage s-1
        @.. utmp = uᵢ₋₂ + μₛ₋₂ * yₛ₋₂ + C₁ * yₛ₋₃ + δ₅ * W.dW * Xₛ₋₃ + δ₄ * W.dW * Xₛ₋₂
        ttmp = tᵢ₋₂ + μₛ₋₂ + C₁
        integrator.f(yₛ₋₁, utmp, p, ttmp)

        @.. utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + δ₁ * yₛ₋₂ + 1 // 12 * W.dW * Xₛ₋₃ + 1 // 4 * W.dW * Xₛ₋₂
        ttmp = tᵢ₋₁ + μₛ₋₃ + δ₁
        integrator.f.g(Xₛ₋₁, utmp, p, ttmp)
        @.. u += (σ - τ) * dt * yₛ₋₁ + 3 // 8 * W.dW * Xₛ₋₁

        #stage s
        @.. utmp = uᵢ₋₂ + C₁ * yₛ₋₃ + μₛ₋₂ * yₛ₋₂ + σ * dt * yₛ₋₁ + δ₆ * W.dW * Xₛ₋₃ + δ₇ * W.dW * Xₛ₋₂ +
            δ₈ * W.dW * Xₛ₋₁
        ttmp = tᵢ₋₂ + C₁ + μₛ₋₂ + σ * dt

        @.. uᵢ₋₁ = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + δ₂ * yₛ₋₂ + δ₃ * yₛ₋₁ - 5 // 4 * W.dW * Xₛ₋₃ + 1 // 4 * W.dW * Xₛ₋₂ +
            2 * W.dW * Xₛ₋₁
        tᵢ₋₁ = tᵢ₋₁ + μₛ₋₃ + δ₂ + δ₃

        integrator.f(yₛ₋₁, utmp, p, ttmp)
        integrator.f.g(Xₛ₋₁, uᵢ₋₁, p, tᵢ₋₁)
        @.. u += (σ + τ) * dt * yₛ₋₁ + 1 // 8 * W.dW * Xₛ₋₁
    else
        # stage s-3
        integrator.f(yₛ₋₃, uᵢ₋₁, p, tᵢ₋₁)
        @.. utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃
        ttmp = tᵢ₋₁ + μₛ₋₃
        integrator.f.g(Xₛ₋₃, utmp, p, ttmp)
        mul!(SXₛ₋₃, Xₛ₋₃, W.dW)
        @.. u += C₁ * yₛ₋₃ + 1 // 8 * SXₛ₋₃

        #stage s-2
        @.. utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + δ₄ * SXₛ₋₃
        ttmp = tᵢ₋₁ + μₛ₋₃
        integrator.f(yₛ₋₂, utmp, p, ttmp)
        for i in 1:length(W.dW)
            WikRange .= 1 .* (1:length(W.dW) .== i) .* W.dW
            # @.. utmp = uᵢ₋₂ + C₁*yₛ₋₃ + 2//3*@view(Xₛ₋₃[:,i])*W.dW[i]
            mul!(SXₛ₋₂, Xₛ₋₃, WikRange)
            @.. utmp = uᵢ₋₂ + C₁ * yₛ₋₃ + 2 // 3 * SXₛ₋₂
            ttmp = tᵢ₋₂ + C₁
            integrator.f.g(Gₛ, utmp, p, ttmp)
            WikRange .= 1 .* (1:length(W.dW) .== i)
            # @.. @view(Xₛ₋₂[:,i]) =  @view(Gₛ[:,i])
            @.. Xₛ₋₂ = Gₛ * W.dW
        end
        mul!(SXₛ₋₂, Xₛ₋₂, W.dW)
        @.. u += μₛ₋₂ * yₛ₋₂ + 3 // 8 * SXₛ₋₂

        #stage s-1
        @.. utmp = uᵢ₋₂ + μₛ₋₂ * yₛ₋₂ + C₁ * yₛ₋₃ + δ₅ * SXₛ₋₃ + δ₄ * SXₛ₋₂
        ttmp = tᵢ₋₂ + μₛ₋₂ + C₁
        integrator.f(yₛ₋₁, utmp, p, ttmp)
        for i in 1:length(W.dW)
            WikRange .= 1 .* (1:length(W.dW) .== i) .* W.dW
            # @.. utmp = uᵢ₋₁ + μₛ₋₃*yₛ₋₃ + δ₁*yₛ₋₂ - 1//6*W.dW[i]*@view(Xₛ₋₃[:,i]) - 1//2*W.dW[i]*@view(Xₛ₋₂[:,i]) + 1//4*SXₛ₋₃ + 3//4*SXₛ₋₂
            @.. utmp = uᵢ₋₁ + μₛ₋₃ * yₛ₋₃ + δ₁ * yₛ₋₂ + 1 // 4 * SXₛ₋₃ + 3 // 4 * SXₛ₋₂
            mul!(SXₛ₋₁, Xₛ₋₃, WikRange)
            @.. utmp += 1 // 6 * SXₛ₋₁
            mul!(SXₛ₋₁, Xₛ₋₂, WikRange)
            @.. utmp += 1 // 2 * SXₛ₋₁
            ttmp = tᵢ₋₁ + μₛ₋₃ + δ₁
            integrator.f.g(Gₛ, utmp, p, ttmp)
            WikRange .= 1 .* (1:length(W.dW) .== i)
            # @.. @view(Xₛ₋₁[:,i]) = @view(Gₛ[:,i])
            @.. Xₛ₋₁ = Gₛ * WikRange
        end
        mul!(SXₛ₋₁, Xₛ₋₁, W.dW)
        @.. u += (σ - τ) * dt * yₛ₋₁ + 3 // 8 * SXₛ₋₁

        #stage s
        @.. utmp = uᵢ₋₂ + C₁ * yₛ₋₃ + μₛ₋₂ * yₛ₋₂ + σ * dt * yₛ₋₁ + δ₆ * SXₛ₋₃ + δ₇ * SXₛ₋₂ + δ₈ * SXₛ₋₁
        ttmp = tᵢ₋₂ + C₁ + μₛ₋₂ + σ * dt
        # utmp = integrator.f(utmp,p,ttmp)
        # u += (τ/σ)*dt*utmp
        #
        # SXₛ₋₁ = zero(uprev)
        # for i in 1:length(W.dW)
        #   utmp = uᵢ₋₁ + θₛ₋₃*(uᵢ₋₁ - uᵢ₋₂) + μₛ₋₃*yₛ₋₃ + δ₂*yₛ₋₂ + δ₃*yₛ₋₁ - 3//2*W.dW[i]*@view(Xₛ₋₃[:,i]) -
        #           1//2*W.dW[i]*@view(Xₛ₋₂[:,i]) + 2*W.dW[i]*@view(Xₛ₋₁[:,i]) + 1//4*SXₛ₋₃ + 3//4*SXₛ₋₂
        #   ttmp = tᵢ₋₁ + θₛ₋₃*(tᵢ₋₁ - tᵢ₋₂) + μₛ₋₃ + δ₂ + δ₃
        #   Gₛ = integrator.f.g(utmp,p,ttmp)
        #   SXₛ₋₁ += W.dW[i]*@view(Gₛ[:,i])
        # end
        # u += 1//8*SXₛ₋₁

        # memory optimisation
        @.. uᵢ₋₁ += μₛ₋₃ * yₛ₋₃
        tᵢ₋₁ += μₛ₋₃
        #now we have uᵢ₋₂ and yₛ₋₃ free
        integrator.f(yₛ₋₃, utmp, p, ttmp)
        @.. SXₛ₋₁ = zero(uprev)
        for i in 1:length(W.dW)
            WikRange .= 1 .* (1:length(W.dW) .== i) .* W.dW
            # @.. uᵢ₋₂ = uᵢ₋₁ + δ₂*yₛ₋₂ + δ₃*yₛ₋₁ - 3//2*W.dW[i]*@view(Xₛ₋₃[:,i]) - 1//2*W.dW[i]*@view(Xₛ₋₂[:,i]) + 2*W.dW[i]*@view(Xₛ₋₁[:,i]) + 1//4*SXₛ₋₃ + 3//4*SXₛ₋₂
            @.. uᵢ₋₂ = uᵢ₋₁ + δ₂ * yₛ₋₂ + δ₃ * yₛ₋₁ + 1 // 4 * SXₛ₋₃ + 3 // 4 * SXₛ₋₂
            mul!(utmp, Xₛ₋₃, WikRange)
            @.. uᵢ₋₂ -= 3 // 2 * utmp
            mul!(utmp, Xₛ₋₂, WikRange)
            @.. uᵢ₋₂ -= 1 // 2 * utmp
            mul!(utmp, Xₛ₋₁, WikRange)
            @.. uᵢ₋₂ += 2 * utmp

            tᵢ₋₂ = tᵢ₋₁ + δ₂ + δ₃
            integrator.f.g(Gₛ, uᵢ₋₂, p, tᵢ₋₂)
            # @.. SXₛ₋₁ += W.dW[i]*@view(Gₛ[:,i])
            mul!(utmp, Gₛ, WikRange)
            @.. SXₛ₋₁ += utmp
        end
        @.. u += (σ + τ) * dt * yₛ₋₃ + 1 // 8 * SXₛ₋₁

        for i in 1:length(W.dW)
            @.. SXₛ₋₁ = zero(uprev)
            for j in 1:length(W.dW)
                WikRange .= 1 .* (1:length(W.dW) .== j)
                mul!(utmp, Xₛ₋₃, WikRange)
                if j != i
                    # (i > j) && (@.. SXₛ₋₁ += -vec_χ[i]*W.dW[j]*@view(Xₛ₋₃[:,j]))
                    # (i < j) && (@.. SXₛ₋₁ += vec_χ[j]*W.dW[i]*@view(Xₛ₋₃[:,j]))
                    (i > j) && (@.. SXₛ₋₁ += -vec_χ[i] * W.dW[j] * utmp)
                    (i < j) && (@.. SXₛ₋₁ += vec_χ[j] * W.dW[i] * utmp)
                end
            end
            ttmp = tᵢ₋₁
            @.. utmp = uᵢ₋₁ - 1 // 4 * SXₛ₋₁
            integrator.f.g(Gₛ, utmp, p, ttmp)
            @.. utmp = uᵢ₋₁ + 1 // 4 * SXₛ₋₁
            integrator.f.g(Xₛ₋₁, utmp, p, ttmp)
            sqrt_dt *= (length(W.dW) - 1)
            # @.. u += sqrt_dt*(@view(Xₛ₋₁[:,i]) - @view(Gₛ[:,i]))
            WikRange .= 1 .* (1:length(W.dW) .== i)
            mul!(utmp, Xₛ₋₁, WikRange)
            @.. u += sqrt_dt * utmp
            mul!(utmp, Gₛ, WikRange)
            @.. u -= sqrt_dt * utmp
        end
    end

    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SROCKC2ConstantCache)
    (; t, dt, uprev, u, W, p, f) = integrator
    (; recf, mσ, mτ) = cache

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    cache.mdeg = Int(floor(sqrt((2 * abs(dt) * integrator.opts.internalnorm(integrator.eigen_est, t) + 1.5) / 0.811) + 1))
    cache.mdeg = max(3, min(cache.mdeg, 200)) - 2
    choose_deg!(integrator, cache)

    mdeg = cache.mdeg
    start = cache.start
    deg_index = cache.deg_index
    σ = mσ[deg_index]
    τ = mτ[deg_index]

    sqrt_dt = sqrt(dt)

    μ = recf[start]  # here κ = 0
    tᵢ = t + dt * μ
    tᵢ₋₁ = tᵢ
    tᵢ₋₂ = t

    # stage 1
    uᵢ₋₂ = uprev
    uᵢ = integrator.f(uprev, p, t)
    uᵢ₋₁ = uprev + dt * μ * uᵢ

    # stages 2 upto s-2
    for i in 2:mdeg
        μ, κ = recf[start + 2 * (i - 2) + 1], recf[start + 2 * (i - 2) + 2]
        ν = 1 + κ
        uᵢ = integrator.f(uᵢ₋₁, p, tᵢ₋₁)
        uᵢ = dt * μ * uᵢ + ν * uᵢ₋₁ - κ * uᵢ₋₂
        uᵢ₋₂ = uᵢ₋₁
        uᵢ₋₁ = uᵢ
        tᵢ = dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂
        tᵢ₋₂ = tᵢ₋₁
        tᵢ₋₁ = tᵢ
    end

    #2 stage-finishing procedure
    #stage s-1
    uᵢ = integrator.f(uᵢ₋₁, p, tᵢ₋₁)
    uᵢ₋₂ = uᵢ₋₁ + dt * σ * uᵢ
    u = uᵢ₋₁ + dt * (σ - τ) * uᵢ
    tᵢ₋₂ = tᵢ₋₁ + dt * σ

    #stage s
    uᵢ₋₂ = integrator.f(uᵢ₋₂, p, tᵢ₋₂)
    u += dt * (σ + τ) * uᵢ₋₂

    if (W.dW isa Number) || (length(W.dW) == 1) || is_diagonal_noise(integrator.sol.prob)
        Gₛ = integrator.f.g(uᵢ₋₁, p, tᵢ₋₁)
        u += Gₛ .* W.dW

        uᵢ₋₂ = uᵢ₋₁ .+ ((W.dW .^ 2 .- dt) ./ (2 .* sqrt_dt)) .* Gₛ
        Gₛ₁ = integrator.f.g(uᵢ₋₂, p, tᵢ₋₁)
        u += 1 // 2 .* sqrt_dt .* Gₛ₁

        uᵢ₋₂ = uᵢ₋₁ .- ((W.dW .^ 2 .- dt) ./ (2 .* sqrt_dt)) .* Gₛ
        Gₛ₁ = integrator.f.g(uᵢ₋₂, p, tᵢ₋₁)
        u -= 1 // 2 .* sqrt_dt .* Gₛ₁

    else
        Gₛ = integrator.f.g(uᵢ₋₁, p, tᵢ₋₁)
        uᵢ₋₂ = Gₛ * W.dW
        u += uᵢ₋₂

        for i in 1:length(W.dW)
            WikRange = 1 .* (1:length(W.dW) .== i)
            uᵢ = Gₛ * WikRange
            tmp = uᵢ₋₁ - 1 // 2 * sqrt_dt * uᵢ + 1 // 2 * (W.dW[i] / sqrt_dt) * uᵢ₋₂
            Gₛ₁ = integrator.f.g(tmp, p, tᵢ₋₁)
            u += 1 // 2 * sqrt_dt * (Gₛ₁ * WikRange)
            tmp = uᵢ₋₁ + 1 // 2 * sqrt_dt * uᵢ - 1 // 2 * (W.dW[i] / sqrt_dt) * uᵢ₋₂
            Gₛ₁ = integrator.f.g(tmp, p, tᵢ₋₁)
            u -= 1 // 2 * sqrt_dt * (Gₛ₁ * WikRange)
        end
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::SROCKC2Cache)
    (; uᵢ, tmp, uᵢ₋₁, uᵢ₋₂, k, Gₛ, Gₛ₁, WikRange) = cache
    (; t, dt, uprev, u, W, p, f) = integrator

    (; recf, mσ, mτ) = cache.constantcache
    ccache = cache.constantcache

    alg = unwrap_alg(integrator, true)
    alg.eigen_est === nothing ? maxeig!(integrator, cache) : alg.eigen_est(integrator)
    ccache.mdeg = Int(
        floor(
            sqrt(
                (
                    2 * abs(dt) * integrator.opts.internalnorm(
                        integrator.eigen_est, t
                    ) + 1.5
                ) / 0.811
            ) + 1
        )
    )
    ccache.mdeg = max(3, min(ccache.mdeg, 200)) - 2
    choose_deg!(integrator, cache)

    mdeg = ccache.mdeg
    start = ccache.start
    deg_index = ccache.deg_index
    σ = mσ[deg_index]
    τ = mτ[deg_index]

    sqrt_dt = sqrt(dt)
    μ = recf[start]  # here κ = 0
    tᵢ = t + dt * μ
    tᵢ₋₁ = tᵢ
    tᵢ₋₂ = t

    # stage 1
    @.. uᵢ₋₂ = uprev
    integrator.f(k, uprev, p, t)
    @.. uᵢ₋₁ = uprev + dt * μ * k

    # stages 2 upto s-2
    for i in 2:mdeg
        μ, κ = recf[start + 2 * (i - 2) + 1], recf[start + 2 * (i - 2) + 2]
        ν = 1 + κ
        integrator.f(k, uᵢ₋₁, p, t)
        @.. uᵢ = dt * μ * k + ν * uᵢ₋₁ - κ * uᵢ₋₂
        @.. uᵢ₋₂ = uᵢ₋₁
        @.. uᵢ₋₁ = uᵢ
        tᵢ = dt * μ + ν * tᵢ₋₁ - κ * tᵢ₋₂
        tᵢ₋₂ = tᵢ₋₁
        tᵢ₋₁ = tᵢ
    end

    #2 stage-finishing procedure
    #stage s-1
    integrator.f(k, uᵢ₋₁, p, tᵢ₋₁)
    @.. uᵢ₋₂ = uᵢ₋₁ + dt * σ * k
    @.. u = uᵢ₋₁ + dt * (σ - τ) * k
    tᵢ₋₂ = tᵢ₋₁ + dt * σ

    #stage s
    integrator.f(k, uᵢ₋₂, p, tᵢ₋₂)
    @.. u += dt * (σ + τ) * k

    if (W.dW isa Number) || (length(W.dW) == 1) || is_diagonal_noise(integrator.sol.prob)
        integrator.f.g(Gₛ, uᵢ₋₁, p, tᵢ₋₁)
        @.. u += Gₛ * W.dW

        @.. uᵢ₋₂ = uᵢ₋₁ + ((W.dW^2 - dt) / (2 * sqrt_dt)) * Gₛ
        integrator.f.g(Gₛ₁, uᵢ₋₂, p, tᵢ₋₁)
        @.. u += 1 // 2 * sqrt_dt * Gₛ₁

        @.. uᵢ₋₂ = uᵢ₋₁ - ((W.dW^2 - dt) / (2 * sqrt_dt)) * Gₛ
        integrator.f.g(Gₛ₁, uᵢ₋₂, p, tᵢ₋₁)
        @.. u -= 1 // 2 * sqrt_dt * Gₛ₁
    else
        integrator.f.g(Gₛ, uᵢ₋₁, p, tᵢ₋₁)
        mul!(uᵢ₋₂, Gₛ, W.dW)
        @.. u += uᵢ₋₂

        for i in 1:length(W.dW)
            WikRange .= 1 .* (1:length(W.dW) .== i)
            mul!(uᵢ, Gₛ, WikRange)

            @.. tmp = uᵢ₋₁ - 1 // 2 * sqrt_dt * uᵢ + 1 // 2 * (W.dW[i] / sqrt_dt) * uᵢ₋₂
            integrator.f.g(Gₛ₁, tmp, p, tᵢ₋₁)
            mul!(tmp, Gₛ₁, WikRange)
            @.. u += 1 // 2 * sqrt_dt * tmp

            @.. tmp = uᵢ₋₁ + 1 // 2 * sqrt_dt * uᵢ - 1 // 2 * (W.dW[i] / sqrt_dt) * uᵢ₋₂
            integrator.f.g(Gₛ₁, tmp, p, tᵢ₋₁)
            mul!(tmp, Gₛ₁, WikRange)
            @.. u -= 1 // 2 * sqrt_dt * tmp
        end
    end

    integrator.u = u
end
