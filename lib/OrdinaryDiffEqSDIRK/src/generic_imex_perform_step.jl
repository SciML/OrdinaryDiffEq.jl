mutable struct ESDIRKIMEXConstantCache{Tab, N} <: OrdinaryDiffEqConstantCache
    nlsolver::N
    tab::Tab
end

mutable struct ESDIRKIMEXCache{uType, rateType, uNoUnitsType, N, Tab, kType, StepLimiter} <:
    SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    zs::Vector{uType}
    ks::Vector{kType}
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
    step_limiter!::StepLimiter
end

function full_cache(c::ESDIRKIMEXCache)
    base = (c.u, c.uprev, c.fsalfirst, c.zs..., c.atmp)
    if eltype(c.ks) !== Nothing
        return tuple(base..., c.ks...)
    end
    return base
end

function alg_cache(
        alg::OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = ESDIRKIMEXTableau(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ = tab.Ai[2, 2]
    c = tab.c[2]
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    return ESDIRKIMEXConstantCache(nlsolver, tab)
end

function alg_cache(
        alg::OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = ESDIRKIMEXTableau(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ = tab.Ai[2, 2]
    c = tab.c[2]
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    s = tab.s
    if f isa SplitFunction
        ks = [zero(u) for _ in 1:s]
    else
        ks = Vector{Nothing}(nothing, s)
    end

    zs = [zero(u) for _ in 1:(s - 1)]
    push!(zs, nlsolver.z)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    return ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, atmp, nlsolver, tab, alg.step_limiter!
    )
end

function initialize!(integrator, cache::ESDIRKIMEXConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return nothing
end

function initialize!(integrator, cache::ESDIRKIMEXCache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    return nothing
end

@muladd function perform_step!(
        integrator, cache::ESDIRKIMEXConstantCache, repeat_step = false
    )
    (; t, dt, uprev, u, p) = integrator
    nlsolver = cache.nlsolver
    tab = cache.tab
    (; Ai, bi, Ae, be, c, btilde, ebtilde, α, s) = tab
    alg = unwrap_alg(integrator, true)
    γ = Ai[2, 2]

    f2 = nothing
    k = Vector{typeof(u)}(undef, s)
    if integrator.f isa SplitFunction
        f_impl = integrator.f.f1
        f2 = integrator.f.f2
    else
        f_impl = integrator.f
    end

    markfirststage!(nlsolver)

    z = Vector{typeof(u)}(undef, s)

    if integrator.f isa SplitFunction
        z[1] = dt * f_impl(uprev, p, t)
    else
        z[1] = dt * integrator.fsalfirst
    end

    if integrator.f isa SplitFunction
        k[1] = dt * integrator.fsalfirst - z[1]
    end

    for i in 2:s
        tmp = uprev
        for j in 1:(i - 1)
            tmp = tmp + Ai[i, j] * z[j]
        end

        if integrator.f isa SplitFunction
            for j in 1:(i - 1)
                tmp = tmp + Ae[i, j] * k[j]
            end
        end

        if integrator.f isa SplitFunction
            z_guess = z[1]
        elseif α !== nothing && !iszero(α[i, 1])
            z_guess = zero(u)
            for j in 1:(i - 1)
                z_guess = z_guess + α[i, j] * z[j]
            end
        else
            z_guess = zero(u)
        end

        nlsolver.z = z_guess
        nlsolver.tmp = tmp
        nlsolver.c = c[i]
        nlsolver.γ = γ
        z[i] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return

        if integrator.f isa SplitFunction && i < s
            u_stage = tmp + γ * z[i]
            k[i] = dt * f2(u_stage, p, t + c[i] * dt)
            integrator.stats.nf2 += 1
        end
    end

    u = nlsolver.tmp + γ * z[s]
    if integrator.f isa SplitFunction
        k[s] = dt * f2(u, p, t + dt)
        integrator.stats.nf2 += 1
        u = uprev
        for i in 1:s
            u = u + bi[i] * z[i] + be[i] * k[i]
        end
    end

    if integrator.opts.adaptive
        tmp = zero(u)
        for i in 1:s
            tmp = tmp + btilde[i] * z[i]
        end
        if integrator.f isa SplitFunction && ebtilde !== nothing
            for i in 1:s
                tmp = tmp + ebtilde[i] * k[i]
            end
        end
        if isnewton(nlsolver) && alg.smooth_est
            integrator.stats.nsolve += 1
            est = _reshape(get_W(nlsolver) \ _vec(tmp), axes(tmp))
        else
            est = tmp
        end
        atmp = calculate_residuals(
            est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.f isa SplitFunction
        integrator.k[1] = integrator.fsalfirst
        integrator.fsallast = integrator.f(u, p, t + dt)
        integrator.k[2] = integrator.fsallast
    else
        integrator.fsallast = z[s] ./ dt
        integrator.k[1] = integrator.fsalfirst
        integrator.k[2] = integrator.fsallast
    end
    integrator.u = u
end

@muladd function perform_step!(integrator, cache::ESDIRKIMEXCache, repeat_step = false)
    (; t, dt, uprev, u, p) = integrator
    (; zs, ks, atmp, nlsolver, step_limiter!) = cache
    (; tmp) = nlsolver
    tab = cache.tab
    (; Ai, bi, Ae, be, c, btilde, ebtilde, α, s) = tab
    alg = unwrap_alg(integrator, true)
    γ = Ai[2, 2]

    f2 = nothing
    if integrator.f isa SplitFunction
        f_impl = integrator.f.f1
        f2 = integrator.f.f2
    else
        f_impl = integrator.f
    end

    markfirststage!(nlsolver)

    if integrator.f isa SplitFunction && !repeat_step && !integrator.last_stepfail
        f_impl(zs[1], integrator.uprev, p, integrator.t)
        zs[1] .*= dt
    else
        @..zs[1] = dt * integrator.fsalfirst
    end

    if integrator.f isa SplitFunction
        @..ks[1] = dt * integrator.fsalfirst - zs[1]
    end

    for i in 2:s
        @..tmp = uprev
        for j in 1:(i - 1)
            @..tmp += Ai[i, j] * zs[j]
        end

        if integrator.f isa SplitFunction
            for j in 1:(i - 1)
                @..tmp += Ae[i, j] * ks[j]
            end
        end

        if integrator.f isa SplitFunction
            copyto!(zs[i], zs[1])
        elseif α !== nothing && !iszero(α[i, 1])
            fill!(zs[i], zero(eltype(u)))
            for j in 1:(i - 1)
                @..zs[i] += α[i, j] * zs[j]
            end
        else
            fill!(zs[i], zero(eltype(u)))
        end

        nlsolver.z = zs[i]
        nlsolver.c = c[i]
        nlsolver.γ = γ
        zs[i] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
        if i > 2
            isnewton(nlsolver) && set_new_W!(nlsolver, false)
        end

        if integrator.f isa SplitFunction && i < s
            @..u = tmp + γ * zs[i]
            f2(ks[i], u, p, t + c[i] * dt)
            ks[i] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    @..u = tmp + γ * zs[s]
    if integrator.f isa SplitFunction
        f2(ks[s], u, p, t + dt)
        ks[s] .*= dt
        integrator.stats.nf2 += 1
        @..u = uprev
        for i in 1:s
            @..u += bi[i] * zs[i] + be[i] * ks[i]
        end
    end

    step_limiter!(u, integrator, p, t + dt)

    if integrator.opts.adaptive
        @..tmp = zero(eltype(u))
        for i in 1:s
            @..tmp += btilde[i] * zs[i]
        end
        if integrator.f isa SplitFunction && ebtilde !== nothing
            for i in 1:s
                @..tmp += ebtilde[i] * ks[i]
            end
        end
        if isnewton(nlsolver) && alg.smooth_est
            est = nlsolver.cache.dz
            linres = dolinsolve(
                integrator, nlsolver.cache.linsolve; b = _vec(tmp),
                linu = _vec(est)
            )
            integrator.stats.nsolve += 1
        else
            est = tmp
        end
        calculate_residuals!(
            atmp, est, uprev, u, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t
        )
        integrator.EEst = integrator.opts.internalnorm(atmp, t)
    end

    if integrator.f isa SplitFunction
        integrator.f(integrator.fsallast, u, p, t + dt)
    else
        @..integrator.fsallast = zs[s] / dt
    end
end
