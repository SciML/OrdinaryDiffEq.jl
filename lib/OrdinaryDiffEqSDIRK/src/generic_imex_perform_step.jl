mutable struct ESDIRKIMEXConstantCache{Tab, N} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

mutable struct ESDIRKIMEXCache{uType, rateType, N, Tab, kType, StepLimiter} <:
    SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    zs::Vector{uType}
    ks::Vector{kType}
    nlsolver::N
    tab::Tab
    step_limiter!::StepLimiter
end

function full_cache(c::ESDIRKIMEXCache)
    base = (c.u, c.uprev, c.fsalfirst, c.zs...)
    if eltype(c.ks) !== Nothing
        return tuple(base..., c.ks...)
    end
    return base
end

const ESDIRKIMEXAlgorithm = Union{ARS222, ARS232, ARS443, BHR553}

function alg_cache(
        alg::ESDIRKIMEXAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
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
        alg::ESDIRKIMEXAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
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

    return ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, nlsolver, tab, alg.step_limiter!
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
    (; Ai, bi, Ae, be, c, s) = tab
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

    # Stage 1: explicit (ESDIRK: a₁₁ = 0)
    if integrator.f isa SplitFunction
        z[1] = dt * f_impl(uprev, p, t)
    else
        z[1] = dt * integrator.fsalfirst
    end

    if integrator.f isa SplitFunction
        k[1] = dt * integrator.fsalfirst - z[1]
    end

    # Stages 2..s
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

    # Compute solution
    u = nlsolver.tmp + γ * z[s]
    if integrator.f isa SplitFunction
        k[s] = dt * f2(u, p, t + dt)
        integrator.stats.nf2 += 1
        u = uprev
        for i in 1:s
            u = u + bi[i] * z[i] + be[i] * k[i]
        end
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
    (; zs, ks, nlsolver, step_limiter!) = cache
    (; tmp) = nlsolver
    tab = cache.tab
    (; Ai, bi, Ae, be, c, s) = tab
    γ = Ai[2, 2]

    f2 = nothing
    if integrator.f isa SplitFunction
        f_impl = integrator.f.f1
        f2 = integrator.f.f2
    else
        f_impl = integrator.f
    end

    markfirststage!(nlsolver)

    # Stage 1: explicit (ESDIRK: a₁₁ = 0)
    if integrator.f isa SplitFunction && !repeat_step && !integrator.last_stepfail
        f_impl(zs[1], integrator.uprev, p, integrator.t)
        zs[1] .*= dt
    else
        @.. broadcast=false zs[1] = dt * integrator.fsalfirst
    end

    if integrator.f isa SplitFunction
        @.. broadcast=false ks[1] = dt * integrator.fsalfirst - zs[1]
    end

    # Stages 2..s
    for i in 2:s
        @.. broadcast=false tmp = uprev
        for j in 1:(i - 1)
            @.. broadcast=false tmp += Ai[i, j] * zs[j]
        end

        if integrator.f isa SplitFunction
            for j in 1:(i - 1)
                @.. broadcast=false tmp += Ae[i, j] * ks[j]
            end
        end

        if integrator.f isa SplitFunction
            copyto!(zs[i], zs[1])
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
            @.. broadcast=false u = tmp + γ * zs[i]
            f2(ks[i], u, p, t + c[i] * dt)
            ks[i] .*= dt
            integrator.stats.nf2 += 1
        end
    end

    # Compute solution
    @.. broadcast=false u = tmp + γ * zs[s]
    if integrator.f isa SplitFunction
        f2(ks[s], u, p, t + dt)
        ks[s] .*= dt
        integrator.stats.nf2 += 1
        @.. broadcast=false u = uprev
        for i in 1:s
            @.. broadcast=false u += bi[i] * zs[i] + be[i] * ks[i]
        end
    end

    step_limiter!(u, integrator, p, t + dt)

    if integrator.f isa SplitFunction
        integrator.f(integrator.fsallast, u, p, t + dt)
    else
        @.. broadcast=false integrator.fsallast = zs[s] / dt
    end
end
