# Generic implementation for purely implicit SDIRK methods (SFSDIRK4-8, Hairer4/42).
# Kept separate from generic_imex_perform_step.jl because those methods use ESDIRKIMEXTableau
# (2D Butcher arrays, explicit+implicit parts, dense-output r-matrices, step_limiter!),
# whereas pure SDIRK methods have only implicit stages, no explicit component, and no
# interpolation coefficients in the perform_step path. Merging the two would require
# special-casing every IMEX-specific feature in this simpler context.
struct PureSDIRKTableau{T, T2}
    γ::T2
    A::Vector{Vector{T}}
    c::Vector{T2}
    b::Vector{T}
    btilde::Vector{T}
    α::Vector{Vector{T}}
    s::Int
end

function _sfsdirk_α(::Type{T}, s::Int) where {T}
    α = Vector{Vector{T}}(undef, s)
    α[1] = T[]
    α[2] = T[]
    if s >= 3
        α[3] = T[one(T), zero(T)]
    end
    for i in 4:s
        αᵢ = zeros(T, i - 1)
        αᵢ[i - 1] = one(T)
        α[i] = αᵢ
    end
    return α
end

function PureSDIRKTableau(::SFSDIRK4, ::Type{T}, ::Type{T2}) where {T, T2}
    tab = SFSDIRK4Tableau(T, T2)
    s = 4
    A = [T[], T[tab.a21], T[tab.a31, tab.a32], T[tab.a41, tab.a42, tab.a43]]
    c = T2[tab.γ, tab.c2, tab.c3, tab.c4]
    b = T[tab.a51, tab.a52, tab.a53, tab.a54]
    return PureSDIRKTableau(tab.γ, A, c, b, T[], _sfsdirk_α(T, s), s)
end

function PureSDIRKTableau(::SFSDIRK5, ::Type{T}, ::Type{T2}) where {T, T2}
    tab = SFSDIRK5Tableau(T, T2)
    s = 5
    A = [T[], T[tab.a21], T[tab.a31, tab.a32], T[tab.a41, tab.a42, tab.a43],
        T[tab.a51, tab.a52, tab.a53, tab.a54]]
    c = T2[tab.γ, tab.c2, tab.c3, tab.c4, tab.c5]
    b = T[tab.a61, tab.a62, tab.a63, tab.a64, tab.a65]
    return PureSDIRKTableau(tab.γ, A, c, b, T[], _sfsdirk_α(T, s), s)
end

function PureSDIRKTableau(::SFSDIRK6, ::Type{T}, ::Type{T2}) where {T, T2}
    tab = SFSDIRK6Tableau(T, T2)
    s = 6
    A = [T[], T[tab.a21], T[tab.a31, tab.a32], T[tab.a41, tab.a42, tab.a43],
        T[tab.a51, tab.a52, tab.a53, tab.a54],
        T[tab.a61, tab.a62, tab.a63, tab.a64, tab.a65]]
    c = T2[tab.γ, tab.c2, tab.c3, tab.c4, tab.c5, tab.c6]
    b = T[tab.a71, tab.a72, tab.a73, tab.a74, tab.a75, tab.a76]
    return PureSDIRKTableau(tab.γ, A, c, b, T[], _sfsdirk_α(T, s), s)
end

function PureSDIRKTableau(::SFSDIRK7, ::Type{T}, ::Type{T2}) where {T, T2}
    tab = SFSDIRK7Tableau(T, T2)
    s = 7
    A = [T[], T[tab.a21], T[tab.a31, tab.a32], T[tab.a41, tab.a42, tab.a43],
        T[tab.a51, tab.a52, tab.a53, tab.a54],
        T[tab.a61, tab.a62, tab.a63, tab.a64, tab.a65],
        T[tab.a71, tab.a72, tab.a73, tab.a74, tab.a75, tab.a76]]
    c = T2[tab.γ, tab.c2, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7]
    b = T[tab.a81, tab.a82, tab.a83, tab.a84, tab.a85, tab.a86, tab.a87]
    return PureSDIRKTableau(tab.γ, A, c, b, T[], _sfsdirk_α(T, s), s)
end

function PureSDIRKTableau(::SFSDIRK8, ::Type{T}, ::Type{T2}) where {T, T2}
    tab = SFSDIRK8Tableau(T, T2)
    s = 8
    A = [T[], T[tab.a21], T[tab.a31, tab.a32], T[tab.a41, tab.a42, tab.a43],
        T[tab.a51, tab.a52, tab.a53, tab.a54],
        T[tab.a61, tab.a62, tab.a63, tab.a64, tab.a65],
        T[tab.a71, tab.a72, tab.a73, tab.a74, tab.a75, tab.a76],
        T[tab.a81, tab.a82, tab.a83, tab.a84, tab.a85, tab.a86, tab.a87]]
    c = T2[tab.γ, tab.c2, tab.c3, tab.c4, tab.c5, tab.c6, tab.c7, tab.c8]
    b = T[tab.a91, tab.a92, tab.a93, tab.a94, tab.a95, tab.a96, tab.a97, tab.a98]
    return PureSDIRKTableau(tab.γ, A, c, b, T[], _sfsdirk_α(T, s), s)
end

function PureSDIRKTableau(alg::Union{Hairer4, Hairer42}, ::Type{T}, ::Type{T2}) where {T, T2}
    tab = alg isa Hairer4 ? Hairer4Tableau(T, T2) : Hairer42Tableau(T, T2)
    s = 5
    A = [T[], T[tab.a21], T[tab.a31, tab.a32], T[tab.a41, tab.a42, tab.a43],
        T[tab.a51, tab.a52, tab.a53, tab.a54]]
    c = T2[tab.γ, tab.c2, tab.c3, tab.c4, one(T2)]
    b = T[tab.a51, tab.a52, tab.a53, tab.a54, T(tab.γ)]
    btilde = T[tab.btilde1, tab.btilde2, tab.btilde3, tab.btilde4, tab.btilde5]
    α = [T[], T[T(tab.α21)], T[T(tab.α31), T(tab.α32)],
        T[T(tab.α41), zero(T), T(tab.α43)],
        T[tab.bhat1, tab.bhat2, tab.bhat3, tab.bhat4]]
    return PureSDIRKTableau(tab.γ, A, c, b, btilde, α, s)
end

mutable struct PureSDIRKConstantCache{N, Tab} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

mutable struct PureSDIRKMutableCache{uType, rateType, uNoUnitsType, N, Tab} <:
    SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    zs::Vector{uType}
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
end

function full_cache(c::PureSDIRKMutableCache)
    return (c.u, c.uprev, c.fsalfirst, c.zs..., c.atmp)
end

const _PureSDIRKAlg = Union{
    OrdinaryDiffEqNewtonNonAdaptiveSDIRKAlgorithm,
    OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm,
}

function alg_cache(
        alg::_PureSDIRKAlg, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = PureSDIRKTableau(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, tab.γ, tab.c[1], Val(false), verbose
    )
    return PureSDIRKConstantCache(nlsolver, tab)
end

function alg_cache(
        alg::_PureSDIRKAlg, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = PureSDIRKTableau(alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, tab.γ, tab.c[1], Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)
    zs = [zero(u) for _ in 1:(tab.s - 1)]
    push!(zs, nlsolver.z)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    return PureSDIRKMutableCache(u, uprev, fsalfirst, zs, atmp, nlsolver, tab)
end

_smooth_est(alg::OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm) = alg.smooth_est
_smooth_est(alg) = false

function _pure_sdirk_predict_oop(αᵢ, zs, u)
    isempty(αᵢ) && return zero(u)
    zᵢ = zero(u)
    for j in eachindex(αᵢ)
        zᵢ = zᵢ + αᵢ[j] * zs[j]
    end
    return zᵢ
end

function _pure_sdirk_predict_iip!(zᵢ, αᵢ, zs)
    if isempty(αᵢ)
        zᵢ .= zero(eltype(zᵢ))
        return nothing
    end
    nonzero_idx = 0
    nonzero_count = 0
    for j in eachindex(αᵢ)
        if !iszero(αᵢ[j])
            nonzero_idx = j
            nonzero_count += 1
        end
    end
    if nonzero_count == 1 && αᵢ[nonzero_idx] == one(eltype(αᵢ))
        @.. broadcast=false zᵢ=zs[nonzero_idx]
        return nothing
    end
    zᵢ .= zero(eltype(zᵢ))
    for j in eachindex(αᵢ)
        @.. broadcast=false zᵢ=zᵢ + αᵢ[j] * zs[j]
    end
    return nothing
end

function _pure_sdirk_set_z1_iip!(zs, integrator, alg::OrdinaryDiffEqNewtonAdaptiveSDIRKAlgorithm)
    u, uprev = integrator.u, integrator.uprev
    dt = integrator.dt
    if integrator.success_iter > 0 && !integrator.reeval_fsal &&
            alg.extrapolant == :interpolant
        current_extrapolant!(u, integrator.t + dt, integrator)
        @.. broadcast=false zs[1]=u - uprev
    elseif alg.extrapolant == :linear
        @.. broadcast=false zs[1]=dt * integrator.fsalfirst
    else
        zs[1] .= zero(eltype(zs[1]))
    end
    return nothing
end

function _pure_sdirk_set_z1_iip!(zs, integrator, alg)
    zs[1] .= zero(eltype(zs[1]))
    return nothing
end

@muladd function perform_step!(
        integrator, cache::PureSDIRKConstantCache, repeat_step = false
    )
    (; dt, uprev, u, t) = integrator
    (; nlsolver) = cache
    (; A, c, b, btilde, α, s) = cache.tab
    alg = unwrap_alg(integrator, true)

    markfirststage!(nlsolver)

    zs = Vector{typeof(u)}(undef, s)
    zs[1] = zero(u)
    nlsolver.z = zs[1]
    nlsolver.c = c[1]
    nlsolver.tmp = uprev
    zs[1] = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    for i in 2:s
        zs[i] = _pure_sdirk_predict_oop(α[i], zs, u)
        nlsolver.z = zs[i]

        tmp = uprev
        Aᵢ = A[i]
        for j in eachindex(Aᵢ)
            tmp = tmp + Aᵢ[j] * zs[j]
        end
        nlsolver.tmp = tmp
        nlsolver.c = c[i]

        zs[i] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
    end

    u = uprev
    for i in 1:s
        u = u + b[i] * zs[i]
    end

    if integrator.opts.adaptive && !isempty(btilde)
        tmp = zero(u)
        for i in 1:s
            tmp = tmp + btilde[i] * zs[i]
        end
        if isnewton(nlsolver) && _smooth_est(alg)
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

    integrator.fsallast = zs[s] ./ dt
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u
end

@muladd function perform_step!(
        integrator, cache::PureSDIRKMutableCache, repeat_step = false
    )
    (; dt, uprev, u, t) = integrator
    (; zs, atmp, nlsolver) = cache
    (; tmp) = nlsolver
    (; A, c, b, btilde, α, s) = cache.tab
    alg = unwrap_alg(integrator, true)

    markfirststage!(nlsolver)

    _pure_sdirk_set_z1_iip!(zs, integrator, alg)
    nlsolver.z = zs[1]
    nlsolver.c = c[1]
    nlsolver.tmp = uprev
    zs[1] = nlsolve!(nlsolver, integrator, cache, repeat_step)
    nlsolvefail(nlsolver) && return

    for i in 2:s
        _pure_sdirk_predict_iip!(zs[i], α[i], zs)
        nlsolver.z = zs[i]

        @.. broadcast=false tmp=uprev
        Aᵢ = A[i]
        for j in eachindex(Aᵢ)
            @.. broadcast=false tmp=tmp + Aᵢ[j] * zs[j]
        end
        nlsolver.tmp = tmp
        nlsolver.c = c[i]

        if i == 2
            isnewton(nlsolver) && set_new_W!(nlsolver, false)
        end

        zs[i] = nlsolve!(nlsolver, integrator, cache, repeat_step)
        nlsolvefail(nlsolver) && return
    end

    @.. broadcast=false u=uprev
    for i in 1:s
        @.. broadcast=false u=u + b[i] * zs[i]
    end

    if integrator.opts.adaptive && !isempty(btilde)
        @.. broadcast=false tmp=zero(eltype(tmp))
        for i in 1:s
            @.. broadcast=false tmp=tmp + btilde[i] * zs[i]
        end
        if _smooth_est(alg) && isnewton(nlsolver)
            est = nlsolver.cache.dz
            dolinsolve(
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

    @.. broadcast=false integrator.fsallast=zs[s] / dt
end
