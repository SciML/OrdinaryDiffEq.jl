mutable struct Kvaerno3ConstantCache{Tab, N} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

function alg_cache(alg::Kvaerno3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Kvaerno3Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, 2tab.γ
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))
    Kvaerno3ConstantCache(nlsolver, tab)
end

@cache mutable struct Kvaerno3Cache{uType, rateType, uNoUnitsType, Tab, N, StepLimiter} <:
                      SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    z₁::uType
    z₂::uType
    z₃::uType
    z₄::uType
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
    step_limiter!::StepLimiter
end

function alg_cache(alg::Kvaerno3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Kvaerno3Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, 2tab.γ
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))
    fsalfirst = zero(rate_prototype)

    z₁ = zero(u)
    z₂ = zero(u)
    z₃ = zero(u)
    z₄ = nlsolver.z
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    Kvaerno3Cache(
        u, uprev, fsalfirst, z₁, z₂, z₃, z₄, atmp, nlsolver, tab, alg.step_limiter!)
end

@cache mutable struct KenCarp3ConstantCache{N, Tab} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

function alg_cache(alg::KenCarp3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = KenCarp3Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))

    KenCarp3ConstantCache(nlsolver, tab)
end

@cache mutable struct KenCarp3Cache{
    uType, rateType, uNoUnitsType, N, Tab, kType, StepLimiter} <:
                      SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    z₁::uType
    z₂::uType
    z₃::uType
    z₄::uType
    ks::Vector{kType}
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
    step_limiter!::StepLimiter
end

function alg_cache(alg::KenCarp3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = KenCarp3Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))
    fsalfirst = zero(rate_prototype)

    if f isa SplitFunction
        ks = [zero(u) for _ in 1:4] 
    else
        ks = [nothing for _ in 1:4]
        uf = UJacobianWrapper(f, t, p)
    end

    z₁ = zero(u)
    z₂ = zero(u)
    z₃ = zero(u)
    z₄ = nlsolver.z
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    KenCarp3Cache(u, uprev, fsalfirst, z₁, z₂, z₃, z₄, ks, atmp, nlsolver, tab, alg.step_limiter!)
end

@cache mutable struct CFNLIRK3ConstantCache{N, Tab} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

function alg_cache(alg::CFNLIRK3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = CFNLIRK3Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))

    CFNLIRK3ConstantCache(nlsolver, tab)
end

@cache mutable struct CFNLIRK3Cache{uType, rateType, uNoUnitsType, N, Tab, kType} <:
                      SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    z₁::uType
    z₂::uType
    z₃::uType
    z₄::uType
    ks::Vector{kType}
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
end

function alg_cache(alg::CFNLIRK3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = CFNLIRK3Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))
    fsalfirst = zero(rate_prototype)

    ks = [zero(u) for _ in 1:4]

    z₁ = zero(u)
    z₂ = zero(u)
    z₃ = zero(u)
    z₄ = nlsolver.z
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    CFNLIRK3Cache(u, uprev, fsalfirst, z₁, z₂, z₃, z₄, ks, atmp, nlsolver, tab)
end

@cache mutable struct Kvaerno4ConstantCache{N, Tab} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

function alg_cache(alg::Kvaerno4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Kvaerno4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))
    Kvaerno4ConstantCache(nlsolver, tab)
end

@cache mutable struct Kvaerno4Cache{uType, rateType, uNoUnitsType, N, Tab, StepLimiter} <:
                      SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    z₁::uType
    z₂::uType
    z₃::uType
    z₄::uType
    z₅::uType
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
    step_limiter!::StepLimiter
end

function alg_cache(alg::Kvaerno4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Kvaerno4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))
    fsalfirst = zero(rate_prototype)

    z₁ = zero(u)
    z₂ = zero(u)
    z₃ = zero(u)
    z₄ = zero(u)
    z₅ = nlsolver.z
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    Kvaerno4Cache(
        u, uprev, fsalfirst, z₁, z₂, z₃, z₄, z₅, atmp, nlsolver, tab, alg.step_limiter!)
end

@cache mutable struct KenCarp4ConstantCache{N, Tab} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

function alg_cache(alg::KenCarp4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = KenCarp4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))
    KenCarp4ConstantCache(nlsolver, tab)
end

@cache mutable struct KenCarp4Cache{
    uType, rateType, uNoUnitsType, N, Tab, kType, StepLimiter} <:
                      SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    z₁::uType
    z₂::uType
    z₃::uType
    z₄::uType
    z₅::uType
    z₆::uType
    ks::Vector{kType}
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
    step_limiter!::StepLimiter
end

TruncatedStacktraces.@truncate_stacktrace KenCarp4Cache 1

function alg_cache(alg::KenCarp4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = KenCarp4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))
    fsalfirst = zero(rate_prototype)

    if f isa SplitFunction
        ks = [zero(u) for _ in 1:6]
    else
        ks = [nothing for _ in 1:6]
        uf = UJacobianWrapper(f, t, p)
    end

    z₁ = zero(u)
    z₂ = zero(u)
    z₃ = zero(u)
    z₄ = zero(u)
    z₅ = zero(u)
    z₆ = nlsolver.z
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    KenCarp4Cache(
        u, uprev, fsalfirst, z₁, z₂, z₃, z₄, z₅, z₆, ks, atmp,
        nlsolver, tab, alg.step_limiter!)
end

@cache mutable struct Kvaerno5ConstantCache{N, Tab} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

function alg_cache(alg::Kvaerno5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Kvaerno5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))

    Kvaerno5ConstantCache(nlsolver, tab)
end

@cache mutable struct Kvaerno5Cache{uType, rateType, uNoUnitsType, N, Tab, StepLimiter} <:
                      SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    z₁::uType
    z₂::uType
    z₃::uType
    z₄::uType
    z₅::uType
    z₆::uType
    z₇::uType
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
    step_limiter!::StepLimiter
end

function alg_cache(alg::Kvaerno5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Kvaerno5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))
    fsalfirst = zero(rate_prototype)

    z₁ = zero(u)
    z₂ = zero(u)
    z₃ = zero(u)
    z₄ = zero(u)
    z₅ = zero(u)
    z₆ = zero(u)
    z₇ = nlsolver.z
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    Kvaerno5Cache(u, uprev, fsalfirst, z₁, z₂, z₃, z₄, z₅, z₆,
        z₇, atmp, nlsolver, tab, alg.step_limiter!)
end

@cache mutable struct KenCarp5ConstantCache{N, Tab} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

function alg_cache(alg::KenCarp5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = KenCarp5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))

    KenCarp5ConstantCache(nlsolver, tab)
end

@cache mutable struct KenCarp5Cache{
    uType, rateType, uNoUnitsType, N, Tab, kType, StepLimiter} <:
                      SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    z₁::uType
    z₂::uType
    z₃::uType
    z₄::uType
    z₅::uType
    z₆::uType
    z₇::uType
    z₈::uType
    ks::Vector{kType}
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
    step_limiter!::StepLimiter
end

function alg_cache(alg::KenCarp5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = KenCarp5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))
    fsalfirst = zero(rate_prototype)

    if f isa SplitFunction
        ks = [zero(u) for _ in 1:8]
    else
        ks = [nothing for _ in 1:8]
    end

    z₁ = zero(u)
    z₂ = zero(u)
    z₃ = zero(u)
    z₄ = zero(u)
    z₅ = zero(u)
    z₆ = zero(u)
    z₇ = zero(u)
    z₈ = nlsolver.z
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    KenCarp5Cache(u, uprev, fsalfirst, z₁, z₂, z₃, z₄, z₅, z₆, z₇, z₈,
        ks, atmp, nlsolver, tab, alg.step_limiter!)
end

@cache mutable struct KenCarp47ConstantCache{N, Tab} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

function alg_cache(alg::KenCarp47, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = KenCarp47Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))

    KenCarp47ConstantCache(nlsolver, tab)
end

@cache mutable struct KenCarp47Cache{uType, rateType, uNoUnitsType, N, Tab, kType} <:
                      SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    z₁::uType
    z₂::uType
    z₃::uType
    z₄::uType
    z₅::uType
    z₆::uType
    z₇::uType
    ks::Vector{kType}
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
end
TruncatedStacktraces.@truncate_stacktrace KenCarp47Cache 1

function alg_cache(alg::KenCarp47, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = KenCarp47Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))
    fsalfirst = zero(rate_prototype)

    if f isa SplitFunction
        ks = [zero(u) for _ in 1:7]
    else
        ks = [nothing for _ in 1:7]
    end

    z₁ = zero(u)
    z₂ = zero(u)
    z₃ = zero(u)
    z₄ = zero(u)
    z₅ = zero(u)
    z₆ = zero(u)
    z₇ = nlsolver.z
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    KenCarp47Cache(u, uprev, fsalfirst, z₁, z₂, z₃, z₄, z₅, z₆, z₇,
        ks, atmp, nlsolver, tab)
end

@cache mutable struct KenCarp58ConstantCache{N, Tab} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

function alg_cache(alg::KenCarp58, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = KenCarp58Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))

    KenCarp58ConstantCache(nlsolver, tab)
end

@cache mutable struct KenCarp58Cache{uType, rateType, uNoUnitsType, N, Tab, kType} <:
                      SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    z₁::uType
    z₂::uType
    z₃::uType
    z₄::uType
    z₅::uType
    z₆::uType
    z₇::uType
    z₈::uType
    ks::Vector{kType}
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
end

TruncatedStacktraces.@truncate_stacktrace KenCarp58Cache 1

function alg_cache(alg::KenCarp58, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = KenCarp58Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))
    fsalfirst = zero(rate_prototype)

    if f isa SplitFunction
        ks = [zero(u) for _ in 1:8]
    else
        ks = [nothing for _ in 1:8]
    end

    z₁ = zero(u)
    z₂ = zero(u)
    z₃ = zero(u)
    z₄ = zero(u)
    z₅ = zero(u)
    z₆ = zero(u)
    z₇ = zero(u)
    z₈ = nlsolver.z
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    KenCarp58Cache(u, uprev, fsalfirst, z₁, z₂, z₃, z₄, z₅, z₆, z₇, z₈,
        ks, atmp, nlsolver, tab)
end
