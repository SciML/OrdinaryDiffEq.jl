@cache struct EulerCache{uType, rateType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
    fsalfirst::rateType
end

@cache struct SplitEulerCache{uType, rateType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
    fsalfirst::rateType
end

function alg_cache(alg::SplitEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    SplitEulerCache(u, uprev, zero(u), zero(rate_prototype), zero(rate_prototype))
end

struct SplitEulerConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SplitEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    SplitEulerConstantCache()
end

function alg_cache(alg::Euler, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    EulerCache(u, uprev, zero(u), zero(rate_prototype), zero(rate_prototype))
end

struct EulerConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::Euler, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    EulerConstantCache()
end

@cache struct HeunCache{uType, rateType, uNoUnitsType, StageLimiter, StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    atmp::uNoUnitsType
    k::rateType
    fsalfirst::rateType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

@cache struct RalstonCache{uType, rateType, uNoUnitsType, StageLimiter, StepLimiter, Thread
                           } <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    atmp::uNoUnitsType
    k::rateType
    fsalfirst::rateType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::Heun, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    HeunCache(u, uprev, zero(u), atmp, zero(rate_prototype),
              zero(rate_prototype), alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

function alg_cache(alg::Ralston, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    RalstonCache(u, uprev, zero(u), atmp, zero(rate_prototype),
                 zero(rate_prototype), alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

struct HeunConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::Heun, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    HeunConstantCache()
end

struct RalstonConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::Ralston, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    RalstonConstantCache()
end

@cache struct MidpointCache{uType, rateType, uNoUnitsType, StageLimiter, StepLimiter, Thread
                            } <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct MidpointConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::Midpoint, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    MidpointCache(u, uprev, k, tmp, atmp, fsalfirst, alg.stage_limiter!, alg.step_limiter!,
                  alg.thread)
end

function alg_cache(alg::Midpoint, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    MidpointConstantCache()
end

@cache struct RK4Cache{uType, rateType, uNoUnitsType, StageLimiter, StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k₂::rateType
    k₃::rateType
    k₄::rateType
    k::rateType
    tmp::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct RK4ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::RK4, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k₁ = zero(rate_prototype)
    k₂ = zero(rate_prototype)
    k₃ = zero(rate_prototype)
    k₄ = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    RK4Cache(u, uprev, k₁, k₂, k₃, k₄, k, tmp, atmp, alg.stage_limiter!, alg.step_limiter!,
             alg.thread)
end

function alg_cache(alg::RK4, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    RK4ConstantCache()
end

@cache struct BS3Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter, StepLimiter,
                       Thread} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::BS3, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = BS3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    BS3Cache(u, uprev, k1, k2, k3, k4, utilde, tmp, atmp, tab, alg.stage_limiter!,
             alg.step_limiter!, alg.thread)
end

function alg_cache(alg::BS3, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    BS3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct OwrenZen3Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter,
                             StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::OwrenZen3, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = OwrenZen3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    OwrenZen3Cache(u, uprev, k1, k2, k3, k4, utilde, tmp, atmp, tab, alg.stage_limiter!,
                   alg.step_limiter!, alg.thread)
end

function alg_cache(alg::OwrenZen3, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    OwrenZen3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct OwrenZen4Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter,
                             StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::OwrenZen4, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = OwrenZen4ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    OwrenZen4Cache(u, uprev, k1, k2, k3, k4, k5, k6, utilde, tmp, atmp, tab,
                   alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

function alg_cache(alg::OwrenZen4, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    OwrenZen4ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct OwrenZen5Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter,
                             StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::OwrenZen5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = OwrenZen5ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    OwrenZen5Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, utilde, tmp, atmp, tab,
                   alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

function alg_cache(alg::OwrenZen5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    OwrenZen5ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct BS5Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter, StepLimiter,
                       Thread} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::BS5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = BS5ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    BS5Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, utilde, tmp, atmp, tab,
             alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

function alg_cache(alg::BS5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    BS5ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct Tsit5Cache{uType, rateType, uNoUnitsType, StageLimiter, StepLimiter,
                         Thread} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end
if isdefined(Base, :Experimental) && isdefined(Base.Experimental, :silence!)
    Base.Experimental.silence!(Tsit5Cache)
end

@cache struct RK46NLCache{uType, rateType, TabType, StageLimiter, StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    tmp::uType
    fsalfirst::rateType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct RK46NLConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    α2::T
    α3::T
    α4::T
    α5::T
    α6::T
    β1::T
    β2::T
    β3::T
    β4::T
    β5::T
    β6::T
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2

    function RK46NLConstantCache(T, T2)
        α2 = T(-0.737101392796)
        α3 = T(-1.634740794343)
        α4 = T(-0.744739003780)
        α5 = T(-1.469897351522)
        α6 = T(-2.813971388035)
        β1 = T(0.032918605146)
        β2 = T(0.823256998200)
        β3 = T(0.381530948900)
        β4 = T(0.200092213184)
        β5 = T(1.718581042715)
        β6 = T(0.27)
        c2 = T2(0.032918605146)
        c3 = T2(0.249351723343)
        c4 = T2(0.466911705055)
        c5 = T2(0.582030414044)
        c6 = T2(0.847252983783)
        new{T, T2}(α2, α3, α4, α5, α6, β1, β2, β3, β4, β5, β6, c2, c3, c4, c5, c6)
    end
end

function alg_cache(alg::RK46NL, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    tab = RK46NLConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    RK46NLCache(u, uprev, k, tmp, fsalfirst, tab, alg.stage_limiter!, alg.step_limiter!,
                alg.thread)
end

function alg_cache(alg::RK46NL, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    RK46NLConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

function alg_cache(alg::Tsit5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    Tsit5Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, utilde, tmp, atmp,
               alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

TruncatedStacktraces.@truncate_stacktrace Tsit5Cache 1

function alg_cache(alg::Tsit5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    Tsit5ConstantCache()
end

@cache struct DP5Cache{uType, rateType, uNoUnitsType, StageLimiter, StepLimiter,
                       Thread} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    dense_tmp3::rateType
    dense_tmp4::rateType
    update::rateType
    bspl::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::DP5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = k2
    k7 = zero(rate_prototype) # This is FSAL'd to k1
    dense_tmp3 = k2
    dense_tmp4 = k5
    bspl = k3

    tmp = zero(u) # has to be separate for FSAL
    utilde = tmp

    if eltype(u) != uEltypeNoUnits || calck
        update = zero(rate_prototype)
        atmp = similar(u, uEltypeNoUnits)
        recursivefill!(atmp, false)
    else
        update = k7
        atmp = k3
    end

    cache = DP5Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, dense_tmp3, dense_tmp4, update,
                     bspl, utilde, tmp, atmp, alg.stage_limiter!, alg.step_limiter!,
                     alg.thread)
    cache
end

function alg_cache(alg::DP5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    DP5ConstantCache()
end

@cache struct Anas5Cache{uType, rateType, uNoUnitsType, TabType} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    tab::TabType
end

function alg_cache(alg::Anas5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Anas5ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    Anas5Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, utilde, tmp, atmp, tab)
end

function alg_cache(alg::Anas5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    Anas5ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct KYK2014DGSSPRK_3S2_Cache{uType, rateType, TabType, StageLimiter, StepLimiter,
                                       Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    tab::TabType
    #temporary values for Shu-Osher
    u_1::uType
    u_2::uType
    kk_1::rateType
    kk_2::rateType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct KYK2014DGSSPRK_3S2_ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    #These are not α and β for RK but for Shu-Osher
    #see top of page 317 in
    #Optimal Strong-Stability-Preserving Runge–Kutta Time Discretizations for
    #Discontinuous Garlekin Methods, Kubatko, Yaeger, Ketcheson 2014
    α_10::T
    α_20::T
    α_21::T
    α_30::T
    α_32::T
    β_10::T
    β_21::T
    β_30::T
    β_32::T
    #Shu-Osher is normally stated for autonomous systems, the times
    #are calculated by hand for this scheme
    c_1::T
    c_2::T

    function KYK2014DGSSPRK_3S2_ConstantCache(T, T2)
        α_10 = T(1.0)
        α_20 = T(0.087353119859156)
        α_21 = T(0.912646880140844)
        α_30 = T(0.344956917166841)
        α_32 = T(0.655043082833159)
        β_10 = T(0.528005024856522)
        β_21 = T(0.481882138633993)
        β_30 = T(0.022826837460491)
        β_32 = T(0.345866039233415)
        c_1 = β_10
        c_2 = α_21 * β_10 + β_21 # ==0.96376427726
        new{T, T2}(α_10, α_20, α_21, α_30, α_32, β_10, β_21, β_30, β_32, c_1, c_2)
    end
end

function alg_cache(alg::KYK2014DGSSPRK_3S2, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    u_1 = zero(u)
    u_2 = zero(u)
    kk_1 = zero(rate_prototype)
    kk_2 = zero(rate_prototype)
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    tab = KYK2014DGSSPRK_3S2_ConstantCache(constvalue(uBottomEltypeNoUnits),
                                           constvalue(tTypeNoUnits))
    KYK2014DGSSPRK_3S2_Cache(u, uprev, k, fsalfirst, tab, u_1, u_2, kk_1, kk_2,
                             alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

function alg_cache(alg::KYK2014DGSSPRK_3S2, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    KYK2014DGSSPRK_3S2_ConstantCache(constvalue(uBottomEltypeNoUnits),
                                     constvalue(tTypeNoUnits))
end

@cache struct RKO65Cache{uType, rateType, TabType, StageLimiter, StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    tmp::uType
    fsalfirst::rateType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

u_cache(c::RKO65Cache) = ()
du_cache(c::RKO65Cache) = (c.k, c.du, c.fsalfirst)

struct RKO65ConstantCache{T1, T2} <: OrdinaryDiffEqConstantCache
    α21::T1
    α31::T1
    α41::T1
    α51::T1

    α32::T1
    α42::T1
    α52::T1
    α62::T1

    α43::T1
    α53::T1
    α63::T1

    α54::T1
    α64::T1

    α65::T1

    β2::T1
    β3::T1
    β4::T1
    β5::T1
    β6::T1

    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2

    function RKO65ConstantCache(T1, T2)
        # elements of Butcher Table
        α21 = T1(1 // 6)
        α31 = T1(-15 // 8)
        α41 = T1(-9 // 1)
        α51 = T1(-3 // 1)

        α32 = T1(21 // 8)
        α42 = T1(75 // 7)
        α52 = T1(34257 // 8750)
        α62 = T1(123 // 380)

        α43 = T1(-5 // 7)
        α53 = T1(-114 // 875)
        α63 = T1(5 // 2)

        α54 = T1(19 // 1250)
        α64 = T1(3 // 20)

        α65 = T1(-75 // 38)

        β2 = T1(54 // 133)
        β3 = T1(32 // 21)
        β4 = T1(1 // 18)
        β5 = T1(-125 // 114)
        β6 = T1(1 // 9)

        c1 = T2(2 // 3)
        c2 = T2(1 // 6)
        c3 = T2(3 // 4)
        c4 = T2(1 // 1)
        c5 = T2(4 // 5)
        c6 = T2(1 // 1)
        new{T1, T2}(α21, α31, α41, α51, α32, α42, α52, α62, α43, α53, α63, α54, α64, α65,
                    β2, β3, β4, β5, β6, c1, c2, c3, c4, c5, c6)
    end
end

function alg_cache(alg::RKO65, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    RKO65ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits)) # why not real(tTypeNoUnits)?
end

function alg_cache(alg::RKO65, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)

    k = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)

    fsalfirst = zero(rate_prototype)

    tab = RKO65ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    RKO65Cache(u, uprev, k, k1, k2, k3, k4, k5, k6, tmp, fsalfirst, tab, alg.stage_limiter!,
               alg.step_limiter!, alg.thread)
end

@cache struct FRK65Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter, StepLimiter,
                         Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    utilde::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    k9::rateType
    tmp::uType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct FRK65ConstantCache{T1, T2} <: OrdinaryDiffEqConstantCache
    α21::T1
    α31::T1
    α41::T1
    α51::T1
    α61::T1
    α71::T1
    α81::T1
    α91::T1

    α32::T1

    α43::T1
    α53::T1
    α63::T1
    α73::T1
    α83::T1

    α54::T1
    α64::T1
    α74::T1
    α84::T1
    α94::T1

    α65::T1
    α75::T1
    α85::T1
    α95::T1

    α76::T1
    α86::T1
    α96::T1

    α87::T1
    α97::T1

    α98::T1

    β1::T1
    # β4::T1
    # β5::T1
    # β6::T1
    β7::T1
    β8::T1

    β1tilde::T1
    β4tilde::T1
    β5tilde::T1
    β6tilde::T1
    β7tilde::T1
    β8tilde::T1
    β9tilde::T1

    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2
    c7::T2
    c8::T2
    c9::T2

    d1::T1
    d2::T1
    d3::T1
    d4::T1
    d5::T1
    d6::T1
    d7::T1
    d8::T1
    d9::T1
    d10::T1
    d11::T1
    d12::T1
    d13::T1

    e1::T1
    e2::T1
    e3::T1
    e4::T1
    e5::T1
    e6::T1
    e7::T1
    e8::T1
    e9::T1
    e10::T1
    e11::T1

    f1::T1
    f2::T1
    f3::T1
    f4::T1
    f5::T1
    f6::T1
    f7::T1
    f8::T1
    f9::T1
    f10::T1
    f11::T1

    function FRK65ConstantCache(T1, T2)

        # elements of Butcher Table
        α21 = T1(1 // 89)
        α31 = T1(-38624 // 142129)
        α41 = T1(51 // 1508)
        α51 = T1(3259284578 // 3517556363)
        α61 = T1(-108363632681 // 45875676369)
        α71 = T1(7137368591 // 11299833148)
        α81 = T1(8898824396 // 9828950919)
        α91 = T1(1026331676 // 33222204855)

        α32 = T1(51442 // 142129)

        α43 = T1(153 // 1508)
        α53 = T1(-69727055112 // 19553806387)
        α63 = T1(80902506271 // 8700424616)
        α73 = T1(-33088067061 // 10572251159)
        α83 = T1(25673454973 // 11497947835)

        α54 = T1(36230363390 // 11788838981)
        α64 = T1(-120088218786 // 17139312481)
        α74 = T1(11481363823 // 3650030081)
        α84 = T1(-74239028301 // 15737704666)
        α94 = T1(1450675392 // 5936579813)

        α65 = T1(4533285649 // 6676940598)
        α75 = T1(-4096673444 // 7349814937)
        α85 = T1(222688842816 // 44196813415)
        α95 = T1(4617877550 // 16762182457)

        α76 = T1(9911918171 // 12847192605)
        α86 = T1(-105204445705 // 30575217706)
        α96 = T1(1144867463 // 6520294355)

        α87 = T1(8799291910 // 8966990271)
        α97 = T1(1822809703 // 7599996644)

        α98 = T1(79524953 // 2351253316)

        β1 = T1(1026331676 // 33222204855)
        #β4 = T1(1450675392//5936579813)
        #β5 = T1(4617877550//16762182457)
        #β6 = T1(1144867463//6520294355)
        β7 = T1(1822809703 // 7599996644)
        β8 = T1(79524953 // 2351253316)

        β1tilde = T1(413034411 // 13925408836)
        β4tilde = T1(1865954212 // 7538591735)
        β5tilde = T1(4451980162 // 16576017119)
        β6tilde = T1(1157843020 // 6320223511)
        β7tilde = T1(802708729 // 3404369569)
        β8tilde = T1(-251398161 // 17050111121)
        β9tilde = T1(1 // 20)

        c2 = T2(1 // 89)
        c3 = T2(34 // 377)
        c4 = T2(51 // 377)
        c5 = T2(14497158 // 33407747)
        c6 = T2(9744566553 // 16002998914)
        c7 = T2(330 // 383)
        c8 = T2(1 // 1)
        c9 = T2(1 // 1)

        d1 = T1(140209127 // 573775965)
        d2 = T1(-8530039 // 263747097)
        d3 = T1(-308551 // 104235790)
        d4 = T1(233511 // 333733259)
        d5 = T1(9126 // 184950985)
        d6 = T1(22 // 50434083)
        d7 = T1(19 // 424427471)
        d8 = T1(-28711 // 583216934059)
        d9 = T1(-3831531 // 316297807)
        d10 = T1(551767 // 187698280)
        d11 = T1(9205 // 210998423)
        d12 = T1(-250 // 519462673)
        d13 = T1(67 // 327513887)

        e1 = T1(437217689 // 1587032700)
        e2 = T1(-15824413 // 592362279)
        e3 = T1(-1563775 // 341846569)
        e4 = T1(270497 // 369611210)
        e5 = T1(-26623 // 453099487)
        e6 = T1(-616297487849)
        e7 = T1(-47682337 // 491732789)
        e8 = T1(-4778275 // 287766311)
        e9 = T1(641177 // 265376522)
        e10 = T1(44633 // 291742143)
        e11 = T1(611 // 223639880)

        f1 = T1(44861261 // 255495624)
        f2 = T1(-11270940 // 352635157)
        f3 = T1(-182222 // 232874507)
        f4 = T1(164263 // 307215200)
        f5 = T1(32184 // 652060417)
        f6 = T1(-352 // 171021903)
        f7 = T1(-18395427 // 101056291)
        f8 = T1(-621686 // 139501937)
        f9 = T1(2030024 // 612171255)
        f10 = T1(-711049 // 7105160932)
        f11 = T1(267 // 333462710)

        new{T1, T2}(α21, α31, α41, α51, α61, α71, α81, α91, α32, α43, α53, α63, α73, α83,
                    α54, α64, α74, α84, α94, α65, α75, α85, α95, α76, α86, α96, α87, α97,
                    α98, β1, β7, β8, β1tilde, β4tilde, β5tilde, β6tilde, β7tilde, β8tilde,
                    β9tilde, c2, c3, c4, c5, c6, c7, c8, c9, d1, d2, d3, d4, d5, d6, d7, d8,
                    d9, d10, d11, d12, d13, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11,
                    f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11)
    end
end

function alg_cache(alg::FRK65, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    FRK65ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

function alg_cache(alg::FRK65, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = FRK65ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    k9 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    FRK65Cache(u, uprev, utilde, k1, k2, k3, k4, k5, k6, k7, k8, k9, tmp, atmp, tab,
               alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

@cache struct RKMCache{uType, rateType, TabType, StageLimiter, StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k::rateType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    tmp::uType
    fsalfirst::rateType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct RKMConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
    α2::T
    α3::T
    α4::T
    α5::T
    α6::T
    β1::T
    β2::T
    β3::T
    β4::T
    β6::T
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2

    function RKMConstantCache(::Type{T}, ::Type{T2}) where {T, T2}
        # α2 = T(0.16791846623918)
        # α3 = T(0.48298439719700)
        # α4 = T(0.70546072965982)
        # α5 = T(0.09295870406537)
        # α6 = T(0.76210081248836)
        # β1 = T(-0.15108370762927)
        # β2 = T(0.75384683913851)
        # β3 = T(-0.36016595357907)
        # β4 = T(0.52696773139913)
        # β6 = T(0.23043509067071)
        # c2 = T2(0.16791846623918)
        # c3 = T2(0.48298439719700)
        # c4 = T2(0.70546072965982)
        # c5 = T2(0.09295870406537)
        # c6 = T2(0.76210081248836)

        α2 = T(0.167266187050662)
        α3 = T(0.484574582244783)
        α4 = T(0.536909403373491)
        α5 = T(0.082069535961948)
        α6 = T(0.853923000035347)
        β1 = T(-0.028289441132839)
        β2 = T(0.463968918564710)
        β3 = T(-0.434414348751899)
        β4 = T(0.693796229087598)
        β6 = T(0.304938642232430)
        c2 = T2(0.167266187050662)
        c3 = T2(0.484574582244783)
        c4 = T2(0.536909403373491)
        c5 = T2(0.082069535961948)
        c6 = T2(0.853923000035347)

        new{T, T2}(α2, α3, α4, α5, α6, β1, β2, β3, β4, β6, c2, c3, c4, c5, c6)
    end
end

function alg_cache(alg::RKM, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    RKMConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

function alg_cache(alg::RKM, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = RKMConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    tmp = zero(u)
    fsalfirst = zero(rate_prototype)
    RKMCache(u, uprev, k, k1, k2, k3, k4, k5, k6, tmp, fsalfirst, tab, alg.stage_limiter!,
             alg.step_limiter!, alg.thread)
end

@cache struct MSRK5Cache{uType, rateType, TabType, StageLimiter, StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    fsalfirst::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    k9::rateType
    k::rateType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::MSRK5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MSRK5ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

function alg_cache(alg::MSRK5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    k9 = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    fsalfirst = zero(u)
    tab = MSRK5ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    MSRK5Cache(u, uprev, tmp, fsalfirst, k1, k2, k3, k4, k5, k6, k7, k8, k9, k, tab,
               alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

@cache struct MSRK6Cache{uType, rateType, TabType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    fsalfirst::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    k9::rateType
    k::rateType
    tab::TabType
end

function alg_cache(alg::MSRK6, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MSRK6ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

function alg_cache(alg::MSRK6, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    k9 = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    fsalfirst = zero(u)
    tab = MSRK6ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    MSRK6Cache(u, uprev, tmp, fsalfirst, k1, k2, k3, k4, k5, k6, k7, k8, k9, k, tab)
end

@cache struct Stepanov5Cache{uType, rateType, TabType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    fsalfirst::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k::rateType
    tab::TabType
end

function alg_cache(alg::Stepanov5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Stepanov5ConstantCache(constvalue(uBottomEltypeNoUnits),
                                  constvalue(tTypeNoUnits))
end

function alg_cache(alg::Stepanov5, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    fsalfirst = zero(u)
    tab = Stepanov5ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    Stepanov5Cache(u, uprev, tmp, fsalfirst, k1, k2, k3, k4, k5, k6, k7, k, tab)
end

@cache struct SIR54Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter, StepLimiter,
                         Thread} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::SIR54, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SIR54ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

function alg_cache(alg::SIR54, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    tab = SIR54ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    SIR54Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, k8, utilde, tmp, atmp, tab,
               alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

@cache struct Alshina2Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter,
                            StepLimiter, Thread} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    utilde::uType
    k1::rateType
    k2::rateType
    atmp::uNoUnitsType
    tmp::uType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::Alshina2, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Alshina2ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

function alg_cache(alg::Alshina2, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    utilde = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tab = Alshina2ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    Alshina2Cache(u, uprev, utilde, k1, k2, atmp, tmp, tab, alg.stage_limiter!,
                  alg.step_limiter!, alg.thread)
end

@cache struct Alshina3Cache{uType, rateType, uNoUnitsType, TabType, StageLimiter,
                            StepLimiter, Thread} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    utilde::uType
    k1::rateType
    k2::rateType
    k3::rateType
    atmp::uNoUnitsType
    tmp::uType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::Alshina3, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Alshina3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

function alg_cache(alg::Alshina3, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    utilde = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tab = Alshina3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    Alshina3Cache(u, uprev, utilde, k1, k2, k3, atmp, tmp, tab, alg.stage_limiter!,
                  alg.step_limiter!, alg.thread)
end

@cache struct Alshina6Cache{uType, rateType, TabType, StageLimiter, StepLimiter, Thread} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    tmp::uType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(alg::Alshina6, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return Alshina6ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

function alg_cache(alg::Alshina6, u, rate_prototype, ::Type{uEltypeNoUnits},
                   ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
                   dt, reltol, p, calck,
                   ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    tmp = zero(u)
    tab = Alshina6ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    Alshina6Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, tmp, tab, alg.stage_limiter!,
                  alg.step_limiter!, alg.thread)
end
