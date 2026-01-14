abstract type SSPRKMutableCache <: OrdinaryDiffEqMutableCache end
abstract type SSPRKConstantCache <: OrdinaryDiffEqConstantCache end
get_fsalfirstlast(cache::SSPRKMutableCache, u) = (cache.fsalfirst, cache.k)

@cache struct SSPRK22Cache{uType, rateType, StageLimiter, StepLimiter, Thread} <:
    SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK22ConstantCache <: SSPRKConstantCache end

function alg_cache(
        alg::SSPRK22, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    return SSPRK22Cache(u, uprev, k, fsalfirst, alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

function alg_cache(
        alg::SSPRK22, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK22ConstantCache()
end

@cache struct SSPRK33Cache{uType, rateType, StageLimiter, StepLimiter, Thread} <:
    SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK33ConstantCache <: SSPRKConstantCache end

function alg_cache(
        alg::SSPRK33, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    return SSPRK33Cache(u, uprev, k, fsalfirst, alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

function alg_cache(
        alg::SSPRK33, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK33ConstantCache()
end

@cache struct KYKSSPRK42Cache{
        uType,
        rateType,
        TabType,
        StageLimiter,
        StepLimiter,
        Thread,
    } <: SSPRKMutableCache
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

struct KYKSSPRK42ConstantCache{T, T2} <: SSPRKConstantCache
    α20::T
    α21::T
    α30::T
    α32::T
    α40::T
    α43::T
    β10::T
    β21::T
    β30::T
    β32::T
    β40::T
    β43::T
    c1::T2
    c2::T2
    c3::T2
end

function KYKSSPRK42ConstantCache(T, T2)
    α20 = T(0.394806441339829)
    α21 = T(0.605193558660171)
    α30 = T(0.00279730708739)
    α32 = T(0.99720269291261)
    α40 = T(0.252860909354373)
    α43 = T(0.747139090645627)
    β10 = T(0.406584463657504)
    β21 = T(0.246062298456822)
    β30 = T(0.013637216641451)
    β32 = T(0.405447122055692)
    β40 = T(0.016453567333598)
    β43 = T(0.303775146447707)
    c1 = T2(0.406584463657504)
    c2 = T2(0.4921245969136438)
    c3 = T2(0.9098323119879613)
    return KYKSSPRK42ConstantCache(
        α20, α21, α30, α32, α40, α43, β10, β21, β30, β32, β40, β43, c1,
        c2, c3
    )
end

function alg_cache(
        alg::KYKSSPRK42, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    tab = KYKSSPRK42ConstantCache(
        constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits)
    )
    return KYKSSPRK42Cache(
        u, uprev, k, tmp, fsalfirst, tab, alg.stage_limiter!, alg.step_limiter!,
        alg.thread
    )
end

function alg_cache(
        alg::KYKSSPRK42, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return KYKSSPRK42ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct SSPRK53Cache{uType, rateType, TabType, StageLimiter, StepLimiter, Thread} <:
    SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    tmp::uType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK53ConstantCache{T, T2} <: SSPRKConstantCache
    α30::T
    α32::T
    α40::T
    α43::T
    α52::T
    α54::T
    β10::T
    β21::T
    β32::T
    β43::T
    β54::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2

    function SSPRK53ConstantCache(T, T2)
        α30 = T(0.355909775063327)
        α32 = T(0.644090224936674)
        α40 = T(0.367933791638137)
        α43 = T(0.632066208361863)
        α52 = T(0.237593836598569)
        α54 = T(0.762406163401431)
        β10 = T(0.377268915331368)
        β21 = T(0.377268915331368)
        β32 = T(0.242995220537396)
        β43 = T(0.23845893284629)
        β54 = T(0.287632146308408)
        c1 = T2(0.377268915331368)
        c2 = T2(0.754537830662736)
        c3 = T2(0.728985661612188)
        c4 = T2(0.69922613593167)

        return new{T, T2}(α30, α32, α40, α43, α52, α54, β10, β21, β32, β43, β54, c1, c2, c3, c4)
    end
end

function alg_cache(
        alg::SSPRK53, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    tab = SSPRK53ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return SSPRK53Cache(
        u, uprev, k, fsalfirst, tmp, tab, alg.stage_limiter!, alg.step_limiter!,
        alg.thread
    )
end

function alg_cache(
        alg::SSPRK53, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK53ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct SSPRK53_2N1Cache{
        uType,
        rateType,
        TabType,
        StageLimiter,
        StepLimiter,
        Thread,
    } <: SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK53_2N1ConstantCache{T, T2} <: SSPRKConstantCache
    α40::T
    α43::T
    β10::T
    β21::T
    β32::T
    β43::T
    β54::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2

    function SSPRK53_2N1ConstantCache(T, T2)
        α40 = T(0.571403511494104)
        α43 = T(0.428596488505896)
        β10 = T(0.443568244942995)
        β21 = T(0.291111420073766)
        β32 = T(0.270612601278217)
        β43 = T(0.110577759392786)
        β54 = T(0.458557505351052)
        c1 = T2(0.443568244942995)
        c2 = T2(0.734679665016762)
        c3 = T2(1.005292266294979)
        c4 = T2(0.541442494648948)

        return new{T, T2}(α40, α43, β10, β21, β32, β43, β54, c1, c2, c3, c4)
    end
end

function alg_cache(
        alg::SSPRK53_2N1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    tab = SSPRK53_2N1ConstantCache(
        constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits)
    )
    return SSPRK53_2N1Cache(
        u, uprev, k, fsalfirst, tab, alg.stage_limiter!, alg.step_limiter!,
        alg.thread
    )
end

function alg_cache(
        alg::SSPRK53_2N1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK53_2N1ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct SSPRK53_2N2Cache{
        uType,
        rateType,
        TabType,
        StageLimiter,
        StepLimiter,
        Thread,
    } <: SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK53_2N2ConstantCache{T, T2} <: SSPRKConstantCache
    α30::T
    α32::T
    α50::T
    α54::T
    β10::T
    β21::T
    β32::T
    β43::T
    β54::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2

    function SSPRK53_2N2ConstantCache(T, T2)
        α30 = T(0.682342861037239)
        α32 = T(0.317657138962761)
        α50 = T(0.0452309744824)
        α54 = T(0.9547690255176)
        β10 = T(0.465388589249323)
        β21 = T(0.465388589249323)
        β32 = T(0.124745797313998)
        β43 = T(0.465388589249323)
        β54 = T(0.154263303748666)
        c1 = T2(0.465388589249323)
        c2 = T2(0.930777178498646)
        c3 = T2(0.42041381284771)
        c4 = T2(0.885802402097033)

        return new{T, T2}(α30, α32, α50, α54, β10, β21, β32, β43, β54, c1, c2, c3, c4)
    end
end

function alg_cache(
        alg::SSPRK53_2N2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    tab = SSPRK53_2N2ConstantCache(
        constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits)
    )
    return SSPRK53_2N2Cache(
        u, uprev, k, fsalfirst, tab, alg.stage_limiter!, alg.step_limiter!,
        alg.thread
    )
end

function alg_cache(
        alg::SSPRK53_2N2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK53_2N2ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct SSPRK53_HCache{uType, rateType, TabType, StageLimiter, StepLimiter, Thread} <:
    SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    tmp::uType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK53_HConstantCache{T, T2} <: SSPRKConstantCache
    α30::T
    α32::T
    α40::T
    α41::T
    α43::T
    β10::T
    β21::T
    β32::T
    β43::T
    β54::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2

    function SSPRK53_HConstantCache(T, T2)
        α30 = T(0.308684154602513)
        α32 = T(0.691315845397487)
        α40 = T(0.280514990468574)
        α41 = T(0.270513101776498)
        α43 = T(0.448971907754928)
        β10 = T(0.377268915331368)
        β21 = T(0.377268915331368)
        β32 = T(0.260811979144498)
        β43 = T(0.169383144652957)
        β54 = T(0.377268915331368)
        c1 = T2(0.377268915331368)
        c2 = T2(0.754537830662737)
        c3 = T2(0.782435937433493)
        c4 = T2(0.622731084668631)

        return new{T, T2}(α30, α32, α40, α41, α43, β10, β21, β32, β43, β54, c1, c2, c3, c4)
    end
end

function alg_cache(
        alg::SSPRK53_H, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    tab = SSPRK53_HConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return SSPRK53_HCache(
        u, uprev, k, fsalfirst, tmp, tab, alg.stage_limiter!, alg.step_limiter!,
        alg.thread
    )
end

function alg_cache(
        alg::SSPRK53_H, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK53_HConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct SSPRK63Cache{uType, rateType, TabType, StageLimiter, StepLimiter, Thread} <:
    SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    tmp::uType
    u₂::uType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK63ConstantCache{T, T2} <: SSPRKConstantCache
    α40::T
    α41::T
    α43::T
    α62::T
    α65::T
    β10::T
    β21::T
    β32::T
    β43::T
    β54::T
    β65::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2

    function SSPRK63ConstantCache(T, T2)
        α40 = T(0.476769811285196)
        α41 = T(0.098511733286064)
        α43 = T(0.42471845542874)
        α62 = T(0.155221702560091)
        α65 = T(0.844778297439909)
        β10 = T(0.284220721334261)
        β21 = T(0.284220721334261)
        β32 = T(0.284220721334261)
        β43 = T(0.12071378576593)
        β54 = T(0.284220721334261)
        β65 = T(0.2401034970659)
        c1 = T2(0.284220721334261)
        c2 = T2(0.568441442668522)
        c3 = T2(0.852662164002783)
        c4 = T2(0.510854218958172)
        c5 = T2(0.795074940292433)

        return new{T, T2}(
            α40, α41, α43, α62, α65, β10, β21, β32, β43, β54, β65, c1, c2, c3, c4,
            c5
        )
    end
end

function alg_cache(
        alg::SSPRK63, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    u₂ = zero(u)
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    tab = SSPRK63ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return SSPRK63Cache(
        u, uprev, k, fsalfirst, tmp, u₂, tab, alg.stage_limiter!,
        alg.step_limiter!, alg.thread
    )
end

function alg_cache(
        alg::SSPRK63, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK63ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct SSPRK73Cache{uType, rateType, TabType, StageLimiter, StepLimiter, Thread} <:
    SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    tmp::uType
    u₁::uType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK73ConstantCache{T, T2} <: SSPRKConstantCache
    α40::T
    α43::T
    α50::T
    α51::T
    α54::T
    α73::T
    α76::T
    β10::T
    β21::T
    β32::T
    β43::T
    β54::T
    β65::T
    β76::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2

    function SSPRK73ConstantCache(T, T2)
        α40 = T(0.184962588071072)
        α43 = T(0.815037411928928)
        α50 = T(0.18071865657038)
        α51 = T(0.314831034403793)
        α54 = T(0.504450309025826)
        α73 = T(0.120199)
        α76 = T(0.879801)
        β10 = T(0.233213863663009)
        β21 = T(0.233213863663009)
        β32 = T(0.233213863663009)
        β43 = T(0.190078023865845)
        β54 = T(0.117644805593912)
        β65 = T(0.233213863663009)
        β76 = T(0.205181790464579)
        c1 = T2(0.233213863663009)
        c2 = T2(0.466427727326018)
        c3 = T2(0.699641590989027)
        c4 = T2(0.760312095463379)
        c5 = T2(0.574607439040817)
        c6 = T2(0.807821302703826)

        return new{T, T2}(
            α40, α43, α50, α51, α54, α73, α76, β10, β21, β32, β43, β54, β65, β76, c1,
            c2, c3, c4, c5, c6
        )
    end
end

function alg_cache(
        alg::SSPRK73, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    u₁ = zero(u)
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    tab = SSPRK73ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return SSPRK73Cache(
        u, uprev, k, fsalfirst, tmp, u₁, tab, alg.stage_limiter!,
        alg.step_limiter!, alg.thread
    )
end

function alg_cache(
        alg::SSPRK73, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK73ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct SSPRK83Cache{uType, rateType, TabType, StageLimiter, StepLimiter, Thread} <:
    SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    tmp::uType
    u₂::uType
    u₃::uType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK83ConstantCache{T, T2} <: SSPRKConstantCache
    α50::T
    α51::T
    α54::T
    α61::T
    α65::T
    α72::T
    α73::T
    α76::T
    β10::T
    β21::T
    β32::T
    β43::T
    β54::T
    β65::T
    β76::T
    β87::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2
    c5::T2
    c6::T2
    c7::T2

    function SSPRK83ConstantCache(T, T2)
        α50 = T(0.421366967085359)
        α51 = T(0.005949401107575)
        α54 = T(0.572683631807067)
        α61 = T(0.004254010666365)
        α65 = T(0.995745989333635)
        α72 = T(0.104380143093325)
        α73 = T(0.243265240906726)
        α76 = T(0.65235461599995)
        β10 = T(0.195804015330143)
        β21 = T(0.195804015330143)
        β32 = T(0.195804015330143)
        β43 = T(0.195804015330143)
        β54 = T(0.112133754621673)
        β65 = T(0.194971062960412)
        β76 = T(0.127733653231944)
        β87 = T(0.195804015330143)
        c1 = T2(0.195804015330143)
        c2 = T2(0.391608030660286)
        c3 = T2(0.587412045990429)
        c4 = T2(0.783216061320572)
        c5 = T2(0.561833689734037)
        c6 = T2(0.755247658555329)
        c7 = T2(0.804195984669857)

        return new{T, T2}(
            α50, α51, α54, α61, α65, α72, α73, α76, β10, β21, β32, β43, β54, β65,
            β76, β87, c1, c2, c3, c4, c5, c6, c7
        )
    end
end

function alg_cache(
        alg::SSPRK83, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    u₂ = zero(u)
    u₃ = zero(u)
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    tab = SSPRK83ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return SSPRK83Cache(
        u, uprev, k, fsalfirst, tmp, u₂, u₃, tab, alg.stage_limiter!,
        alg.step_limiter!, alg.thread
    )
end

function alg_cache(
        alg::SSPRK83, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK83ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct SSPRK43Cache{
        uType, rateType, uNoUnitsType, TabType, StageLimiter,
        StepLimiter, Thread,
    } <: SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    utilde::uType
    atmp::uNoUnitsType
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK43ConstantCache{T, T2} <: SSPRKConstantCache
    one_third_u::T
    two_thirds_u::T
    half_u::T
    half_t::T2

    function SSPRK43ConstantCache(T, T2)
        one_third_u = inv(T(3))
        two_thirds_u = 2 * one_third_u
        half_u = T(0.5)
        half_t = T2(0.5)

        return new{T, T2}(one_third_u, two_thirds_u, half_u, half_t)
    end
end

function alg_cache(
        alg::SSPRK43, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tab = SSPRK43ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return SSPRK43Cache(
        u, uprev, k, fsalfirst, utilde, atmp, tab, alg.stage_limiter!,
        alg.step_limiter!, alg.thread
    )
end

function alg_cache(
        alg::SSPRK43, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK43ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct SSPRK432Cache{
        uType,
        rateType,
        uNoUnitsType,
        StageLimiter,
        StepLimiter,
        Thread,
    } <: SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    utilde::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK432ConstantCache <: SSPRKConstantCache end

function alg_cache(
        alg::SSPRK432, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    return SSPRK432Cache(
        u, uprev, k, fsalfirst, utilde, atmp, alg.stage_limiter!,
        alg.step_limiter!, alg.thread
    )
end

function alg_cache(
        alg::SSPRK432, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK432ConstantCache()
end

@cache mutable struct SSPRKMSVS32Cache{
        uType, rateType, dtArrayType, dtType, StageLimiter,
        StepLimiter, Thread,
    } <: SSPRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    u_2::uType
    u_1::uType
    k::rateType
    tmp::uType
    dts::dtArrayType
    dtf::dtArrayType
    μ::dtType
    v_n::Float64
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
    step::Int
end

@cache mutable struct SSPRKMSVS32ConstantCache{uType, dtArrayType, dtType} <:
    OrdinaryDiffEqConstantCache
    u_2::uType
    u_1::uType
    dts::dtArrayType
    dtf::dtArrayType
    μ::dtType
    v_n::Float64
    step::Int
end

function alg_cache(
        alg::SSPRKMSVS32, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    fsalfirst = zero(rate_prototype)
    dts = fill(zero(dt), 3)
    dtf = fill(zero(dt), 2)
    μ = zero(dt)
    u_2 = zero(u)
    u_1 = zero(u)
    k = zero(rate_prototype)
    tmp = zero(u)
    return SSPRKMSVS32Cache(
        u, uprev, fsalfirst, u_2, u_1, k, tmp, dts, dtf, μ, 0.5,
        alg.stage_limiter!, alg.step_limiter!, alg.thread, 1
    )
end

function alg_cache(
        alg::SSPRKMSVS32, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dts = fill(zero(dt), 3)
    dtf = fill(zero(dt), 2)
    μ = zero(dt)
    u_2 = u
    u_1 = u
    return SSPRKMSVS32ConstantCache(u_2, u_1, dts, dtf, μ, 0.5, 1)
end

@cache mutable struct SSPRKMSVS43Cache{
        uType,
        rateType,
        StageLimiter,
        StepLimiter,
        Thread,
    } <: SSPRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    u_3::uType
    u_2::uType
    u_1::uType
    k::rateType
    k1::rateType
    k2::rateType
    k3::rateType
    tmp::uType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
    step::Int
end

@cache mutable struct SSPRKMSVS43ConstantCache{uType, rateType} <:
    OrdinaryDiffEqConstantCache
    u_3::uType
    u_2::uType
    u_1::uType
    k1::rateType
    k2::rateType
    k3::rateType
    step::Int
end

function alg_cache(
        alg::SSPRKMSVS43, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    fsalfirst = zero(rate_prototype)
    u_3 = zero(u)
    u_2 = zero(u)
    u_1 = zero(u)
    k = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    tmp = zero(u)
    return SSPRKMSVS43Cache(
        u, uprev, fsalfirst, u_3, u_2, u_1, k, k1, k2, k3, tmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread, 1
    )
end

function alg_cache(
        alg::SSPRKMSVS43, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    u_3 = u
    u_2 = u
    u_1 = u
    k1 = rate_prototype
    k2 = rate_prototype
    k3 = rate_prototype
    return SSPRKMSVS43ConstantCache(u_3, u_2, u_1, k1, k2, k3, 1)
end

@cache struct SSPRK932Cache{
        uType,
        rateType,
        uNoUnitsType,
        StageLimiter,
        StepLimiter,
        Thread,
    } <: SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    utilde::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK932ConstantCache <: SSPRKConstantCache end

function alg_cache(
        alg::SSPRK932, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    return SSPRK932Cache(
        u, uprev, k, fsalfirst, utilde, atmp, alg.stage_limiter!,
        alg.step_limiter!, alg.thread
    )
end

function alg_cache(
        alg::SSPRK932, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK932ConstantCache()
end

@cache struct SSPRK54Cache{uType, rateType, TabType, StageLimiter, StepLimiter, Thread} <:
    SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    k₃::rateType
    u₂::uType
    u₃::uType
    tmp::uType # should be u₄, but tmp is needed for callbacks
    tab::TabType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK54ConstantCache{T, T2} <: SSPRKConstantCache
    β10::T
    α20::T
    α21::T
    β21::T
    α30::T
    α32::T
    β32::T
    α40::T
    α43::T
    β43::T
    α52::T
    α53::T
    β53::T
    α54::T
    β54::T
    c1::T2
    c2::T2
    c3::T2
    c4::T2

    function SSPRK54ConstantCache(T, T2)
        β10 = T(0.39175222657189)
        α20 = T(0.444370493651235)
        α21 = T(0.555629506348765)
        β21 = T(0.368410593050371)
        α30 = T(0.620101851488403)
        α32 = T(0.379898148511597)
        β32 = T(0.251891774271694)
        α40 = T(0.178079954393132)
        α43 = T(0.821920045606868)
        β43 = T(0.544974750228521)
        α52 = T(0.517231671970585)
        α53 = T(0.096059710526147)
        β53 = T(0.06369246866629)
        α54 = T(0.386708617503269)
        β54 = T(0.226007483236906)
        c1 = T2(0.39175222657189)
        c2 = T2(0.58607968931154)
        c3 = T2(0.4745423631214)
        c4 = T2(0.935010630967653)

        return new{T, T2}(
            β10, α20, α21, β21, α30, α32, β32, α40, α43, β43, α52, α53, β53, α54,
            β54, c1, c2, c3, c4
        )
    end
end

function alg_cache(
        alg::SSPRK54, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    u₂ = zero(u)
    u₃ = zero(u)
    tmp = zero(u)
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    k₃ = zero(rate_prototype)
    tab = SSPRK54ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return SSPRK54Cache(
        u, uprev, k, fsalfirst, k₃, u₂, u₃, tmp, tab, alg.stage_limiter!,
        alg.step_limiter!, alg.thread
    )
end

function alg_cache(
        alg::SSPRK54, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK54ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct SSPRK104Cache{uType, rateType, StageLimiter, StepLimiter, Thread} <:
    SSPRKMutableCache
    u::uType
    uprev::uType
    k::rateType
    fsalfirst::rateType
    k₄::rateType
    tmp::uType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

struct SSPRK104ConstantCache <: SSPRKConstantCache end

function alg_cache(
        alg::SSPRK104, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    k = zero(rate_prototype)
    if calck
        fsalfirst = zero(k)
    else
        fsalfirst = k
    end
    k₄ = zero(rate_prototype)
    return SSPRK104Cache(
        u, uprev, k, fsalfirst, k₄, tmp, alg.stage_limiter!, alg.step_limiter!,
        alg.thread
    )
end

function alg_cache(
        alg::SSPRK104, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SSPRK104ConstantCache()
end

@cache struct KYK2014DGSSPRK_3S2_Cache{
        uType, rateType, TabType, StageLimiter, StepLimiter,
        Thread,
    } <:
    SSPRKMutableCache
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
        return new{T, T2}(α_10, α_20, α_21, α_30, α_32, β_10, β_21, β_30, β_32, c_1, c_2)
    end
end

function alg_cache(
        alg::KYK2014DGSSPRK_3S2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    u_1 = zero(u)
    u_2 = zero(u)
    kk_1 = zero(rate_prototype)
    kk_2 = zero(rate_prototype)
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    tab = KYK2014DGSSPRK_3S2_ConstantCache(
        constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits)
    )
    return KYK2014DGSSPRK_3S2_Cache(
        u, uprev, k, fsalfirst, tab, u_1, u_2, kk_1, kk_2,
        alg.stage_limiter!, alg.step_limiter!, alg.thread
    )
end

function alg_cache(
        alg::KYK2014DGSSPRK_3S2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose = ODEVerbosity()
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return KYK2014DGSSPRK_3S2_ConstantCache(
        constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits)
    )
end
