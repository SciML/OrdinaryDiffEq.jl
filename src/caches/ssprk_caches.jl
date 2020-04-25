@cache struct SSPRK22Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType # superfluous, only needed for callbacks...
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

struct SSPRK22ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK22,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  SSPRK22Cache(u,uprev,k,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK22,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) = SSPRK22ConstantCache()


@cache struct SSPRK33Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType # superfluous, only needed for callbacks...
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

struct SSPRK33ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK33,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  SSPRK33Cache(u,uprev,k,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK33,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) = SSPRK33ConstantCache()

@cache struct KYKSSPRK42Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
   u::uType
   uprev::uType
   k::rateType
   tmp::uType
   fsalfirst::rateType
   tab::TabType
 end

struct KYKSSPRK42ConstantCache{T, T2} <: OrdinaryDiffEqConstantCache
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
  α30 = T(0.002797307087390)
  α32 = T(0.997202692912610)
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
  KYKSSPRK42ConstantCache(α20, α21, α30, α32, α40, α43, β10, β21, β30, β32, β40, β43, c1, c2, c3)
end

function alg_cache(alg::KYKSSPRK42,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = KYKSSPRK42ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
  KYKSSPRK42Cache(u,uprev,k,tmp, fsalfirst,tab)
end

function alg_cache(alg::KYKSSPRK42,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  KYKSSPRK42ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct SSPRK53Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK53ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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
    β43 = T(0.238458932846290)
    β54 = T(0.287632146308408)
    c1 = T2(0.377268915331368)
    c2 = T2(0.754537830662736)
    c3 = T2(0.728985661612188)
    c4 = T2(0.699226135931670)

    new{T,T2}(α30, α32, α40, α43, α52, α54, β10, β21, β32, β43, β54, c1, c2, c3, c4)
  end
end

function alg_cache(alg::SSPRK53,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = SSPRK53ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK53Cache(u,uprev,k,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK53,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  SSPRK53ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct SHLDDRK52Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct SHLDDRK52ConstantCache{T1,T2} <: OrdinaryDiffEqConstantCache
  α2::T1
  α3::T1
  α4::T1
  α5::T1
  β1::T1
  β2::T1
  β3::T1
  β4::T1
  β5::T1
  c2::T2
  c3::T2
  c4::T2
  c5::T2
end

function SHLDDRK52ConstantCache(T1,T2)
  α2 = T1(-0.6913065)
  α3 = T1(-2.655155)
  α4 = T1(-0.8147688)
  α5 = T1(-0.6686587)
  β1 = T1(0.1)
  β2 = T1(0.75)
  β3 = T1(0.7)
  β4 = T1(0.479313)
  β5 = T1(0.310392)
  c2 = T2(0.1)
  c3 = T2(0.3315201)
  c4 = T2(0.4577796)
  c5 = T2(0.8666528)
  SHLDDRK52ConstantCache(α2, α3, α4, α5, β1, β2, β3, β4, β5, c2, c3, c4, c5)
end

function alg_cache(alg::SHLDDRK52,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  SHLDDRK52ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end

function alg_cache(alg::SHLDDRK52,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = SHLDDRK52ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SHLDDRK52Cache(u,uprev,k,tmp,fsalfirst,tab)
end

@cache mutable struct SHLDDRK_2NCache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
  step::Int
end

mutable struct SHLDDRK_2NConstantCache{T1,T2} <: OrdinaryDiffEqConstantCache
  α21::T1
  α31::T1
  α41::T1
  α51::T1
  β11::T1
  β21::T1
  β31::T1
  β41::T1
  β51::T1
  c21::T2
  c31::T2
  c41::T2
  c51::T2

  α22::T1
  α32::T1
  α42::T1
  α52::T1
  α62::T1
  β12::T1
  β22::T1
  β32::T1
  β42::T1
  β52::T1
  β62::T1
  c22::T2
  c32::T2
  c42::T2
  c52::T2
  c62::T2

  step::Int
end

function SHLDDRK_2NConstantCache(T1,T2)
  α21 = T1(-0.6051226)
  α31 = T1(-2.0437564)
  α41 = T1(-0.7406999)
  α51 = T1(-4.4231765)
  β11 = T1(0.2687454)
  β21 = T1(0.8014706)
  β31 = T1(0.5051570)
  β41 = T1(0.5623568)
  β51 = T1(0.0590065)
  c21 = T2(0.2687454)
  c31 = T2(0.5852280)
  c41 = T2(0.6827066)
  c51 = T2(1.1646854)

  α22 = T1(-0.4412737)
  α32 = T1(-1.0739820)
  α42 = T1(-1.7063570)
  α52 = T1(-2.7979293)
  α62 = T1(-4.0913537)
  β12 = T1(0.1158488)
  β22 = T1(0.3728769)
  β32 = T1(0.7379536)
  β42 = T1(0.5798110)
  β52 = T1(1.0312849)
  β62 = T1(0.15)
  c22 = T2(0.1158485)
  c32 = T2(0.3241850)
  c42 = T2(0.6193208)
  c52 = T2(0.8034472)
  c62 = T2(0.9184166)
  SHLDDRK_2NConstantCache(α21, α31, α41, α51, β11, β21, β31, β41, β51, c21, c31, c41, c51, α22, α32, α42, α52, α62, β12, β22, β32, β42, β52, β62, c22, c32, c42, c52, c62, 1)
end

function alg_cache(alg::SHLDDRK_2N,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  SHLDDRK_2NConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
end

function alg_cache(alg::SHLDDRK_2N,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = SHLDDRK_2NConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SHLDDRK_2NCache(u,uprev,k,tmp,fsalfirst,tab,1)
end

@cache struct SSPRK53_2N1Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK53_2N1ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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

    new{T,T2}(α40, α43,β10, β21, β32, β43, β54, c1, c2, c3, c4)
  end
end

function alg_cache(alg::SSPRK53_2N1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = SSPRK53_2N1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK53_2N1Cache(u,uprev,k,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK53_2N1,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  SSPRK53_2N1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


@cache struct SSPRK53_2N2Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK53_2N2ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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
    α50 = T(0.045230974482400)
    α54 = T(0.954769025517600)
    β10 = T(0.465388589249323)
    β21 = T(0.465388589249323)
    β32 = T(0.124745797313998)
    β43 = T(0.465388589249323)
    β54 = T(0.154263303748666)
    c1 = T2(0.465388589249323)
    c2 = T2(0.930777178498646)
    c3 = T2(0.420413812847710)
    c4 = T2(0.885802402097033)

    new{T,T2}(α30, α32,α50,α54,β10, β21, β32, β43, β54, c1, c2, c3, c4)
  end
end

function alg_cache(alg::SSPRK53_2N2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = SSPRK53_2N2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK53_2N2Cache(u,uprev,k,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK53_2N2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  SSPRK53_2N2ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct SSPRK53_HCache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK53_HConstantCache{T,T2} <: OrdinaryDiffEqConstantCache

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

  function SSPRK53_HConstantCache(T,T2)


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

    new{T,T2}(α30, α32, α40, α41, α43, β10, β21, β32, β43, β54, c1, c2, c3, c4)
  end
end

function alg_cache(alg::SSPRK53_H,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  tab = SSPRK53_HConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK53_HCache(u,uprev,k,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK53_H,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  SSPRK53_HConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end



@cache struct SSPRK63Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  u₂::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK63ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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
    α43 = T(0.424718455428740)
    α62 = T(0.155221702560091)
    α65 = T(0.844778297439909)
    β10 = T(0.284220721334261)
    β21 = T(0.284220721334261)
    β32 = T(0.284220721334261)
    β43 = T(0.120713785765930)
    β54 = T(0.284220721334261)
    β65 = T(0.240103497065900)
    c1 = T2(0.284220721334261)
    c2 = T2(0.568441442668522)
    c3 = T2(0.852662164002783)
    c4 = T2(0.510854218958172)
    c5 = T2(0.795074940292433)

    new{T,T2}(α40, α41, α43, α62, α65, β10, β21, β32, β43, β54, β65, c1, c2, c3, c4, c5)
  end
end

function alg_cache(alg::SSPRK63,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  u₂ = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = SSPRK63ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK63Cache(u,uprev,k,tmp,u₂,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK63,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  SSPRK63ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


@cache struct SSPRK73Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  u₁::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK73ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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
    α50 = T(0.180718656570380)
    α51 = T(0.314831034403793)
    α54 = T(0.504450309025826)
    α73 = T(0.120199000000000)
    α76 = T(0.879801000000000)
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

    new{T,T2}(α40, α43, α50, α51, α54, α73, α76, β10, β21, β32, β43, β54, β65, β76, c1, c2, c3, c4, c5, c6)
  end
end

function alg_cache(alg::SSPRK73,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  u₁ = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = SSPRK73ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK73Cache(u,uprev,k,tmp,u₁,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK73,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  SSPRK73ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


@cache struct SSPRK83Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  u₂::uType
  u₃::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK83ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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
    α76 = T(0.652354615999950)
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

    new{T,T2}(α50, α51, α54, α61, α65, α72, α73, α76, β10, β21, β32, β43, β54, β65, β76, β87, c1, c2, c3, c4, c5, c6, c7)
  end
end

function alg_cache(alg::SSPRK83,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  u₂ = similar(u)
  u₃ = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = SSPRK83ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK83Cache(u,uprev,k,tmp,u₂,u₃,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK83,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  SSPRK83ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


@cache struct SSPRK432Cache{uType,rateType,uNoUnitsType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType # superfluous, only needed for callbacks
  fsalfirst::rateType
  utilde::uType
  atmp::uNoUnitsType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

struct SSPRK432ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK432,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  SSPRK432Cache(u,uprev,k,tmp,fsalfirst,utilde,atmp,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK432,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) = SSPRK432ConstantCache()


@cache mutable struct SSPRKMSVS32Cache{uType,rateType,dtArrayType,dtType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
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
  step::Int
end

@cache mutable struct SSPRKMSVS32ConstantCache{uType,dtArrayType,dtType} <: OrdinaryDiffEqConstantCache
  u_2::uType
  u_1::uType
  dts::dtArrayType
  dtf::dtArrayType
  μ::dtType
  v_n::Float64
  step::Int
end

function alg_cache(alg::SSPRKMSVS32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  fsalfirst = zero(rate_prototype)
  dts = fill(zero(dt),3)
  dtf = fill(zero(dt),2)
  μ = zero(dt)
  u_2 = similar(u)
  u_1 = similar(u)
  k  = zero(rate_prototype)
  tmp = similar(u)
  SSPRKMSVS32Cache(u,uprev,fsalfirst,u_2,u_1,k,tmp,dts,dtf,μ,0.5,alg.stage_limiter!,alg.step_limiter!,1)
end

function alg_cache(alg::SSPRKMSVS32,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  dts = fill(zero(dt),3)
  dtf = fill(zero(dt),2)
  μ = zero(dt)
  u_2 = u
  u_1 = u
  SSPRKMSVS32ConstantCache(u_2,u_1,dts,dtf,μ,0.5,1)
end


@cache mutable struct SSPRKMSVS43Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
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
  step::Int
end

@cache mutable struct SSPRKMSVS43ConstantCache{uType,rateType} <: OrdinaryDiffEqConstantCache
  u_3::uType
  u_2::uType
  u_1::uType
  k1::rateType
  k2::rateType
  k3::rateType
  step::Int
end

function alg_cache(alg::SSPRKMSVS43,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  fsalfirst = zero(rate_prototype)
  u_3 = similar(u)
  u_2 = similar(u)
  u_1 = similar(u)
  k  = zero(rate_prototype)
  k1  = zero(rate_prototype)
  k2  = zero(rate_prototype)
  k3  = zero(rate_prototype)
  tmp = similar(u)
  SSPRKMSVS43Cache(u,uprev,fsalfirst,u_3,u_2,u_1,k,k1,k2,k3,tmp,alg.stage_limiter!,alg.step_limiter!,1)
end

function alg_cache(alg::SSPRKMSVS43,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  u_3 = u
  u_2 = u
  u_1 = u
  k1  = rate_prototype
  k2  = rate_prototype
  k3  = rate_prototype
  SSPRKMSVS43ConstantCache(u_3,u_2,u_1,k1,k2,k3,1)
end


@cache struct SSPRK932Cache{uType,rateType,uNoUnitsType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType # superfluous, only needed for callbacks
  fsalfirst::rateType
  utilde::uType
  atmp::uNoUnitsType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

struct SSPRK932ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK932,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  SSPRK932Cache(u,uprev,k,tmp,fsalfirst,utilde,atmp,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK932,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) = SSPRK932ConstantCache()


@cache struct SSPRK54Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₃::rateType
  u₂::uType
  u₃::uType
  tmp::uType # should be u₄, but tmp is needed for callbacks
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK54ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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
    β10 = T(0.391752226571890)
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
    β53 = T(0.063692468666290)
    α54 = T(0.386708617503269)
    β54 = T(0.226007483236906)
    c1 = T2(0.391752226571890)
    c2 = T2(0.586079689311540)
    c3 = T2(0.474542363121400)
    c4 = T2(0.935010630967653)

    new{T,T2}(β10, α20, α21, β21, α30, α32, β32, α40, α43, β43, α52, α53, β53, α54, β54, c1, c2, c3, c4)
  end
end

function alg_cache(alg::SSPRK54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  u₂ = similar(u)
  u₃ = similar(u)
  tmp = similar(u)
  k = zero(rate_prototype)
  k₃ = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = SSPRK54ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK54Cache(u,uprev,k,k₃,u₂,u₃,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  SSPRK54ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache mutable struct VTSRKCache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  k::rateType
  k₃::rateType
  u₂::uType
  u₃::uType
  y1::uType
  y2::uType
  y3::uType
  y4::uType
  y5::uType
  y6::uType
  y7::uType
  y8::uType
  tmp::uType # should be u₄, but tmp is needed for callbacks
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
  step::Int
end

mutable struct VTSRKConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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
  d7::T
  n2::T
  n3::T
  n6::T
  n8::T
  q20::T
  q21::T
  q30::T
  q32::T
  q41::T
  q43::T
  q51::T
  q54::T
  q61::T
  q65::T
  q70::T
  q71::T
  q76::T
  q80::T
  q81::T
  q82::T
  q87::T
  r::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2
  c8::T2
  c1hat::T2
  c2hat::T2
  c3hat::T2
  c4hat::T2
  step::Int
end
function VTSRKConstantCache(T, T2)
  β10 = convert(T,0.391752226571890)
  α20 = convert(T,0.444370493651235)
  α21 = convert(T,0.555629506348765)
  β21 = convert(T,0.368410593050371)
  α30 = convert(T,0.620101851488403)
  α32 = convert(T,0.379898148511597)
  β32 = convert(T,0.251891774271694)
  α40 = convert(T,0.178079954393132)
  α43 = convert(T,0.821920045606868)
  β43 = convert(T,0.544974750228521)
  α52 = convert(T,0.517231671970585)
  α53 = convert(T,0.096059710526147)
  β53 = convert(T,0.063692468666290)
  α54 = convert(T,0.386708617503269)
  β54 = convert(T,0.226007483236906)
  d7 = convert(T,0.003674184820260)
  n2 = convert(T,0.179502832154858)
  n3 = convert(T,0.073789956884809)
  n6 = convert(T,0.017607159013167)
  n8 = convert(T,0.729100051947166)
  q20 = convert(T,0.085330772947643)
  q21 = convert(T,0.914669227052357)
  q30 = convert(T,0.058121281984411)
  q32 = convert(T,0.941878718015589)
  q41 = convert(T,0.036365639242841)
  q43 = convert(T,0.802870131352638)
  q51 = convert(T,0.491214340660555)
  q54 = convert(T,0.508785659339445)
  q61 = convert(T,0.566135231631241)
  q65 = convert(T,0.433864768368758)
  q70 = convert(T,0.020705281786630)
  q71 = convert(T,0.091646079651566)
  q76 = convert(T,0.883974453741544)
  q80 = convert(T,0.008506650138784)
  q81 = convert(T,0.110261531523242)
  q82 = convert(T,0.030113037742445)
  q87 = convert(T,0.851118780595529)
  r = convert(T,3.5794)
  c2 = convert(T2,0.19404565885657)
  c3 = convert(T2,0.40402262622011853)
  c4 = convert(T2,0.5588403940142084)
  c5 = convert(T2,0.5637064101382472)
  c6 = convert(T2,0.5239487828668272)
  c7 = convert(T2,0.7171278236757)
  c8 = convert(T2,0.887074044732322)
  c1hat = convert(T2,0.391752226571890)
  c2hat = convert(T2,0.586079689311540)
  c3hat = convert(T2,0.474542363121400)
  c4hat = convert(T2,0.935010630967653)

  VTSRKConstantCache(β10, α20, α21, β21, α30, α32, β32, α40, α43, β43, α52, α53, β53, α54, β54, d7, n2, n3, n6, n8, q20, q21, q30, q32, q41, q43, q51, q54, q61, q65, q70, q71, q76, q80, q81, q82, q87,r, c2,c3, c4, c5, c6, c7, c8, c1hat, c2hat, c3hat, c4hat,1)
end


function alg_cache(alg::VTSRK,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  y1 = similar(u)
  y2 = similar(u)
  y3 = similar(u)
  y4 = similar(u)
  y5 = similar(u)
  y6 = similar(u)
  y7 = similar(u)
  y8 = similar(u)
  u₂ = similar(u)
  u₃ = similar(u)
  tmp = similar(u)
  k = zero(rate_prototype)
  k₃ = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = VTSRKConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  VTSRKCache(u,uprev,uprev2,k,k₃,u₂,u₃,y1,y2,y3,y4,y5,y6,y7,y8,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab,1)
end

function alg_cache(alg::VTSRK,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  VTSRKConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct SSPRK104Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₄::rateType
  tmp::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

struct SSPRK104ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK104,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  tmp = similar(u)
  k = zero(rate_prototype)
  k₄ = zero(rate_prototype)
  if calck
    fsalfirst = zero(rate_prototype)
  else
    fsalfirst = k
  end
  SSPRK104Cache(u,uprev,k,k₄,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK104,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false}) = SSPRK104ConstantCache()
