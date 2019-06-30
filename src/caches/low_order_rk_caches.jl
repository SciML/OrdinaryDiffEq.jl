@cache struct EulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
end

@cache struct SplitEulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
end

function alg_cache(alg::SplitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  SplitEulerCache(u,uprev,similar(u),zero(rate_prototype),zero(rate_prototype))
end

struct SplitEulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::SplitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = SplitEulerConstantCache()

function alg_cache(alg::Euler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  EulerCache(u,uprev,similar(u),zero(rate_prototype),zero(rate_prototype))
end

struct EulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Euler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = EulerConstantCache()

@cache struct HeunCache{uType,rateType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  atmp::uNoUnitsType
  k::rateType
  fsalfirst::rateType
end

@cache struct RalstonCache{uType,rateType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  atmp::uNoUnitsType
  k::rateType
  fsalfirst::rateType
end

function alg_cache(alg::Heun,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  HeunCache(u,uprev,similar(u),similar(u, uEltypeNoUnits),zero(rate_prototype),
            zero(rate_prototype))
end

function alg_cache(alg::Ralston,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  RalstonCache(u,uprev,similar(u),similar(u, uEltypeNoUnits),zero(rate_prototype),
               zero(rate_prototype))
end

struct HeunConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Heun,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = HeunConstantCache()

struct RalstonConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Ralston,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = RalstonConstantCache()

@cache struct MidpointCache{uType,rateType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  atmp::uNoUnitsType
  fsalfirst::rateType
end

struct MidpointConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::Midpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u); atmp = similar(u, uEltypeNoUnits)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  MidpointCache(u,uprev,k,tmp,atmp,fsalfirst)
end

alg_cache(alg::Midpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = MidpointConstantCache()

@cache struct RK4Cache{uType,rateType,uNoUnitsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k₂::rateType
  k₃::rateType
  k₄::rateType
  k::rateType
  tmp::uType
  atmp::uNoUnitsType
end

struct RK4ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::RK4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  k₁ = zero(rate_prototype)
  k₂ = zero(rate_prototype)
  k₃ = zero(rate_prototype)
  k₄ = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = similar(u); atmp = similar(u, uEltypeNoUnits)
  RK4Cache(u,uprev,k₁,k₂,k₃,k₄,k,tmp,atmp)
end

alg_cache(alg::RK4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = RK4ConstantCache()

@cache struct BS3Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
end

function alg_cache(alg::BS3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = BS3ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  BS3Cache(u,uprev,k1,k2,k3,k4,utilde,tmp,atmp,tab)
end

alg_cache(alg::BS3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = BS3ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))


@cache struct OwrenZen3Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
end

function alg_cache(alg::OwrenZen3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = OwrenZen3ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  OwrenZen3Cache(u,uprev,k1,k2,k3,k4,utilde,tmp,atmp,tab)
end

alg_cache(alg::OwrenZen3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = OwrenZen3ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

@cache struct OwrenZen4Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
end

function alg_cache(alg::OwrenZen4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = OwrenZen4ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k6 = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  OwrenZen4Cache(u,uprev,k1,k2,k3,k4,k5,k6,utilde,tmp,atmp,tab)
end

alg_cache(alg::OwrenZen4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = OwrenZen4ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

@cache struct OwrenZen5Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
end

function alg_cache(alg::OwrenZen5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = OwrenZen5ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k6 = zero(rate_prototype)
  k7 = zero(rate_prototype)
  k8 = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  OwrenZen5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,utilde,tmp,atmp,tab)
end

alg_cache(alg::OwrenZen5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = OwrenZen5ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

@cache struct BS5Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
end

function alg_cache(alg::BS5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = BS5ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k6 = zero(rate_prototype)
  k7 = zero(rate_prototype)
  k8 = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  BS5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,utilde,tmp,atmp,tab)
end

alg_cache(alg::BS5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = BS5ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

@cache struct Tsit5Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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

@cache struct RK46NLCache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct RK46NLConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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

  function RK46NLConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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
    new{T,T2}(α2, α3, α4, α5, α6, β1, β2, β3, β4, β5, β6, c2, c3, c4, c5, c6)
  end
end

function alg_cache(alg::RK46NL,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = RK46NLConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  RK46NLCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::RK46NL,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  RK46NLConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

function alg_cache(alg::Tsit5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = Tsit5ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k6 = zero(rate_prototype)
  k7 = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  Tsit5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp,tab)
end

alg_cache(alg::Tsit5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = Tsit5ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

@cache struct DP5Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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
  tab::TabType
end

function alg_cache(alg::DP5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
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

  tmp = similar(u) # has to be separate for FSAL
  utilde = tmp

  if eltype(u) != uEltypeNoUnits || calck
    update = zero(rate_prototype)
    atmp = similar(u,uEltypeNoUnits)
  else
    update = k7
    atmp = k3
  end

  tab = DP5ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  cache = DP5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp,tab)
  cache
end

alg_cache(alg::DP5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = DP5ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

@cache struct Anas5Cache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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

function alg_cache(alg::Anas5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tab = Anas5ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k6 = zero(rate_prototype)
  k7 = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  Anas5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp,tab)
end

alg_cache(alg::Anas5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = Anas5ConstantCache(constvalue(uBottomEltypeNoUnits),constvalue(tTypeNoUnits))

@cache struct KYK2014DGSSPRK_3S2_Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
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
end

struct KYK2014DGSSPRK_3S2_ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
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

  function KYK2014DGSSPRK_3S2_ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
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
    new{T,T2}(α_10, α_20, α_21, α_30, α_32, β_10, β_21, β_30, β_32, c_1, c_2)
  end
end

function alg_cache(alg::KYK2014DGSSPRK_3S2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  u_1 = similar(u)
  u_2 = similar(u)
  kk_1 = zero(rate_prototype)
  kk_2 = zero(rate_prototype)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = KYK2014DGSSPRK_3S2_ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  KYK2014DGSSPRK_3S2_Cache(u,uprev,k,fsalfirst,tab, u_1, u_2, kk_1, kk_2)
end

function alg_cache(alg::KYK2014DGSSPRK_3S2,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  KYK2014DGSSPRK_3S2_ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end
