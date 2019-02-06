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
  tab = BS3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  BS3Cache(u,uprev,k1,k2,k3,k4,utilde,tmp,atmp,tab)
end

alg_cache(alg::BS3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = BS3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))


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
  tab = OwrenZen3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  OwrenZen3Cache(u,uprev,k1,k2,k3,k4,utilde,tmp,atmp,tab)
end

alg_cache(alg::OwrenZen3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = OwrenZen3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

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
  tab = OwrenZen4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
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

alg_cache(alg::OwrenZen4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = OwrenZen4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

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
  tab = OwrenZen5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
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

alg_cache(alg::OwrenZen5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = OwrenZen5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

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
  tab = BS5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
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

alg_cache(alg::BS5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = BS5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

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
  tab = Tsit5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
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

alg_cache(alg::Tsit5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = Tsit5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

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

  tab = DP5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  cache = DP5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp,tab)
  cache
end

alg_cache(alg::DP5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = DP5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

@cache struct DP5ThreadedCache{uType,rateType,uNoUnitsType,TabType} <: OrdinaryDiffEqMutableCache
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

function alg_cache(alg::DP5Threaded,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  k1 = zero(rate_prototype)
  k2 = zero(rate_prototype)
  k3 = zero(rate_prototype)
  k4 = zero(rate_prototype)
  k5 = zero(rate_prototype)
  k6 = zero(rate_prototype)
  k7 = zero(rate_prototype)
  dense_tmp3 = zero(rate_prototype)
  dense_tmp4 = zero(rate_prototype)
  update = zero(rate_prototype)
  bspl = zero(rate_prototype)
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  tab = DP5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  DP5ThreadedCache(u,uprev,k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp,tab)
end

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
  tab = Anas5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
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

alg_cache(alg::Anas5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = Anas5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

@cache struct CFRLDDRK64Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct CFRLDDRK64ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  α1::T
  α2::T
  α3::T
  α4::T
  α5::T
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

  function CFRLDDRK64ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    α1 = T(0.17985400977138)
    α2 = T(0.14081893152111)
    α3 = T(0.08255631629428)
    α4 = T(0.65804425034331)
    α5 = T(0.31862993413251)
    β1 = T(0.10893125722541)
    β2 = T(0.13201701492152)
    β3 = T(0.38911623225517)
    β4 = T(-0.59203884581148)
    β5 = T(0.47385028714844)
    β6 = T(0.48812405426094)
    c2 = T2(0.28878526699679)
    c3 = T2(0.38176720366804)
    c4 = T2(0.71262082069639)
    c5 = T2(0.69606990893393)
    c6 = T2(0.83050587987157)
    new{T,T2}(α1, α2, α3, α4, α5, β1, β2, β3, β4, β5, β6, c2, c3, c4, c5, c6)
  end
end

function alg_cache(alg::CFRLDDRK64,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = CFRLDDRK64ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  CFRLDDRK64Cache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::CFRLDDRK64,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  CFRLDDRK64ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct TSLDDRK74Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct TSLDDRK74ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  α1::T
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
  β7::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2

  function TSLDDRK74ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    α1 = T(0.241566650129646868)
    α2 = T(0.0423866513027719953)
    α3 = T(0.215602732678803776)
    α4 = T(0.232328007537583987)
    α5 = T(0.256223412574146438)
    α6 = T(0.0978694102142697230)
    β1 = T(0.0941840925477795334)
    β2 = T(0.149683694803496998)
    β3 = T(0.285204742060440058)
    β4 = T(-0.122201846148053668)
    β5 = T(0.0605151571191401122)
    β6 = T(0.345986987898399296)
    β7 = T(0.186627171718797670)
    c2 = T2(0.335750742677426401)
    c3 = T2(0.286254438654048527)
    c4 = T2(0.744675262090520366)
    c5 = T2(0.639198690801246909)
    c6 = T2(0.723609252956949472)
    c7 = T2(0.91124223849547205)
    new{T,T2}(α1, α2, α3, α4, α5, α6, β1, β2, β3, β4, β5, β6, β7, c2, c3, c4, c5, c6, c7)
  end
end

function alg_cache(alg::TSLDDRK74,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = TSLDDRK74ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  TSLDDRK74Cache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::TSLDDRK74,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  TSLDDRK74ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


@cache struct NDBLSRK124Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct NDBLSRK124ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  α2::T
  α3::T
  α4::T
  α5::T
  α6::T
  α7::T
  α8::T
  α9::T
  α10::T
  α11::T
  α12::T
  β1::T
  β2::T
  β3::T
  β4::T
  β5::T
  β6::T
  β7::T
  β8::T
  β9::T
  β10::T
  β11::T
  β12::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2
  c8::T2
  c9::T2
  c10::T2
  c11::T2
  c12::T2

  function NDBLSRK124ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    α2  = T(-0.0923311242368072)
    α3  = T(-0.9441056581158819)
    α4  = T(-4.3271273247576394)
    α5  = T(-2.1557771329026072)
    α6  = T(-0.9770727190189062)
    α7  = T(-0.7581835342571139)
    α8  = T(-1.7977525470825499)
    α9  = T(-2.6915667972700770)
    α10 = T(-4.6466798960268143)
    α11 = T(-0.1539613783825189)
    α12 = T(-0.5943293901830616)
    β1  = T(0.0650008435125904)
    β2  = T(0.0161459902249842)
    β3  = T(0.5758627178358159)
    β4  = T(0.1649758848361671)
    β5  = T(0.3934619494248182)
    β6  = T(0.0443509641602719)
    β7  = T(0.2074504268408778)
    β8  = T(0.6914247433015102)
    β9  = T(0.3766646883450449)
    β10 = T(0.0757190350155483)
    β11 = T(0.2027862031054088)
    β12 = T(0.2167029365631842)
    c2  = T2(0.0650008435125904)
    c3  = T2(0.0796560563081853)
    c4  = T2(0.1620416710085376)
    c5  = T2(0.2248877362907778)
    c6  = T2(0.2952293985641261)
    c7  = T2(0.3318332506149405)
    c8  = T2(0.4094724050198658)
    c9  = T2(0.6356954475753369)
    c10 = T2(0.6806551557645497)
    c11 = T2(0.7143773712418350)
    c12 = T2(0.9032588871651854)
    new{T,T2}(α2, α3, α4, α5, α6, α7, α8, α9, α10, α11, α12, β1, β2, β3, β4, β5, β6, β7, β8, β9, β10, β11, β12, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12)
  end
end

function alg_cache(alg::NDBLSRK124,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = NDBLSRK124ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  NDBLSRK124Cache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::NDBLSRK124,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  NDBLSRK124ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct NDBLSRK134Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
    u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct NDBLSRK134ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  α2::T
  α3::T
  α4::T
  α5::T
  α6::T
  α7::T
  α8::T
  α9::T
  α10::T
  α11::T
  α12::T
  α13::T
  β1::T
  β2::T
  β3::T
  β4::T
  β5::T
  β6::T
  β7::T
  β8::T
  β9::T
  β10::T
  β11::T
  β12::T
  β13::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2
  c8::T2
  c9::T2
  c10::T2
  c11::T2
  c12::T2
  c13::T2

  function NDBLSRK134ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    α2  = T(-0.6160178650170565)
    α3  = T(-0.4449487060774118)
    α4  = T(-1.0952033345276178)
    α5  = T(-1.2256030785959187)
    α6  = T(-0.2740182222332805)
    α7  = T(-0.0411952089052647)
    α8  = T(-0.1797084899153560)
    α9  = T(-1.1771530652064288)
    α10 = T(-0.4078831463120878)
    α11 = T(-0.8295636426191777)
    α12 = T(-4.7895970584252288)
    α13 = T(-0.6606671432964504)
    β1  = T(0.0271990297818803)
    β2  = T(0.1772488819905108)
    β3  = T(0.0378528418949694)
    β4  = T(0.6086431830142991)
    β5  = T(0.2154313974316100)
    β6  = T(0.2066152563885843)
    β7  = T(0.0415864076069797)
    β8  = T(0.0219891884310925)
    β9  = T(0.9893081222650993)
    β10 = T(0.0063199019859826)
    β11 = T(0.3749640721105318)
    β12 = T(1.6080235151003195)
    β13 = T(0.0961209123818189)
    c2  = T2(0.0271990297818803)
    c3  = T2(0.0952594339119365)
    c4  = T2(0.1266450286591127)
    c5  = T2(0.1825883045699772)
    c6  = T2(0.3737511439063931)
    c7  = T2(0.5301279418422206)
    c8  = T2(0.5704177433952291)
    c9  = T2(0.5885784947099155)
    c10 = T2(0.6160769826246714)
    c11 = T2(0.6223252334314046)
    c12 = T2(0.6897593128753419)
    c13 = T2(0.9126827615920843)
    new{T,T2}(α2, α3, α4, α5, α6, α7, α8, α9, α10, α11, α12, α13, β1, β2, β3, β4, β5, β6, β7, β8, β9, β10, β11, β12, β13, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13)
  end
end

function alg_cache(alg::NDBLSRK134,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = NDBLSRK134ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  NDBLSRK134Cache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::NDBLSRK134,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  NDBLSRK134ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct NDBLSRK144Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct NDBLSRK144ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  α2::T
  α3::T
  α4::T
  α5::T
  α6::T
  α7::T
  α8::T
  α9::T
  α10::T
  α11::T
  α12::T
  α13::T
  α14::T
  β1::T
  β2::T
  β3::T
  β4::T
  β5::T
  β6::T
  β7::T
  β8::T
  β9::T
  β10::T
  β11::T
  β12::T
  β13::T
  β14::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2
  c8::T2
  c9::T2
  c10::T2
  c11::T2
  c12::T2
  c13::T2
  c14::T2

  function NDBLSRK144ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    α2  = T(-0.7188012108672410)
    α3  = T(-0.7785331173421570)
    α4  = T(-0.0053282796654044)
    α5  = T(-0.8552979934029281)
    α6  = T(-3.9564138245774565)
    α7  = T(-1.5780575380587385)
    α8  = T(-2.0837094552574054)
    α9  = T(-0.7483334182761610)
    α10 = T(-0.7032861106563359)
    α11 = T(0.0013917096117681)
    α12 = T(-0.0932075369637460)
    α13 = T(-0.9514200470875948)
    α14 = T(-7.1151571693922548)
    β1  = T(0.0367762454319673)
    β2  = T(0.3136296607553959)
    β3  = T(0.1531848691869027)
    β4  = T(0.0030097086818182)
    β5  = T(0.3326293790646110)
    β6  = T(0.2440251405350864)
    β7  = T(0.3718879239592277)
    β8  = T(0.6204126221582444)
    β9  = T(0.1524043173028741)
    β10 = T(0.0760894927419266)
    β11 = T(0.0077604214040978)
    β12 = T(0.0024647284755382)
    β13 = T(0.0780348340049386)
    β14 = T(5.5059777270269628)
    c2  = T2(0.0367762454319673)
    c3  = T2(0.1249685262725025)
    c4  = T2(0.2446177702277698)
    c5  = T2(0.2476149531070420)
    c6  = T2(0.2969311120382472)
    c7  = T2(0.3978149645802642)
    c8  = T2(0.5270854589440328)
    c9  = T2(0.6981269994175695)
    c10 = T2(0.8190890835352128)
    c11 = T2(0.8527059887098624)
    c12 = T2(0.8604711817462826)
    c13 = T2(0.8627060376969976)
    c14 = T2(0.8734213127600976)
    new{T,T2}(α2, α3, α4, α5, α6, α7, α8, α9, α10, α11, α12, α13, α14, β1, β2, β3, β4, β5, β6, β7, β8, β9, β10, β11, β12, β13, β14, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14)
  end
end

function alg_cache(alg::NDBLSRK144,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = NDBLSRK144ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  NDBLSRK144Cache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::NDBLSRK144,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  NDBLSRK144ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct DGLDDRK73_CCache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct DGLDDRK73_CConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  α2::T
  α3::T
  α4::T
  α5::T
  α6::T
  α7::T
  β1::T
  β2::T
  β3::T
  β4::T
  β5::T
  β6::T
  β7::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2

  function DGLDDRK73_CConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    α2 = T(-0.8083163874983830)
    α3 = T(-1.503407858773331)
    α4 = T(-1.053064525050744)
    α5 = T(-1.463149119280508)
    α6 = T(-0.6592881281087830)
    α7 = T(-1.667891931891068)
    β1 = T(0.01197052673097840)
    β2 = T(0.8886897793820711)
    β3 = T(0.4578382089261419)
    β4 = T(0.5790045253338471)
    β5 = T(0.3160214638138484)
    β6 = T(0.2483525368264122)
    β7 = T(0.06771230959408840)
    c2 = T2(0.01197052673097840)
    c3 = T2(0.1823177940361990)
    c4 = T2(0.5082168062551849)
    c5 = T2(0.6532031220148590)
    c6 = T2(0.8534401385678250)
    c7 = T2(0.9980466084623790)
    new{T,T2}(α2, α3, α4, α5, α6, α7, β1, β2, β3, β4, β5, β6, β7, c2, c3, c4, c5, c6, c7)
  end
end

function alg_cache(alg::DGLDDRK73_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = DGLDDRK73_CConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  DGLDDRK73_CCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::DGLDDRK73_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  DGLDDRK73_CConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


@cache struct DGLDDRK84_CCache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

struct DGLDDRK84_CConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  α2::T
  α3::T
  α4::T
  α5::T
  α6::T
  α7::T
  α8::T
  β1::T
  β2::T
  β3::T
  β4::T
  β5::T
  β6::T
  β7::T
  β8::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2
  c8::T2

  function DGLDDRK84_CConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    α2  = T(-0.7212962482279240)
    α3  = T(-0.01077336571612980)
    α4  = T(-0.5162584698930970)
    α5  = T(-1.730100286632201)
    α6  = T(-5.200129304403076)
    α7  = T(0.7837058945416420)
    α8  = T(-0.5445836094332190)
    β1  = T(0.2165936736758085)
    β2  = T(0.1773950826411583)
    β3  = T(0.01802538611623290)
    β4  = T(0.08473476372541490)
    β5  = T(0.8129106974622483)
    β6  = T(1.903416030422760)
    β7  = T(0.1314841743399048)
    β8  = T(0.2082583170674149)
    c2  = T2(0.2165936736758085)
    c3  = T2(0.2660343487538170)
    c4  = T2(0.2840056122522720)
    c5  = T2(0.3251266843788570)
    c6  = T2(0.4555149599187530)
    c7  = T2(0.7713219317101170)
    c8  = T2(0.9199028964538660)

    new{T,T2}(α2, α3, α4, α5, α6, α7, α8, β1, β2, β3, β4, β5, β6, β7, β8, c2, c3, c4, c5, c6, c7, c8)
  end
end

function alg_cache(alg::DGLDDRK84_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = DGLDDRK84_CConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  DGLDDRK84_CCache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::DGLDDRK84_C,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  DGLDDRK84_CConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end
