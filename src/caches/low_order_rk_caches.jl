struct EulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
end

struct SplitEulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
end

function alg_cache(alg::SplitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  SplitEulerCache(u,uprev,similar(u),zero(rate_prototype),zero(rate_prototype))
end

u_cache(c::SplitEulerCache) = ()
du_cache(c::SplitEulerCache) = (c.k,c.fsalfirst)

struct SplitEulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::SplitEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = SplitEulerConstantCache()

function alg_cache(alg::Euler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  EulerCache(u,uprev,similar(u),zero(rate_prototype),zero(rate_prototype))
end

u_cache(c::EulerCache) = ()
du_cache(c::EulerCache) = (c.k,c.fsalfirst)

struct EulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Euler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = EulerConstantCache()

struct HeunCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  utilde::rateType
  fsalfirst::rateType
end

struct RalstonCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  utilde::rateType
  fsalfirst::rateType
end

u_cache(c::Union{HeunCache,RalstonCache}) = ()
du_cache(c::Union{HeunCache,RalstonCache}) = (c.k,c.fsalfirst,c.utilde)

function alg_cache(alg::Heun,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  HeunCache(u,uprev,similar(u),zero(rate_prototype),zero(rate_prototype),zero(rate_prototype))
end

function alg_cache(alg::Ralston,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  RalstonCache(u,uprev,similar(u),zero(rate_prototype),zero(rate_prototype),zero(rate_prototype))
end

struct HeunConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Heun,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = HeunConstantCache()

struct RalstonConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Ralston,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = RalstonConstantCache()

struct MidpointCache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  atmp::uEltypeNoUnits
  fsalfirst::rateType
end

u_cache(c::MidpointCache) = (c.atmp,)
du_cache(c::MidpointCache) = (c.k,c.fsalfirst)

struct MidpointConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::Midpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u); atmp = similar(u, uEltypeNoUnits)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  MidpointCache(u,uprev,k,tmp,atmp,fsalfirst)
end

alg_cache(alg::Midpoint,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = MidpointConstantCache()

struct RK4Cache{uType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k₂::rateType
  k₃::rateType
  k₄::rateType
  k::rateType
  tmp::uType
  atmp::uEltypeNoUnits
end

u_cache(c::RK4Cache) = (c.atmp,)
du_cache(c::RK4Cache) = (c.fsalfirst,c.k₂,c.k₃,c.k₄,c.k)

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


struct CarpenterKennedy2N54Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  tab::TabType
end

u_cache(c::CarpenterKennedy2N54Cache) = ()
du_cache(c::CarpenterKennedy2N54Cache) = (c.k,c.fsalfirst)

struct CarpenterKennedy2N54ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  A2::T
  A3::T
  A4::T
  A5::T
  B1::T
  B2::T
  B3::T
  B4::T
  B5::T
  c2::T2
  c3::T2
  c4::T2
  c5::T2

  function CarpenterKennedy2N54ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    A2 = T(-567301805773/1357537059087)
    A3 = T(-2404267990393/2016746695238)
    A4 = T(-3550918686646/2091501179385)
    A5 = T(-1275806237668/842570457699)
    B1 = T(1432997174477/9575080441755)
    B2 = T(5161836677717/13612068292357)
    B3 = T(1720146321549/2090206949498)
    B4 = T(3134564353537/4481467310338)
    B5 = T(2277821191437/14882151754819)
    c2 = T2(1432997174477/9575080441755)
    c3 = T2(2526269341429/6820363962896)
    c4 = T2(2006345519317/3224310063776)
    c5 = T2(2802321613138/2924317926251)
    new{T,T2}(A2, A3, A4, A5, B1, B2, B3, B4, B5, c2, c3, c4, c5)
  end
end

function alg_cache(alg::CarpenterKennedy2N54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = CarpenterKennedy2N54ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  CarpenterKennedy2N54Cache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::CarpenterKennedy2N54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  CarpenterKennedy2N54ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
end



struct BS3Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::BS3Cache) = (c.atmp,c.utilde)
du_cache(c::BS3Cache) = (c.fsalfirst,c.k2,c.k3,c.k4)

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

struct OwrenZen3Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::OwrenZen3Cache) = (c.atmp,c.utilde)
du_cache(c::OwrenZen3Cache) = (c.k1,c.k2,c.k3,c.k4)

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

struct OwrenZen4Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::OwrenZen4Cache) = (c.atmp,c.utilde)
du_cache(c::OwrenZen4Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6)

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

struct OwrenZen5Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::OwrenZen5Cache) = (c.atmp,c.utilde)
du_cache(c::OwrenZen5Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8)

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

struct BS5Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::BS5Cache) = (c.utilde,c.atmp)
du_cache(c::BS5Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8)

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

struct Tsit5Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Tsit5Cache) = (c.utilde,c.atmp)
du_cache(c::Tsit5Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7)

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

struct DP5Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DP5Cache) = (c.utilde,c.atmp)
du_cache(c::DP5Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.dense_tmp3,c.dense_tmp4,c.update,c.bspl)

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

struct DP5ThreadedCache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DP5ThreadedCache) = (c.utilde,c.atmp)
du_cache(c::DP5ThreadedCache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.dense_tmp3,c.dense_tmp4,c.update,c.bspl)

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

struct Anas5Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k1::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Anas5Cache) = (c.atmp,c.utilde)
du_cache(c::Anas5Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7)

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
