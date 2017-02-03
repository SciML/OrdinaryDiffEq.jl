abstract OrdinaryDiffEqCache <: DECache
abstract OrdinaryDiffEqConstantCache <: OrdinaryDiffEqCache
abstract OrdinaryDiffEqMutableCache <: OrdinaryDiffEqCache
immutable ODEEmptyCache <: OrdinaryDiffEqConstantCache end
immutable ODEChunkCache{CS} <: OrdinaryDiffEqConstantCache end

type CompositeCache{T,F} <: OrdinaryDiffEqCache
  caches::T
  choice_function::F
  current::Int
end

function alg_cache{T}(alg::CompositeAlgorithm,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{T}})
  caches = map((x)->alg_cache(x,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,Val{T}),alg.algs)
  CompositeCache(caches,alg.choice_function,1)
end

alg_cache{F}(alg::OrdinaryDiffEqAlgorithm,prob,callback::F) = ODEEmptyCache()

immutable DiscreteCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du::rateType
end

u_cache(c::DiscreteCache) = (c.u,c.uprev)
du_cache(c::DiscreteCache) = (c.du)

function alg_cache(alg::Discrete,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  DiscreteCache(u,uprev,rate_prototype)
end

immutable DiscreteConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Discrete,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = DiscreteConstantCache()

immutable EulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
end

function alg_cache(alg::Euler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  EulerCache(u,uprev,similar(u),similar(rate_prototype),similar(rate_prototype))
end

u_cache(c::EulerCache) = (c.u,c.uprev,c.tmp)
du_cache(c::EulerCache) = (c.k,c.fsalfirst)

immutable EulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Euler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = EulerConstantCache()

immutable MidpointCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  du::rateType
  tmp::uType
  fsalfirst::rateType
end

u_cache(c::MidpointCache) = (c.u,c.uprev,c.tmp)
du_cache(c::MidpointCache) = (c.k,c.du,c.fsalfirst)

immutable MidpointConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::Midpoint,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = similar(rate_prototype)
  du = similar(rate_prototype)
  fsalfirst = similar(rate_prototype)
  MidpointCache(u,uprev,k,du,tmp,fsalfirst)
end

alg_cache(alg::Midpoint,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = MidpointConstantCache()

immutable RK4Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k₂::rateType
  k₃::rateType
  k₄::rateType
  k::rateType
  tmp::uType
end

u_cache(c::RK4Cache) = (c.u,c.uprev,c.tmp)
du_cache(c::RK4Cache) = (c.fsalfirst,c.k₂,c.k₃,c.k₄,c.k)

immutable RK4ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::RK4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k₁ = similar(rate_prototype)
  k₂ = similar(rate_prototype)
  k₃ = similar(rate_prototype)
  k₄ = similar(rate_prototype)
  k  = similar(rate_prototype)
  tmp = similar(u)
  RK4Cache(u,uprev,k₁,k₂,k₃,k₄,k,tmp)
end

alg_cache(alg::RK4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = RK4ConstantCache()

immutable BS3Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::BS3Cache) = (c.u,c.uprev,c.tmp,c.atmp,c.utilde)
du_cache(c::BS3Cache) = (c.fsalfirst,c.k2,c.k3,c.k4)

function alg_cache(alg::BS3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = BS3ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype)
  k3 = similar(rate_prototype)
  k4 = similar(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  BS3Cache(u,uprev,k1,k2,k3,k4,utilde,tmp,atmp,tab)
end

alg_cache(alg::BS3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = BS3ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable BS5Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  uhat::uType
  tmp::uType
  atmp::uEltypeNoUnits
  atmptilde::uEltypeNoUnits
  tab::TabType
end

u_cache(c::BS5Cache) = (c.u,c.uprev,c.utilde,c.uhat,c.tmp,c.atmp,c.atmptilde)
du_cache(c::BS5Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8)

function alg_cache(alg::BS5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = BS5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype)
  k3 = similar(rate_prototype)
  k4 = similar(rate_prototype)
  k5 = similar(rate_prototype)
  k6 = similar(rate_prototype)
  k7 = similar(rate_prototype)
  k8 = similar(rate_prototype)
  utilde = similar(u); uhat = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  atmptilde = similar(u,uEltypeNoUnits)
  BS5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,utilde,uhat,tmp,atmp,atmptilde,tab)
end

alg_cache(alg::BS5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = BS5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Tsit5Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Tsit5Cache) = (c.u,c.uprev,c.utilde,c.tmp,c.atmp)
du_cache(c::Tsit5Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7)

function alg_cache(alg::Tsit5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Tsit5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype)
  k3 = similar(rate_prototype)
  k4 = similar(rate_prototype)
  k5 = similar(rate_prototype)
  k6 = similar(rate_prototype)
  k7 = similar(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  atmptilde = similar(u,uEltypeNoUnits)
  Tsit5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp,tab)
end

alg_cache(alg::Tsit5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Tsit5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable DP5Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DP5Cache) = (c.u,c.uprev,c.utilde,c.tmp,c.atmp)
du_cache(c::DP5Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.dense_tmp3,c.dense_tmp4,c.update,c.bspl)

function alg_cache(alg::DP5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype)
  k3 = similar(rate_prototype)
  k4 = similar(rate_prototype)
  k5 = similar(rate_prototype)
  k6 = similar(rate_prototype)
  k7 = similar(rate_prototype)
  dense_tmp3 = similar(rate_prototype)
  dense_tmp4 = similar(rate_prototype)
  update = similar(rate_prototype)
  bspl = similar(rate_prototype)
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  tab = DP5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  DP5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp,tab)
end

alg_cache(alg::DP5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = DP5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable DP5ThreadedCache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DP5ThreadedCache) = (c.u,c.uprev,c.utilde,c.tmp,c.atmp)
du_cache(c::DP5ThreadedCache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.dense_tmp3,c.dense_tmp4,c.update,c.bspl)

function alg_cache(alg::DP5Threaded,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype)
  k3 = similar(rate_prototype)
  k4 = similar(rate_prototype)
  k5 = similar(rate_prototype)
  k6 = similar(rate_prototype)
  k7 = similar(rate_prototype)
  dense_tmp3 = similar(rate_prototype)
  dense_tmp4 = similar(rate_prototype)
  update = similar(rate_prototype)
  bspl = similar(rate_prototype)
  utilde = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  tab = DP5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  DP5ThreadedCache(u,uprev,k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp,tab)
end

immutable Vern6Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  k9::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern6Cache) = (c.u,c.uprev,c.utilde,c.tmp,c.atmp)
du_cache(c::Vern6Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9)

function alg_cache(alg::Vern6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern6ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype)
  k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype);
  k5 = similar(rate_prototype); k6 = similar(rate_prototype); k7 = similar(rate_prototype);
  k8 = similar(rate_prototype); k9 = similar(rate_prototype);
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits);
  Vern6Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,utilde,tmp,atmp,tab)
end

alg_cache(alg::Vern6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern6ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Vern7Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  k9::rateType
  k10::rateType
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern7Cache) = (c.u,c.uprev,c.utilde,c.update,c.tmp,c.atmp)
du_cache(c::Vern7Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10)

function alg_cache(alg::Vern7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype);
  k5 = similar(rate_prototype); k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype);
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); utilde = similar(u); update = similar(u)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits)
  Vern7Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))


immutable Vern8Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern8Cache) = (c.u,c.uprev,c.utilde,c.update,c.tmp,c.atmp)
du_cache(c::Vern8Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13)

function alg_cache(alg::Vern8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype);
  k4 = similar(rate_prototype);
  k5 = similar(rate_prototype); k6 = similar(rate_prototype); k7 = similar(rate_prototype);
  k8 = similar(rate_prototype); tmp = similar(u)
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); k11 = similar(rate_prototype);
  k12 = similar(rate_prototype); k13 = similar(rate_prototype)
  utilde = similar(u); update = similar(u);
  atmp = similar(u,uEltypeNoUnits)
  Vern8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Vern9Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern9Cache) = (c.u,c.uprev,c.utilde,c.update,c.tmp,c.atmp)
du_cache(c::Vern9Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16)

function alg_cache(alg::Vern9,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern9ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype); k2 = similar(rate_prototype);k3 = similar(rate_prototype);
  k4 = similar(rate_prototype);
  k5 = similar(rate_prototype); k6 = similar(rate_prototype);k7 = similar(rate_prototype);
  k8 = similar(rate_prototype);
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); k11 = similar(rate_prototype);
  k12 = similar(rate_prototype); update = similar(u)
  k13 = similar(rate_prototype); k14 = similar(rate_prototype); k15 = similar(rate_prototype);
  k16 =similar(rate_prototype);
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits);
  Vern9Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern9,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern9ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable TanYam7Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

u_cache(c::TanYam7Cache) = (c.u,c.uprev,c.utilde,c.tmp,c.atmp)
du_cache(c::TanYam7Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k)

function alg_cache(alg::TanYam7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = TanYam7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype); k2 = similar(rate_prototype) ; k3 = similar(rate_prototype); k4 = similar(rate_prototype)
  k5 = similar(rate_prototype); k6 = similar(rate_prototype) ; k7 = similar(rate_prototype); k8 = similar(rate_prototype)
  k9 = similar(rate_prototype); k10= similar(rate_prototype) ;
  utilde = similar(u); tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype)
  TanYam7Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,tmp,atmp,k,tab)
end

alg_cache(alg::TanYam7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = TanYam7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable DP8Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  kupdate::rateType
  update::uType
  udiff::rateType
  bspl::rateType
  dense_tmp3::rateType
  dense_tmp4::rateType
  dense_tmp5::rateType
  dense_tmp6::rateType
  dense_tmp7::rateType
  utilde::uType
  tmp::uType
  atmp::uEltypeNoUnits
  atmp2::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DP8Cache) = (c.u,c.uprev,c.update,c.utilde,c.tmp,c.atmp,c.atmp2)
du_cache(c::DP8Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.kupdate,c.udiff,c.bspl,c.dense_tmp3,c.dense_tmp4,c.dense_tmp5,c.dense_tmp6,c.dense_tmp7)

function alg_cache(alg::DP8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k1 = similar(rate_prototype); k2  = similar(rate_prototype); k3  = similar(rate_prototype);  k4 = similar(rate_prototype)
  k5 = similar(rate_prototype); k6  = similar(rate_prototype); k7  = similar(rate_prototype);  k8 = similar(rate_prototype)
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); k11 = similar(rate_prototype); k12 = similar(rate_prototype)
  kupdate = similar(rate_prototype); utilde = similar(u);
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); atmp2 = similar(u,uEltypeNoUnits); update = similar(u)
  k13 = similar(rate_prototype)
  k14 = similar(rate_prototype)
  k15 = similar(rate_prototype)
  k16 = similar(rate_prototype)
  udiff = similar(rate_prototype)
  bspl = similar(rate_prototype)
  # dense_tmp1 = udiff
  # dense_tmp2 = bspl
  dense_tmp3 = similar(rate_prototype)
  dense_tmp4 = similar(rate_prototype)
  dense_tmp5 = similar(rate_prototype)
  dense_tmp6 = similar(rate_prototype)
  dense_tmp7 = similar(rate_prototype)
  tab = DP8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  DP8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,kupdate,
           update,udiff,bspl,dense_tmp3,dense_tmp4,dense_tmp5,dense_tmp6,dense_tmp7,
           utilde,tmp,atmp,atmp2,tab)
end

alg_cache(alg::DP8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = DP8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable TsitPap8Cache{uType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  utilde::uType
  update::uType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

u_cache(c::TsitPap8Cache) = (c.u,c.uprev,c.utilde,c.tmp,c.atmp,c.update)
du_cache(c::TsitPap8Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k)

function alg_cache(alg::TsitPap8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = TsitPap8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype)
  k5 = similar(rate_prototype); k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype)
  k9 = similar(rate_prototype); k10 = similar(rate_prototype); k11 = similar(rate_prototype); k12 = similar(rate_prototype)
  k13 = similar(rate_prototype); update = similar(u); utilde = similar(u); k = similar(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits);
  TsitPap8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,update,tmp,atmp,k,tab)
end

alg_cache(alg::TsitPap8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = TsitPap8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable ExplicitRKCache{uType,rateType,uEltypeNoUnits,ksEltype,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  utilde::rateType
  uEEst::rateType
  atmp::uEltypeNoUnits
  fsalfirst::ksEltype
  fsallast::ksEltype
  kk::Vector{ksEltype}
  tab::TabType
end

u_cache(c::ExplicitRKCache) = (c.u,c.uprev,c.utilde,c.tmp,c.atmp,c.update)
du_cache(c::ExplicitRKCache) = (c.uEEst,c.kk...)

function alg_cache(alg::ExplicitRK,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  kk = Vector{typeof(rate_prototype)}(0)
  for i = 1:alg.tableau.stages
    push!(kk,similar(rate_prototype))
  end
  fsalfirst = kk[1]
  if isfsal(alg.tableau)
    fsallast = kk[end]
  else
    fsallast = similar(rate_prototype)
  end
  utilde = similar(rate_prototype)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  uEEst = similar(rate_prototype)
  tab = ExplicitRKConstantCache(alg.tableau,rate_prototype)
  ExplicitRKCache(u,uprev,tmp,utilde,uEEst,atmp,fsalfirst,fsallast,kk,tab)
end

immutable ExplicitRKConstantCache{MType,VType,KType} <: OrdinaryDiffEqConstantCache
  A::MType
  c::VType
  α::VType
  αEEst::VType
  stages::Int
  kk::KType
end

function ExplicitRKConstantCache(tableau,rate_prototype)
  @unpack A,c,α,αEEst,stages = tableau
  A = A' # Transpose A to column major looping
  kk = Array{typeof(rate_prototype)}(stages) # Not ks since that's for integrator.opts.dense
  ExplicitRKConstantCache(A,c,α,αEEst,stages,kk)
end

alg_cache(alg::ExplicitRK,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = ExplicitRKConstantCache(alg.tableau,rate_prototype)

immutable Feagin10Cache{uType,uEltypeNoUnits,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  k17::rateType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

u_cache(c::Feagin10Cache) = (c.u,c.uprev,c.tmp,c.atmp)
du_cache(c::Feagin10Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.k17,c.k)

function alg_cache(alg::Feagin10,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Feagin10ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype); k5 = similar(rate_prototype)
  k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype); k9 = similar(rate_prototype); k10 = similar(rate_prototype)
  k11 = similar(rate_prototype); k12 = similar(rate_prototype); k13 = similar(rate_prototype); k14 = similar(rate_prototype)
  k15 = similar(rate_prototype); k16 = similar(rate_prototype); k17 = similar(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype)

  Feagin10Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,tmp,atmp,k,tab)
end

alg_cache(alg::Feagin10,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Feagin10ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Feagin12Cache{uType,uEltypeNoUnits,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  k17::rateType
  k18::rateType
  k19::rateType
  k20::rateType
  k21::rateType
  k22::rateType
  k23::rateType
  k24::rateType
  k25::rateType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

u_cache(c::Feagin12Cache) = (c.u,c.uprev,c.tmp,c.atmp)
du_cache(c::Feagin12Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.k17,c.k18,c.k19,c.k20,c.k21,c.k22,c.k23,c.k24,c.k25,c.k)

function alg_cache(alg::Feagin12,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Feagin12ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype); k5 = similar(rate_prototype)
  k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype); k9 = similar(rate_prototype); k10 = similar(rate_prototype)
  k11 = similar(rate_prototype); k12 = similar(rate_prototype); k13 = similar(rate_prototype); k14 = similar(rate_prototype)
  k15 = similar(rate_prototype); k16 = similar(rate_prototype); k17 = similar(rate_prototype); k18 = similar(rate_prototype)
  k19 = similar(rate_prototype); k20 = similar(rate_prototype); k21 = similar(rate_prototype); k22 = similar(rate_prototype)
  k23 = similar(rate_prototype); k24 = similar(rate_prototype); k25 = similar(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype)

  Feagin12Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,tmp,atmp,k,tab)
end

alg_cache(alg::Feagin12,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Feagin12ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))


immutable Feagin14Cache{uType,uEltypeNoUnits,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k7::rateType
  k8::rateType
  k9::rateType
  k10::rateType
  k11::rateType
  k12::rateType
  k13::rateType
  k14::rateType
  k15::rateType
  k16::rateType
  k17::rateType
  k18::rateType
  k19::rateType
  k20::rateType
  k21::rateType
  k22::rateType
  k23::rateType
  k24::rateType
  k25::rateType
  k26::rateType
  k27::rateType
  k28::rateType
  k29::rateType
  k30::rateType
  k31::rateType
  k32::rateType
  k33::rateType
  k34::rateType
  k35::rateType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

u_cache(c::Feagin14Cache) = (c.u,c.uprev,c.tmp,c.atmp)
du_cache(c::Feagin14Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.k17,c.k18,c.k19,c.k20,c.k21,c.k22,c.k23,c.k24,c.k25,c.k26,c.k27,c.k28,c.k29,c.k30,c.k31,c.k32,c.k33,c.k34,c.k35,c.k)


function alg_cache(alg::Feagin14,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Feagin14ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = similar(rate_prototype); k2 = similar(rate_prototype); k3 = similar(rate_prototype); k4 = similar(rate_prototype); k5 = similar(rate_prototype)
  k6 = similar(rate_prototype); k7 = similar(rate_prototype); k8 = similar(rate_prototype); k9 = similar(rate_prototype); k10 = similar(rate_prototype)
  k11 = similar(rate_prototype); k12 = similar(rate_prototype); k13 = similar(rate_prototype); k14 = similar(rate_prototype)
  k15 = similar(rate_prototype); k16 = similar(rate_prototype); k17 = similar(rate_prototype); k18 = similar(rate_prototype)
  k19 = similar(rate_prototype); k20 = similar(rate_prototype); k21 = similar(rate_prototype); k22 = similar(rate_prototype)
  k23 = similar(rate_prototype); k24 = similar(rate_prototype); k25 = similar(rate_prototype)
  k26 = similar(rate_prototype); k27 = similar(rate_prototype); k28 = similar(rate_prototype)
  k29 = similar(rate_prototype); k30 = similar(rate_prototype); k31 = similar(rate_prototype); k32 = similar(rate_prototype)
  k33 = similar(rate_prototype); k34 = similar(rate_prototype); k35 = similar(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits); k = similar(rate_prototype)

  Feagin14Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,
                k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,
                k31,k32,k33,k34,k35,tmp,atmp,k,tab)
end

alg_cache(alg::Feagin14,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Feagin14ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))


type Rosenbrock23Cache{uType,rateType,vecuType,JType,TabType,TFType,UFType,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du1::rateType
  du2::rateType
  f₁::rateType
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uType
  J::JType
  W::JType
  tmp::uType
  tmp2::uType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::vecuType
  factorization::F
end

u_cache(c::Rosenbrock23Cache) = (c.u,c.uprev,c.tmp,c.tmp2,c.dT)
du_cache(c::Rosenbrock23Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock23Cache) = (c.J,c.W)
vecu_cache(c::Rosenbrock23Cache) = (c.vectmp,c.vectmp2,c.vectmp3)

type Rosenbrock32Cache{uType,rateType,vecuType,JType,TabType,TFType,UFType,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du1::rateType
  du2::rateType
  f₁::rateType
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uType
  J::JType
  W::JType
  tmp::uType
  tmp2::uType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::vecuType
  factorization::F
end

u_cache(c::Rosenbrock32Cache) = (c.u,c.uprev,c.tmp,c.tmp2,c.dT)
du_cache(c::Rosenbrock32Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock32Cache) = (c.J,c.W)
vecu_cache(c::Rosenbrock32Cache) = (c.vectmp,c.vectmp2,c.vectmp3)

function alg_cache(alg::Rosenbrock23,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k₁ = similar(rate_prototype)
  k₂ = similar(rate_prototype)
  k₃ = similar(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  # f₀ = similar(u) fsalfirst
  f₁ = similar(rate_prototype)
  vectmp = similar(vec(u))
  vectmp2 = similar(vec(u))
  vectmp3 = similar(vec(u))
  fsalfirst = similar(rate_prototype)
  fsallast = similar(rate_prototype)
  dT = similar(u)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J); tmp2 = similar(u)
  tmp = reshape(vectmp2,size(u)...)
  tab = Rosenbrock23ConstantCache(uEltypeNoUnits,identity,identity)
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev,du2)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(vec(u))
  Rosenbrock23Cache(u,uprev,k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,tmp2,tab,tf,uf,linsolve_tmp,alg.factorization)
end

function alg_cache(alg::Rosenbrock32,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k₁ = similar(rate_prototype)
  k₂ = similar(rate_prototype)
  k₃ = similar(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  # f₀ = similar(u) fsalfirst
  f₁ = similar(rate_prototype)
  vectmp = similar(vec(u))
  vectmp2 = similar(vec(u))
  vectmp3 = similar(vec(u))
  fsalfirst = similar(rate_prototype)
  fsallast = similar(rate_prototype)
  dT = similar(u)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J); tmp2 = similar(u)
  tmp = reshape(vectmp2,size(u)...)
  tab = Rosenbrock32ConstantCache(uEltypeNoUnits,identity,identity)
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev,du2)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(vec(u))
  Rosenbrock32Cache(u,uprev,k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,tmp2,tab,tf,uf,linsolve_tmp,alg.factorization)
end

immutable Rosenbrock23ConstantCache{T,TF,UF} <: OrdinaryDiffEqConstantCache
  c₃₂::T
  d::T
  tf::TF
  uf::UF
end

function Rosenbrock23ConstantCache(T::Type,tf,uf)
  c₃₂ = T(6 + sqrt(2))
  d = T(1/(2+sqrt(2)))
  Rosenbrock23ConstantCache(c₃₂,d,tf,uf)
end

function alg_cache(alg::Rosenbrock23,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock23ConstantCache(uEltypeNoUnits,tf,uf)
end

immutable Rosenbrock32ConstantCache{T,TF,UF} <: OrdinaryDiffEqConstantCache
  c₃₂::T
  d::T
  tf::TF
  uf::UF
end

function Rosenbrock32ConstantCache(T::Type,tf,uf)
  c₃₂ = T(6 + sqrt(2))
  d = T(1/(2+sqrt(2)))
  Rosenbrock32ConstantCache(c₃₂,d,tf,uf)
end

function alg_cache(alg::Rosenbrock32,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  tf = TimeDerivativeWrapper(f,u)
  uf = UDerivativeWrapper(f,t)
  Rosenbrock32ConstantCache(uEltypeNoUnits,tf,uf)
end

type ImplicitEulerCache{uType,vecuType,DiffCacheType,rateType,rhsType,adfType,CS} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  u_old::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  rhs::rhsType
  adf::adfType
end

u_cache(c::ImplicitEulerCache)    = (c.u,c.uprev,c.uprev2,c.u_old)
du_cache(c::ImplicitEulerCache)   = (c.k,c.fsalfirst)
vecu_cache(c::ImplicitEulerCache) = (c.uhold,)
dual_cache(c::ImplicitEulerCache) = (c.dual_cache,)

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  u_old = similar(u); k = similar(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,alg)})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  rhs = RHS_IE(f,u_old,t,t,dual_cache,size(u),eachindex(u))
  fsalfirst = similar(rate_prototype)
  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  else
    adf = nothing
  end
  tmp = u_old
  ImplicitEulerCache{typeof(u),typeof(uhold),typeof(dual_cache),typeof(k),
                     typeof(rhs),typeof(adf),determine_chunksize(u,alg)}(
                     u,uprev,uprev2,uhold,dual_cache,u_old,tmp,k,fsalfirst,rhs,adf)
end

immutable ImplicitEulerConstantCache{vecuType,rhsType,adfType,CS} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  u_old::vecuType
  rhs::rhsType
  adf::adfType
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  u_old = Vector{typeof(u)}(1)
  rhs = RHS_IE_Scalar(f,u_old,t,t)
  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  else
    adf = nothing
  end
  ImplicitEulerConstantCache{typeof(uhold),typeof(rhs),typeof(adf),1}(uhold,u_old,rhs,adf)
end

type TrapezoidCache{uType,vecuType,DiffCacheType,rateType,rhsType,adfType,CS} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  uhold::vecuType
  u_old::uType
  fsalfirst::rateType
  dual_cache::DiffCacheType
  tmp::uType
  k::rateType
  rhs::rhsType
  adf::adfType
end

u_cache(c::TrapezoidCache)    = (c.u,c.uprev,c.uprev2,c.u_old)
du_cache(c::TrapezoidCache)   = (c.k,c.fsalfirst)
vecu_cache(c::TrapezoidCache) = (c.uhold,)
dual_cache(c::TrapezoidCache) = (c.dual_cache,)

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  u_old = similar(u); k = similar(rate_prototype)
  uhold = vec(u); fsalfirst = similar(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,alg)})
  rhs = RHS_Trap(f,u_old,fsalfirst,t,t,size(u),dual_cache,eachindex(u))
  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  else
    adf = nothing
  end
  tmp = u_old
  TrapezoidCache{typeof(u),typeof(uhold),typeof(dual_cache),typeof(k),
    typeof(rhs),typeof(adf),determine_chunksize(u,alg)}(
    u,uprev,uprev2,uhold,u_old,fsalfirst,dual_cache,tmp,k,rhs,adf)
end


immutable TrapezoidConstantCache{vecuType,rhsType,adfType,CS} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  u_old::vecuType
  rhs::rhsType
  adf::adfType
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  u_old = Vector{typeof(u)}(1)
  rhs = RHS_Trap_Scalar(f,u_old,rate_prototype,t,t)
  if alg_autodiff(alg)
    adf = autodiff_setup(rhs,uhold,alg)
  else
    adf = nothing
  end
  TrapezoidConstantCache{typeof(uhold),typeof(rhs),typeof(adf),1}(uhold,u_old,rhs,adf)
end



get_chunksize(cache::DECache) = error("This cache does not have a chunksize.")
get_chunksize{uType,DiffCacheType,rateType,CS}(cache::ImplicitEulerCache{uType,DiffCacheType,rateType,CS}) = CS
get_chunksize{uType,DiffCacheType,rateType,CS}(cache::TrapezoidCache{uType,DiffCacheType,rateType,CS}) = CS
get_chunksize{CS}(cache::ODEChunkCache{CS}) = CS
