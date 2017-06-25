@compat abstract type OrdinaryDiffEqCache <: DECache end
@compat abstract type OrdinaryDiffEqConstantCache <: OrdinaryDiffEqCache end
@compat abstract type OrdinaryDiffEqMutableCache <: OrdinaryDiffEqCache end
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

u_cache(c::DiscreteCache) = ()
du_cache(c::DiscreteCache) = (c.du)

function alg_cache(alg::Discrete,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  DiscreteCache(u,uprev,discrete_scale_by_time(alg) ? rate_prototype : similar(u))
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
  EulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype))
end

u_cache(c::EulerCache) = ()
du_cache(c::EulerCache) = (c.k,c.fsalfirst)

immutable EulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Euler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = EulerConstantCache()

immutable SplitEulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
end

function alg_cache(alg::SplitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  SplitEulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype))
end

u_cache(c::SplitEulerCache) = ()
du_cache(c::SplitEulerCache) = (c.k,c.fsalfirst)

immutable SplitEulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::SplitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = SplitEulerConstantCache()

immutable SymplecticEulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
end

function alg_cache(alg::SymplecticEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  SymplecticEulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype))
end

u_cache(c::SymplecticEulerCache) = ()
du_cache(c::SymplecticEulerCache) = (c.k,c.fsalfirst)

immutable SymplecticEulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::SymplecticEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = SymplecticEulerConstantCache()

immutable MidpointCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
end

u_cache(c::MidpointCache) = ()
du_cache(c::MidpointCache) = (c.k,c.fsalfirst)

immutable MidpointConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::Midpoint,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  MidpointCache(u,uprev,k,tmp,fsalfirst)
end

alg_cache(alg::Midpoint,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = MidpointConstantCache()


immutable SSPRK22Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
end

u_cache(c::SSPRK22Cache) = ()
du_cache(c::SSPRK22Cache) = (c.k,c.fsalfirst)

immutable SSPRK22ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK22,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  SSPRK22Cache(u,uprev,k,tmp,fsalfirst)
end

alg_cache(alg::SSPRK22,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = SSPRK22ConstantCache()


immutable SSPRK33Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
end

u_cache(c::SSPRK33Cache) = ()
du_cache(c::SSPRK33Cache) = (c.k,c.fsalfirst)

immutable SSPRK33ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK33,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  SSPRK33Cache(u,uprev,k,tmp,fsalfirst)
end

alg_cache(alg::SSPRK33,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = SSPRK33ConstantCache()


immutable SSPRK432Cache{uType,rateType,uArrayType,uEltypeNoUnits} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  utilde::uArrayType
  atmp::uEltypeNoUnits
end

u_cache(c::SSPRK432Cache) = (c.utilde,c.atmp)
du_cache(c::SSPRK432Cache) = (c.k,c.fsalfirst)

immutable SSPRK432ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK432,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u))
  SSPRK432Cache(u,uprev,k,tmp,fsalfirst,utilde,atmp)
end

alg_cache(alg::SSPRK432,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = SSPRK432ConstantCache()


immutable SSPRK104Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₄::rateType
  tmp::uType
  fsalfirst::rateType
end

u_cache(c::SSPRK104Cache) = ()
du_cache(c::SSPRK104Cache) = (c.k,c.fsalfirst,c.k₄)

immutable SSPRK104ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK104,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  k₄ = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  SSPRK104Cache(u,uprev,k,k₄,tmp,fsalfirst)
end

alg_cache(alg::SSPRK104,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = SSPRK104ConstantCache()


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

u_cache(c::RK4Cache) = ()
du_cache(c::RK4Cache) = (c.fsalfirst,c.k₂,c.k₃,c.k₄,c.k)

immutable RK4ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::RK4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  k₄ = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  RK4Cache(u,uprev,k₁,k₂,k₃,k₄,k,tmp)
end

alg_cache(alg::RK4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = RK4ConstantCache()

immutable BS3Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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

function alg_cache(alg::BS3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = BS3ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u))
  tmp = similar(u)
  BS3Cache(u,uprev,k1,k2,k3,k4,utilde,tmp,atmp,tab)
end

alg_cache(alg::BS3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = BS3ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable BS5Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  uhat::uType
  tmp::uType
  atmp::uEltypeNoUnits
  atmptilde::uEltypeNoUnits
  tab::TabType
end

u_cache(c::BS5Cache) = (c.utilde,c.uhat,c.atmp,c.atmptilde)
du_cache(c::BS5Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8)

function alg_cache(alg::BS5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = BS5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype)
  k7 = zeros(rate_prototype)
  k8 = zeros(rate_prototype)
  utilde = similar(u,indices(u)); uhat = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u))
  tmp = similar(u)
  atmptilde = similar(u,uEltypeNoUnits,indices(u))
  BS5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,utilde,uhat,tmp,atmp,atmptilde,tab)
end

alg_cache(alg::BS5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = BS5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Tsit5Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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

function alg_cache(alg::Tsit5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Tsit5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype)
  k7 = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u))
  tmp = similar(u)
  atmptilde = similar(u,uEltypeNoUnits,indices(u))
  Tsit5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp,tab)
end

alg_cache(alg::Tsit5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Tsit5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable DP5Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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

function alg_cache(alg::DP5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype)
  k7 = zeros(rate_prototype)
  dense_tmp3 = zeros(rate_prototype)
  dense_tmp4 = zeros(rate_prototype)
  update = zeros(rate_prototype)
  bspl = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u))
  tab = DP5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  cache = DP5Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp,tab)
  cache
end

alg_cache(alg::DP5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = DP5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable DP5ThreadedCache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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

function alg_cache(alg::DP5Threaded,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype)
  k7 = zeros(rate_prototype)
  dense_tmp3 = zeros(rate_prototype)
  dense_tmp4 = zeros(rate_prototype)
  update = zeros(rate_prototype)
  bspl = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u))
  tab = DP5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  DP5ThreadedCache(u,uprev,k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp,tab)
end

immutable Vern6Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern6Cache) = (c.utilde,c.atmp)
du_cache(c::Vern6Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9)

function alg_cache(alg::Vern6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern6ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype);
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype); k7 = zeros(rate_prototype);
  k8 = zeros(rate_prototype); k9 = zeros(rate_prototype);
  utilde = similar(u,indices(u)); tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u));
  Vern6Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,utilde,tmp,atmp,tab)
end

alg_cache(alg::Vern6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern6ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Vern7Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  update::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern7Cache) = (c.utilde,c.update,c.atmp)
du_cache(c::Vern7Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10)

function alg_cache(alg::Vern7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype);
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype); k7 = zeros(rate_prototype); k8 = zeros(rate_prototype);
  k9 = zeros(rate_prototype); k10 = zeros(rate_prototype); utilde = similar(u,indices(u)); update = similar(u,indices(u))
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u))
  Vern7Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))


immutable Vern8Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  update::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern8Cache) = (c.utilde,c.update,c.atmp)
du_cache(c::Vern8Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13)

function alg_cache(alg::Vern8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype);
  k4 = zeros(rate_prototype);
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype); k7 = zeros(rate_prototype);
  k8 = zeros(rate_prototype); tmp = similar(u)
  k9 = zeros(rate_prototype); k10 = zeros(rate_prototype); k11 = zeros(rate_prototype);
  k12 = zeros(rate_prototype); k13 = zeros(rate_prototype)
  utilde = similar(u,indices(u)); update = similar(u,indices(u));
  atmp = similar(u,uEltypeNoUnits,indices(u))
  Vern8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable Vern9Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  update::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::Vern9Cache) = (c.utilde,c.update,c.atmp)
du_cache(c::Vern9Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16)

function alg_cache(alg::Vern9,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Vern9ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype);k3 = zeros(rate_prototype);
  k4 = zeros(rate_prototype);
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype);k7 = zeros(rate_prototype);
  k8 = zeros(rate_prototype);
  k9 = zeros(rate_prototype); k10 = zeros(rate_prototype); k11 = zeros(rate_prototype);
  k12 = zeros(rate_prototype); update = similar(u,indices(u))
  k13 = zeros(rate_prototype); k14 = zeros(rate_prototype); k15 = zeros(rate_prototype);
  k16 =zeros(rate_prototype);
  utilde = similar(u,indices(u)); tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u));
  Vern9Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,utilde,update,tmp,atmp,tab)
end

alg_cache(alg::Vern9,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Vern9ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable TanYam7Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

u_cache(c::TanYam7Cache) = (c.utilde,c.atmp)
du_cache(c::TanYam7Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k)

function alg_cache(alg::TanYam7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = TanYam7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype) ; k3 = zeros(rate_prototype); k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype) ; k7 = zeros(rate_prototype); k8 = zeros(rate_prototype)
  k9 = zeros(rate_prototype); k10= zeros(rate_prototype) ;
  utilde = similar(u,indices(u)); tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u)); k = zeros(rate_prototype)
  TanYam7Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,utilde,tmp,atmp,k,tab)
end

alg_cache(alg::TanYam7,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = TanYam7ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable DP8Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  update::uArrayType
  udiff::rateType
  bspl::rateType
  dense_tmp3::rateType
  dense_tmp4::rateType
  dense_tmp5::rateType
  dense_tmp6::rateType
  dense_tmp7::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  atmp2::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DP8Cache) = (c.update,c.utilde,c.atmp,c.atmp2)
du_cache(c::DP8Cache) = (c.k1,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.kupdate,c.udiff,c.bspl,c.dense_tmp3,c.dense_tmp4,c.dense_tmp5,c.dense_tmp6,c.dense_tmp7)

function alg_cache(alg::DP8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k1 = zeros(rate_prototype); k2  = zeros(rate_prototype); k3  = zeros(rate_prototype);  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype); k6  = zeros(rate_prototype); k7  = zeros(rate_prototype);  k8 = zeros(rate_prototype)
  k9 = zeros(rate_prototype); k10 = zeros(rate_prototype); k11 = zeros(rate_prototype); k12 = zeros(rate_prototype)
  kupdate = zeros(rate_prototype); utilde = similar(u,indices(u));
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u)); atmp2 = similar(u,uEltypeNoUnits,indices(u)); update = similar(u,indices(u))
  k13 = zeros(rate_prototype)
  k14 = zeros(rate_prototype)
  k15 = zeros(rate_prototype)
  k16 = zeros(rate_prototype)
  udiff = zeros(rate_prototype)
  bspl = zeros(rate_prototype)
  # dense_tmp1 = udiff
  # dense_tmp2 = bspl
  dense_tmp3 = zeros(rate_prototype)
  dense_tmp4 = zeros(rate_prototype)
  dense_tmp5 = zeros(rate_prototype)
  dense_tmp6 = zeros(rate_prototype)
  dense_tmp7 = zeros(rate_prototype)
  tab = DP8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  DP8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,kupdate,
           update,udiff,bspl,dense_tmp3,dense_tmp4,dense_tmp5,dense_tmp6,dense_tmp7,
           utilde,tmp,atmp,atmp2,tab)
end

alg_cache(alg::DP8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = DP8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))

immutable TsitPap8Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  utilde::uArrayType
  update::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  k::rateType
  tab::TabType
end

u_cache(c::TsitPap8Cache) = (c.utilde,c.atmp,c.update)
du_cache(c::TsitPap8Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k)

function alg_cache(alg::TsitPap8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = TsitPap8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype); k6 = zeros(rate_prototype); k7 = zeros(rate_prototype); k8 = zeros(rate_prototype)
  k9 = zeros(rate_prototype); k10 = zeros(rate_prototype); k11 = zeros(rate_prototype); k12 = zeros(rate_prototype)
  k13 = zeros(rate_prototype); update = similar(u,indices(u)); utilde = similar(u,indices(u)); k = zeros(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u));
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

u_cache(c::ExplicitRKCache) = (c.utilde,c.atmp,c.update)
du_cache(c::ExplicitRKCache) = (c.uEEst,c.kk...)

function alg_cache(alg::ExplicitRK,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  kk = Vector{typeof(rate_prototype)}(0)
  for i = 1:alg.tableau.stages
    push!(kk,zeros(rate_prototype))
  end
  fsalfirst = kk[1]
  if isfsal(alg.tableau)
    fsallast = kk[end]
  else
    fsallast = zeros(rate_prototype)
  end
  utilde = zeros(rate_prototype)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits,indices(u))
  uEEst = zeros(rate_prototype)
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

u_cache(c::Feagin10Cache) = (c.atmp,)
du_cache(c::Feagin10Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.k17,c.k)

function alg_cache(alg::Feagin10,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Feagin10ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype); k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype); k7 = zeros(rate_prototype); k8 = zeros(rate_prototype); k9 = zeros(rate_prototype); k10 = zeros(rate_prototype)
  k11 = zeros(rate_prototype); k12 = zeros(rate_prototype); k13 = zeros(rate_prototype); k14 = zeros(rate_prototype)
  k15 = zeros(rate_prototype); k16 = zeros(rate_prototype); k17 = zeros(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u)); k = zeros(rate_prototype)

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

u_cache(c::Feagin12Cache) = (c.atmp,)
du_cache(c::Feagin12Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.k17,c.k18,c.k19,c.k20,c.k21,c.k22,c.k23,c.k24,c.k25,c.k)

function alg_cache(alg::Feagin12,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Feagin12ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype); k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype); k7 = zeros(rate_prototype); k8 = zeros(rate_prototype); k9 = zeros(rate_prototype); k10 = zeros(rate_prototype)
  k11 = zeros(rate_prototype); k12 = zeros(rate_prototype); k13 = zeros(rate_prototype); k14 = zeros(rate_prototype)
  k15 = zeros(rate_prototype); k16 = zeros(rate_prototype); k17 = zeros(rate_prototype); k18 = zeros(rate_prototype)
  k19 = zeros(rate_prototype); k20 = zeros(rate_prototype); k21 = zeros(rate_prototype); k22 = zeros(rate_prototype)
  k23 = zeros(rate_prototype); k24 = zeros(rate_prototype); k25 = zeros(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u)); k = zeros(rate_prototype)

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

u_cache(c::Feagin14Cache) = (c.atmp,)
du_cache(c::Feagin14Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.k17,c.k18,c.k19,c.k20,c.k21,c.k22,c.k23,c.k24,c.k25,c.k26,c.k27,c.k28,c.k29,c.k30,c.k31,c.k32,c.k33,c.k34,c.k35,c.k)


function alg_cache(alg::Feagin14,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tab = Feagin14ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  k1 = zeros(rate_prototype); k2 = zeros(rate_prototype); k3 = zeros(rate_prototype); k4 = zeros(rate_prototype); k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype); k7 = zeros(rate_prototype); k8 = zeros(rate_prototype); k9 = zeros(rate_prototype); k10 = zeros(rate_prototype)
  k11 = zeros(rate_prototype); k12 = zeros(rate_prototype); k13 = zeros(rate_prototype); k14 = zeros(rate_prototype)
  k15 = zeros(rate_prototype); k16 = zeros(rate_prototype); k17 = zeros(rate_prototype); k18 = zeros(rate_prototype)
  k19 = zeros(rate_prototype); k20 = zeros(rate_prototype); k21 = zeros(rate_prototype); k22 = zeros(rate_prototype)
  k23 = zeros(rate_prototype); k24 = zeros(rate_prototype); k25 = zeros(rate_prototype)
  k26 = zeros(rate_prototype); k27 = zeros(rate_prototype); k28 = zeros(rate_prototype)
  k29 = zeros(rate_prototype); k30 = zeros(rate_prototype); k31 = zeros(rate_prototype); k32 = zeros(rate_prototype)
  k33 = zeros(rate_prototype); k34 = zeros(rate_prototype); k35 = zeros(rate_prototype)
  tmp = similar(u); atmp = similar(u,uEltypeNoUnits,indices(u)); k = zeros(rate_prototype)

  Feagin14Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,
                k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,
                k31,k32,k33,k34,k35,tmp,atmp,k,tab)
end

alg_cache(alg::Feagin14,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Feagin14ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))


type Rosenbrock23Cache{uType,uArrayType,rateType,du2Type,LinuType,vecuType,JType,TabType,TFType,UFType,F,JCType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du1::rateType
  du2::du2Type
  f₁::rateType
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uArrayType
  J::JType
  W::JType
  tmp::uArrayType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::LinuType
  linsolve_tmp_vec::vecuType
  linsolve::F
  jac_config::JCType
end

u_cache(c::Rosenbrock23Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock23Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock23Cache) = (c.J,c.W)
vecu_cache(c::Rosenbrock23Cache) = (c.vectmp,c.vectmp2,c.vectmp3)

type Rosenbrock32Cache{uType,uArrayType,rateType,du2Type,LinuType,vecuType,JType,TabType,TFType,UFType,F,JCType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k₁::rateType
  k₂::rateType
  k₃::rateType
  du1::rateType
  du2::du2Type
  f₁::rateType
  vectmp::vecuType
  vectmp2::vecuType
  vectmp3::vecuType
  fsalfirst::rateType
  fsallast::rateType
  dT::uArrayType
  J::JType
  W::JType
  tmp::uArrayType
  tab::TabType
  tf::TFType
  uf::UFType
  linsolve_tmp::LinuType
  linsolve_tmp_vec::vecuType
  linsolve::F
  jac_config::JCType
end

u_cache(c::Rosenbrock32Cache) = (c.dT,c.tmp)
du_cache(c::Rosenbrock32Cache) = (c.k₁,c.k₂,c.k₃,c.du1,c.du2,c.f₁,c.fsalfirst,c.fsallast,c.linsolve_tmp)
jac_cache(c::Rosenbrock32Cache) = (c.J,c.W)
vecu_cache(c::Rosenbrock32Cache) = (c.vectmp,c.vectmp2,c.vectmp3)

function alg_cache(alg::Rosenbrock23,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  # f₀ = similar(u) fsalfirst
  f₁ = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J);
  tmp = similar(u,indices(u))
  tab = Rosenbrock23ConstantCache(uEltypeNoUnits,identity,identity)
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  Rosenbrock23Cache(u,uprev,k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,
                    fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,linsolve_tmp_vec,
                    alg.linsolve,jac_config)
end

function alg_cache(alg::Rosenbrock32,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  du1 = zeros(rate_prototype)
  du2 = zeros(rate_prototype)
  # f₀ = similar(u) fsalfirst
  f₁ = zeros(rate_prototype)
  vectmp = vec(similar(u,indices(u)))
  vectmp2 = vec(similar(u,indices(u)))
  vectmp3 = vec(similar(u,indices(u)))
  fsalfirst = zeros(rate_prototype)
  fsallast = zeros(rate_prototype)
  dT = similar(u,indices(u))
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J); tmp = similar(u,indices(u))
  tab = Rosenbrock32ConstantCache(uEltypeNoUnits,identity,identity)
  vf = VectorF(f,size(u))
  vfr = VectorFReturn(f,size(u))
  tf = TimeGradientWrapper(vf,uprev)
  uf = UJacobianWrapper(vfr,t)
  linsolve_tmp = similar(u,indices(u))
  linsolve_tmp_vec = vec(linsolve_tmp)
  jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  Rosenbrock32Cache(u,uprev,k₁,k₂,k₃,du1,du2,f₁,vectmp,vectmp2,vectmp3,fsalfirst,fsallast,dT,J,W,tmp,tab,tf,uf,linsolve_tmp,linsolve_tmp_vec,alg.linsolve,jac_config)
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

type ImplicitEulerCache{uType,uArrayType,vecuType,DiffCacheType,rateType,rhsType,nl_rhsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  u_old::uArrayType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

u_cache(c::ImplicitEulerCache)    = (c.uprev2,c.u_old)
du_cache(c::ImplicitEulerCache)   = (c.k,c.fsalfirst)
vecu_cache(c::ImplicitEulerCache) = (c.uhold,)
dual_cache(c::ImplicitEulerCache) = (c.dual_cache,)

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  u_old = vec(tmp); k = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  rhs = RHS_IE(f,u_old,t,t,dual_cache,size(u),eachindex(u))
  fsalfirst = zeros(rate_prototype)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)

  ImplicitEulerCache{typeof(u),typeof(u_old),typeof(uhold),typeof(dual_cache),
                     typeof(k),typeof(rhs),typeof(nl_rhs)}(
                     u,uprev,uprev2,uhold,dual_cache,u_old,tmp,k,fsalfirst,rhs,nl_rhs)
end

immutable ImplicitEulerConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  u_old::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

function alg_cache(alg::ImplicitEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  u_old = Vector{typeof(u)}(1)
  rhs = RHS_IE_Scalar(f,u_old,t,t)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  ImplicitEulerConstantCache{typeof(uhold),typeof(rhs),typeof(nl_rhs)}(uhold,u_old,rhs,nl_rhs)
end

type TrapezoidCache{uType,uArrayType,vecuType,DiffCacheType,rateType,rhsType,nl_rhsType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  uhold::vecuType
  u_old::uArrayType
  fsalfirst::rateType
  dual_cache::DiffCacheType
  tmp::uType
  k::rateType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

u_cache(c::TrapezoidCache)    = (c.uprev2,c.u_old)
du_cache(c::TrapezoidCache)   = (c.k,c.fsalfirst)
vecu_cache(c::TrapezoidCache) = (c.uhold,)
dual_cache(c::TrapezoidCache) = (c.dual_cache,)

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  u_old = vec(tmp); k = zeros(rate_prototype)
  uhold = vec(u); fsalfirst = zeros(rate_prototype)
  f_old = vec(fsalfirst)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  rhs = RHS_Trap(f,u_old,f_old,t,t,size(u),dual_cache,eachindex(u))
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  TrapezoidCache{typeof(u),typeof(u_old),typeof(uhold),typeof(dual_cache),typeof(k),
    typeof(rhs),typeof(nl_rhs)}(u,uprev,uprev2,uhold,u_old,fsalfirst,dual_cache,tmp,k,rhs,nl_rhs)
end


immutable TrapezoidConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  u_old::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

function alg_cache(alg::Trapezoid,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  u_old = Vector{typeof(u)}(1)
  rhs = RHS_Trap_Scalar(f,u_old,rate_prototype,t,t)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  TrapezoidConstantCache{typeof(uhold),typeof(rhs),typeof(nl_rhs)}(uhold,u_old,rhs,nl_rhs)
end

immutable IIF1ConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

immutable IIF1Cache{uType,vecuType,DiffCacheType,rhsType,nl_rhsType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  tmp::uType
  rhs::rhsType
  nl_rhs::nl_rhsType
  rtmp1::rateType
  fsalfirst::rateType
  k::rateType
end

u_cache(c::IIF1Cache)    = (c.uprev2,c.u_old)
du_cache(c::IIF1Cache)   = (c.rtmp1,c.tmp,c.fsalfirst,c.k)
vecu_cache(c::IIF1Cache) = (c.uhold,)
dual_cache(c::IIF1Cache) = (c.dual_cache,)

function alg_cache(alg::IIF1,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  tmp = zero(rate_prototype)
  rhs = RHS_IIF1_Scalar(f,tmp,t,t)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  IIF1ConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::IIF1,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u,indices(u)); rtmp1 = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  rhs = RHS_IIF1(f,tmp,t,t,dual_cache,size(u),eachindex(u))
  k = similar(rate_prototype); fsalfirst = similar(rate_prototype)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  IIF1Cache(u,uprev,uhold,dual_cache,tmp,rhs,nl_rhs,rtmp1,fsalfirst,k)
end

immutable IIF2ConstantCache{vecuType,rhsType,nl_rhsType} <: OrdinaryDiffEqConstantCache
  uhold::vecuType
  rhs::rhsType
  nl_rhs::nl_rhsType
end

immutable IIF2Cache{uType,vecuType,DiffCacheType,rhsType,nl_rhsType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uhold::vecuType
  dual_cache::DiffCacheType
  tmp::uType
  rhs::rhsType
  nl_rhs::nl_rhsType
  rtmp1::rateType
  fsalfirst::rateType
  k::rateType
end

u_cache(c::IIF2Cache)    = (c.uprev2,c.u_old)
du_cache(c::IIF2Cache)   = (c.rtmp1,c.tmp,c.fsalfirst,c.k)
vecu_cache(c::IIF2Cache) = (c.uhold,)
dual_cache(c::IIF2Cache) = (c.dual_cache,)

function alg_cache(alg::IIF2,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}})
  uhold = Vector{typeof(u)}(1)
  tmp = zero(rate_prototype)
  rhs = RHS_IIF2_Scalar(f,tmp,t,t)
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  IIF2ConstantCache(uhold,rhs,nl_rhs)
end

function alg_cache(alg::IIF2,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u,indices(u)); rtmp1 = zeros(rate_prototype)
  dual_cache = DiffCache(u,Val{determine_chunksize(u,get_chunksize(alg.nlsolve))})
  uhold = vec(u) # this makes uhold the same values as integrator.u
  k = similar(rate_prototype); fsalfirst = similar(rate_prototype)
  rhs = RHS_IIF2(f,tmp,t,t,dual_cache,size(u),eachindex(u))
  nl_rhs = alg.nlsolve(Val{:init},rhs,uhold)
  IIF2Cache(u,uprev,uhold,dual_cache,tmp,rhs,nl_rhs,rtmp1,fsalfirst,k)
end

immutable LawsonEulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  rtmp::rateType
  fsalfirst::rateType
end

function alg_cache(alg::LawsonEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  LawsonEulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype),zeros(rate_prototype))
end

u_cache(c::LawsonEulerCache) = ()
du_cache(c::LawsonEulerCache) = (c.k,c.fsalfirst,c.rtmp)

immutable LawsonEulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::LawsonEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = LawsonEulerConstantCache()

immutable NorsettEulerCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  rtmp::rateType
  fsalfirst::rateType
end

function alg_cache(alg::NorsettEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  NorsettEulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype),zeros(rate_prototype))
end

u_cache(c::NorsettEulerCache) = ()
du_cache(c::NorsettEulerCache) = (c.k,c.fsalfirst)

immutable NorsettEulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::NorsettEuler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = NorsettEulerConstantCache()

get_chunksize(cache::DECache) = error("This cache does not have a chunksize.")
get_chunksize{uType,DiffCacheType,rateType,CS}(cache::ImplicitEulerCache{uType,DiffCacheType,rateType,CS}) = CS
get_chunksize{uType,DiffCacheType,rateType,CS}(cache::TrapezoidCache{uType,DiffCacheType,rateType,CS}) = CS
get_chunksize{CS}(cache::ODEChunkCache{CS}) = CS
