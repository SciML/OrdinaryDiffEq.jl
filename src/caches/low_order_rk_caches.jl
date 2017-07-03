function alg_cache(alg::Euler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  EulerCache(u,uprev,similar(u),zeros(rate_prototype),zeros(rate_prototype))
end

u_cache(c::EulerCache) = ()
du_cache(c::EulerCache) = (c.k,c.fsalfirst)

immutable EulerConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Euler,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = EulerConstantCache()

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
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
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
