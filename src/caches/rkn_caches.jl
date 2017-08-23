struct Nystrom4Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k₂::rateType
  k₃::rateType
  k₄::rateType
  k::rateType
  tmp::uType
end

u_cache(c::Nystrom4Cache) = ()
du_cache(c::Nystrom4Cache) = (c.fsalfirst,c.k₂,c.k₃,c.k₄,c.k)

# struct Nystrom4ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::Nystrom4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  k₄ = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  Nystrom4Cache(u,uprev,k₁,k₂,k₃,k₄,k,tmp)
end

# alg_cache(alg::Nystrom4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = Nystrom4ConstantCache()

struct Nystrom4VelocityIndependentCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k₂::rateType
  k₃::rateType
  k::rateType
  tmp::uType
end

u_cache(c::Nystrom4VelocityIndependentCache) = ()
du_cache(c::Nystrom4VelocityIndependentCache) = (c.fsalfirst,c.k₂,c.k₃,c.k)

function alg_cache(alg::Nystrom4VelocityIndependent,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  Nystrom4VelocityIndependentCache(u,uprev,k₁,k₂,k₃,k,tmp)
end

struct IRKN3Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  k₂::rateType
  k::rateType
  tmp::uType
  onestep_cache::Nystrom4VelocityIndependentCache
  tab::TabType
end

u_cache(c::IRKN3Cache) = ()
du_cache(c::IRKN3Cache) = (c.fsalfirst,c.k₂,c.k)

function alg_cache(alg::IRKN3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  tab = IRKN3ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))
  IRKN3Cache(u,uprev,uprev2,k₁,k₂,k,tmp,Nystrom4VelocityIndependentCache(u,uprev,k₁,k₂,k₃,k,tmp),tab)
end

struct IRKN4Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  k₂::rateType
  k₃::rateType
  k::rateType
  tmp::uType
  tmp2::uType
  onestep_cache::Nystrom4VelocityIndependentCache
  tab::TabType
end

u_cache(c::IRKN4Cache) = ()
du_cache(c::IRKN4Cache) = (c.fsalfirst,c.k₂,c.k₃,c.k)

function alg_cache(alg::IRKN4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  tmp2 = similar(u)
  tab = IRKN4ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))
  IRKN4Cache(u,uprev,uprev2,k₁,k₂,k₃,k,tmp,tmp2,Nystrom4VelocityIndependentCache(u,uprev,k₁,k₂,k₃,k,tmp),tab)
end

struct Nystrom5VelocityIndependentCache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k₂::rateType
  k₃::rateType
  k₄::rateType
  k::rateType
  tmp::uType
  tab::TabType
end

u_cache(c::Nystrom5VelocityIndependentCache) = ()
du_cache(c::Nystrom5VelocityIndependentCache) = (c.fsalfirst,c.k₂,c.k₃,c.k₄,c.k)

function alg_cache(alg::Nystrom5VelocityIndependent,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  k₄ = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  tab = Nystrom5VelocityIndependentConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))
  Nystrom5VelocityIndependentCache(u,uprev,k₁,k₂,k₃,k₄,k,tmp,tab)
end

struct DPRKN6Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::rateType
  k3::rateType
  k4::rateType
  k5::rateType
  k6::rateType
  k::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DPRKN6Cache) = (c.atmp,c.utilde)
du_cache(c::DPRKN6Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6)

function alg_cache(alg::DPRKN6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tab = DPRKN6ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  DPRKN6Cache(u,uprev,k1,k2,k3,k4,k5,k6,k,utilde,tmp,atmp,tab)
end

struct DPRKN8Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  k::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DPRKN8Cache) = (c.atmp,c.utilde)
du_cache(c::DPRKN8Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9)

function alg_cache(alg::DPRKN8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tab = DPRKN8ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype)
  k7 = zeros(rate_prototype)
  k8 = zeros(rate_prototype)
  k9 = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  DPRKN8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k,utilde,tmp,atmp,tab)
end

struct DPRKN12Cache{uType,uArrayType,rateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
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
  k::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DPRKN12Cache) = (c.atmp,c.utilde)
du_cache(c::DPRKN12Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.k17)

function alg_cache(alg::DPRKN12,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tab = DPRKN12ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype)
  k7 = zeros(rate_prototype)
  k8 = zeros(rate_prototype)
  k9 = zeros(rate_prototype)
  k10 = zeros(rate_prototype)
  k11 = zeros(rate_prototype)
  k12 = zeros(rate_prototype)
  k13 = zeros(rate_prototype)
  k14 = zeros(rate_prototype)
  k15 = zeros(rate_prototype)
  k16 = zeros(rate_prototype)
  k17 = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  DPRKN12Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k,utilde,tmp,atmp,tab)
end
