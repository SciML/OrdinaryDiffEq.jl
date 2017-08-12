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

function alg_cache(alg::Nystrom4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  k₄ = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  Nystrom4Cache(u,uprev,k₁,k₂,k₃,k₄,k,tmp)
end

# alg_cache(alg::Nystrom4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{false}}) = Nystrom4ConstantCache()

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

function alg_cache(alg::Nystrom4VelocityIndependent,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  Nystrom4VelocityIndependentCache(u,uprev,k₁,k₂,k₃,k,tmp)
end

struct IRKN4Cache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  k₂::rateType
  k₃::rateType
  k::rateType
  tmp::uType
  onestep_cache::Nystrom4VelocityIndependentCache
end

u_cache(c::IRKN4Cache) = ()
du_cache(c::IRKN4Cache) = (c.fsalfirst,c.k₂,c.k₃,c.k,c.k_₁,c.k_₂,c.k_₃)

function alg_cache(alg::IRKN4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  IRKN4Cache(u,uprev,uprev2,k₁,k₂,k₃,k,tmp,Nystrom4VelocityIndependentCache(u,uprev,k₁,k₂,k₃,k,tmp))
end

struct Nystrom5VelocityIndependentCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k₂::rateType
  k₃::rateType
  k₄::rateType
  k::rateType
  tmp::uType
end

u_cache(c::Nystrom5VelocityIndependentCache) = ()
du_cache(c::Nystrom5VelocityIndependentCache) = (c.fsalfirst,c.k₂,c.k₃,c.k₄,c.k)

function alg_cache(alg::Nystrom5VelocityIndependent,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})
  k₁ = zeros(rate_prototype)
  k₂ = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  k₄ = zeros(rate_prototype)
  k  = zeros(rate_prototype)
  tmp = similar(u)
  Nystrom5VelocityIndependentCache(u,uprev,k₁,k₂,k₃,k₄,k,tmp)
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
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DPRKN6Cache) = (c.atmp,c.utilde)
du_cache(c::DPRKN6Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6)

function alg_cache(alg::DPRKN6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})
  tab = DPRKN6ConstantCache(real(uEltypeNoUnits),real(tTypeNoUnits))
  k1 = zeros(rate_prototype)
  k2 = zeros(rate_prototype)
  k3 = zeros(rate_prototype)
  k4 = zeros(rate_prototype)
  k5 = zeros(rate_prototype)
  k6 = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  DPRKN6Cache(u,uprev,k1,k2,k3,k4,k5,k6,utilde,tmp,atmp,tab)
end

