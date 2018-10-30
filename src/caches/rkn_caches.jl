struct Nystrom4Cache{uType,rateType,reducedRateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k₂::reducedRateType
  k₃::reducedRateType
  k₄::reducedRateType
  k::rateType
  tmp::uType
end

u_cache(c::Nystrom4Cache) = ()
du_cache(c::Nystrom4Cache) = (c.fsalfirst,c.k₂,c.k₃,c.k₄,c.k)

# struct Nystrom4ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::Nystrom4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  reduced_rate_prototype = rate_prototype.x[2]
  k₁ = zero(rate_prototype)
  k₂ = zero(reduced_rate_prototype)
  k₃ = zero(reduced_rate_prototype)
  k₄ = zero(reduced_rate_prototype)
  k  = zero(rate_prototype)
  tmp = similar(u)
  Nystrom4Cache(u,uprev,k₁,k₂,k₃,k₄,k,tmp)
end

struct Nystrom4ConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Nystrom4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = Nystrom4ConstantCache()

# alg_cache(alg::Nystrom4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = Nystrom4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

struct Nystrom4VelocityIndependentCache{uType,rateType,reducedRateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k₂::reducedRateType
  k₃::reducedRateType
  k::rateType
  tmp::uType
end

u_cache(c::Nystrom4VelocityIndependentCache) = ()
du_cache(c::Nystrom4VelocityIndependentCache) = (c.fsalfirst,c.k₂,c.k₃,c.k)

function alg_cache(alg::Nystrom4VelocityIndependent,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  reduced_rate_prototype = rate_prototype.x[2]
  k₁ = zero(rate_prototype)
  k₂ = zero(reduced_rate_prototype)
  k₃ = zero(reduced_rate_prototype)
  k  = zero(rate_prototype)
  tmp = similar(u)
  Nystrom4VelocityIndependentCache(u,uprev,k₁,k₂,k₃,k,tmp)
end

struct Nystrom4VelocityIndependentConstantCache <: OrdinaryDiffEqConstantCache end

alg_cache(alg::Nystrom4VelocityIndependent,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = Nystrom4VelocityIndependentConstantCache()

struct IRKN3Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  k₂::rateType
  k::rateType
  tmp::uType
  tmp2::rateType
  onestep_cache::Nystrom4VelocityIndependentCache
  tab::TabType
end

u_cache(c::IRKN3Cache) = ()
du_cache(c::IRKN3Cache) = (c.fsalfirst,c.k₂,c.k)

function alg_cache(alg::IRKN3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  k₁ = zero(rate_prototype)
  k₂ = zero(rate_prototype)
  k₃ = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = similar(u)
  tab = IRKN3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  IRKN3Cache(u,uprev,uprev2,k₁,k₂,k,tmp,k₃,Nystrom4VelocityIndependentCache(u,uprev,k₁,k₂.x[2],k₃.x[2],k,tmp),tab)
end

alg_cache(alg::IRKN3,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = IRKN3ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

struct IRKN4Cache{uType,rateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  fsalfirst::rateType
  k₂::rateType
  k₃::rateType
  k::rateType
  tmp::uType
  tmp2::rateType
  onestep_cache::Nystrom4VelocityIndependentCache
  tab::TabType
end

u_cache(c::IRKN4Cache) = ()
du_cache(c::IRKN4Cache) = (c.fsalfirst,c.k₂,c.k₃,c.k)

function alg_cache(alg::IRKN4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  k₁ = zero(rate_prototype)
  k₂ = zero(rate_prototype)
  k₃ = zero(rate_prototype)
  k  = zero(rate_prototype)
  tmp = similar(u)
  tmp2 = similar(rate_prototype)
  tab = IRKN4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  IRKN4Cache(u,uprev,uprev2,k₁,k₂,k₃,k,tmp,tmp2,Nystrom4VelocityIndependentCache(u,uprev,k₁,k₂.x[2],k₃.x[2],k,tmp),tab)
end

alg_cache(alg::IRKN4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = IRKN4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

struct Nystrom5VelocityIndependentCache{uType,rateType,reducedRateType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k₂::reducedRateType
  k₃::reducedRateType
  k₄::reducedRateType
  k::rateType
  tmp::uType
  tab::TabType
end

u_cache(c::Nystrom5VelocityIndependentCache) = ()
du_cache(c::Nystrom5VelocityIndependentCache) = (c.fsalfirst,c.k₂,c.k₃,c.k₄,c.k)

function alg_cache(alg::Nystrom5VelocityIndependent,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  reduced_rate_prototype = rate_prototype.x[2]
  k₁ = zero(rate_prototype)
  k₂ = zero(reduced_rate_prototype)
  k₃ = zero(reduced_rate_prototype)
  k₄ = zero(reduced_rate_prototype)
  k  = zero(rate_prototype)
  tmp = similar(u)
  tab = Nystrom5VelocityIndependentConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  Nystrom5VelocityIndependentCache(u,uprev,k₁,k₂,k₃,k₄,k,tmp,tab)
end

alg_cache(alg::Nystrom5VelocityIndependent,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = Nystrom5VelocityIndependentConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

struct DPRKN6Cache{uType,uArrayType,rateType,reducedRateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::reducedRateType
  k3::reducedRateType
  k4::reducedRateType
  k5::reducedRateType
  k6::reducedRateType
  k::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DPRKN6Cache) = (c.atmp,c.utilde)
du_cache(c::DPRKN6Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6)

function alg_cache(alg::DPRKN6,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  reduced_rate_prototype = rate_prototype.x[2]
  tab = DPRKN6ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(reduced_rate_prototype)
  k3 = zero(reduced_rate_prototype)
  k4 = zero(reduced_rate_prototype)
  k5 = zero(reduced_rate_prototype)
  k6 = zero(reduced_rate_prototype)
  k  = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  DPRKN6Cache(u,uprev,k1,k2,k3,k4,k5,k6,k,utilde,tmp,atmp,tab)
end

alg_cache(alg::DPRKN6,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = DPRKN6ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))


struct DPRKN8Cache{uType,uArrayType,rateType,reducedRateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::reducedRateType
  k3::reducedRateType
  k4::reducedRateType
  k5::reducedRateType
  k6::reducedRateType
  k7::reducedRateType
  k8::reducedRateType
  k9::reducedRateType
  k::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DPRKN8Cache) = (c.atmp,c.utilde)
du_cache(c::DPRKN8Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9)

function alg_cache(alg::DPRKN8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  reduced_rate_prototype = rate_prototype.x[2]
  tab = DPRKN8ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(reduced_rate_prototype)
  k3 = zero(reduced_rate_prototype)
  k4 = zero(reduced_rate_prototype)
  k5 = zero(reduced_rate_prototype)
  k6 = zero(reduced_rate_prototype)
  k7 = zero(reduced_rate_prototype)
  k8 = zero(reduced_rate_prototype)
  k9 = zero(reduced_rate_prototype)
  k  = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  DPRKN8Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k,utilde,tmp,atmp,tab)
end

alg_cache(alg::DPRKN8,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = DPRKN8ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

struct DPRKN12Cache{uType,uArrayType,rateType,reducedRateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::reducedRateType
  k3::reducedRateType
  k4::reducedRateType
  k5::reducedRateType
  k6::reducedRateType
  k7::reducedRateType
  k8::reducedRateType
  k9::reducedRateType
  k10::reducedRateType
  k11::reducedRateType
  k12::reducedRateType
  k13::reducedRateType
  k14::reducedRateType
  k15::reducedRateType
  k16::reducedRateType
  k17::reducedRateType
  k::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::DPRKN12Cache) = (c.atmp,c.utilde)
du_cache(c::DPRKN12Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k5,c.k6,c.k7,c.k8,c.k9,c.k10,c.k11,c.k12,c.k13,c.k14,c.k15,c.k16,c.k17)

function alg_cache(alg::DPRKN12,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  reduced_rate_prototype = rate_prototype.x[2]
  tab = DPRKN12ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(reduced_rate_prototype)
  k3 = zero(reduced_rate_prototype)
  k4 = zero(reduced_rate_prototype)
  k5 = zero(reduced_rate_prototype)
  k6 = zero(reduced_rate_prototype)
  k7 = zero(reduced_rate_prototype)
  k8 = zero(reduced_rate_prototype)
  k9 = zero(reduced_rate_prototype)
  k10 = zero(reduced_rate_prototype)
  k11 = zero(reduced_rate_prototype)
  k12 = zero(reduced_rate_prototype)
  k13 = zero(reduced_rate_prototype)
  k14 = zero(reduced_rate_prototype)
  k15 = zero(reduced_rate_prototype)
  k16 = zero(reduced_rate_prototype)
  k17 = zero(reduced_rate_prototype)
  k  = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  DPRKN12Cache(u,uprev,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k,utilde,tmp,atmp,tab)
end

alg_cache(alg::DPRKN12,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = DPRKN12ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

struct ERKN4Cache{uType,uArrayType,rateType,reducedRateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::reducedRateType
  k3::reducedRateType
  k4::reducedRateType
  k::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::ERKN4Cache) = (c.atmp,c.utilde)
du_cache(c::ERKN4Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k)

function alg_cache(alg::ERKN4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  reduced_rate_prototype = rate_prototype.x[2]
  tab = ERKN4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(reduced_rate_prototype)
  k3 = zero(reduced_rate_prototype)
  k4 = zero(reduced_rate_prototype)
  k  = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  ERKN4Cache(u,uprev,k1,k2,k3,k4,k,utilde,tmp,atmp,tab)
end

alg_cache(alg::ERKN4,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = ERKN4ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

struct ERKN5Cache{uType,uArrayType,rateType,reducedRateType,uEltypeNoUnits,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  k2::reducedRateType
  k3::reducedRateType
  k4::reducedRateType
  k::rateType
  utilde::uArrayType
  tmp::uType
  atmp::uEltypeNoUnits
  tab::TabType
end

u_cache(c::ERKN5Cache) = (c.atmp,c.utilde)
du_cache(c::ERKN5Cache) = (c.fsalfirst,c.k2,c.k3,c.k4,c.k)

function alg_cache(alg::ERKN5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  reduced_rate_prototype = rate_prototype.x[2]
  tab = ERKN5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  k1 = zero(rate_prototype)
  k2 = zero(reduced_rate_prototype)
  k3 = zero(reduced_rate_prototype)
  k4 = zero(reduced_rate_prototype)
  k  = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  tmp = similar(u)
  ERKN5Cache(u,uprev,k1,k2,k3,k4,k,utilde,tmp,atmp,tab)
end

alg_cache(alg::ERKN5,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = ERKN5ConstantCache(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
