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

immutable VelocityVerletCache{uType,rateType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
end

u_cache(c::VelocityVerletCache) = ()
du_cache(c::VelocityVerletCache) = (c.k,c.fsalfirst)

immutable VelocityVerletConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::VelocityVerlet,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = zeros(rate_prototype)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  VelocityVerletCache(u,uprev,k,tmp,fsalfirst)
end

alg_cache(alg::VelocityVerlet,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = VelocityVerletConstantCache()

immutable Symplectic3Cache{uType,rateType,tableauType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tableauType
end

u_cache(c::Symplectic3Cache) = ()
du_cache(c::Symplectic3Cache) = (c.k,c.fsalfirst)

function alg_cache(alg::Ruth3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = Ruth3ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic3Cache(u,uprev,k,tmp,fsalfirst,tab)
end

function alg_cache(alg::McAte3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = McAte3ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic3Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::McAte3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = McAte3ConstantCache()

immutable Symplectic4Cache{uType,rateType,tableauType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tableauType
end

u_cache(c::Symplectic4Cache) = ()
du_cache(c::Symplectic4Cache) = (c.k,c.fsalfirst)

function alg_cache(alg::McAte4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = McAte4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic4Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::McAte4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = McAte4ConstantCache()

function alg_cache(alg::CandyRoz4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = CandyRoz4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic4Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::CandyRoz4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = McAte4ConstantCache()

immutable CalvoSanz4Cache{uType,rateType,tableauType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tableauType
end

u_cache(c::CalvoSanz4) = ()
du_cache(c::CalvoSanz4) = (c.k,c.fsalfirst)

function alg_cache(alg::CalvoSanz4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = CalvoSanz4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  CalvoSanz4Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::CalvoSanz4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = CalvoSanz4ConstantCache()

immutable Symplectic6Cache{uType,rateType,tableauType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tableauType
end

u_cache(c::Symplectic6Cache) = ()
du_cache(c::Symplectic6Cache) = (c.k,c.fsalfirst)

function alg_cache(alg::Yoshida6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = Yoshida6ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic6Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::Yoshida6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Yoshida6ConstantCache()
