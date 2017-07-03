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

immutable Symplectic2Cache{uType,rateType,tableauType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tableauType
end

u_cache(c::Symplectic2Cache) = ()
du_cache(c::Symplectic2Cache) = (c.k,c.fsalfirst)

function alg_cache(alg::VerletLeapfrog,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = VerletLeapfrogConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic2Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::VerletLeapfrog,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = VerletLeapfrogConstantCache()

function alg_cache(alg::PseudoVerletLeapfrog,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = PseudoVerletLeapfrogConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic2Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::PseudoVerletLeapfrog,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = PseudoVerletLeapfrogConstantCache()

function alg_cache(alg::McAte2,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = McAte2ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic2Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::McAte2,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = McAte2ConstantCache()

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

alg_cache(alg::Ruth3,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = Ruth3ConstantCache()

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

immutable Symplectic45Cache{uType,rateType,tableauType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tableauType
end

u_cache(c::Symplectic45Cache) = ()
du_cache(c::Symplectic45Cache) = (c.k,c.fsalfirst)

function alg_cache(alg::CalvoSanz4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = CalvoSanz4ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic45Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::CalvoSanz4,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = CalvoSanz4ConstantCache()

function alg_cache(alg::McAte42,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = McAte42ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic45Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::McAte42,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = McAte42ConstantCache()

immutable Symplectic5Cache{uType,rateType,tableauType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tableauType
end

u_cache(c::Symplectic5Cache) = ()
du_cache(c::Symplectic5Cache) = (c.k,c.fsalfirst)

function alg_cache(alg::McAte5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = McAte5ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic5Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::McAte5,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = McAte5ConstantCache()

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

immutable Symplectic62Cache{uType,rateType,tableauType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tableauType
end

u_cache(c::Symplectic62Cache) = ()
du_cache(c::Symplectic62Cache) = (c.k,c.fsalfirst)

function alg_cache(alg::KahanLi6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = KahanLi6ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  Symplectic62Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::KahanLi6,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = KahanLi6ConstantCache()

immutable McAte8Cache{uType,rateType,tableauType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tableauType
end

u_cache(c::McAte8Cache) = ()
du_cache(c::McAte8Cache) = (c.k,c.fsalfirst)

function alg_cache(alg::McAte8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = McAte8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  McAte8Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::McAte8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = McAte8ConstantCache()

immutable KahanLi8Cache{uType,rateType,tableauType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tableauType
end

u_cache(c::KahanLi8Cache) = ()
du_cache(c::KahanLi8Cache) = (c.k,c.fsalfirst)

function alg_cache(alg::KahanLi8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = KahanLi8ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  KahanLi8Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::KahanLi8,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = KahanLi8ConstantCache()

immutable SofSpa10Cache{uType,rateType,tableauType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tableauType
end

u_cache(c::SofSpa10Cache) = ()
du_cache(c::SofSpa10Cache) = (c.k,c.fsalfirst)

function alg_cache(alg::SofSpa10,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  tab = SofSpa10ConstantCache(realtype(uEltypeNoUnits),realtype(tTypeNoUnits))
  SofSpa10Cache(u,uprev,k,tmp,fsalfirst,tab)
end

alg_cache(alg::SofSpa10,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,::Type{Val{false}}) = SofSpa10ConstantCache()
