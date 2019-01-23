@cache struct RichardsonEulerCache{uType,rateType,tabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  tab::tabType
  T::Array{uType,2}
end

struct RichardsonEulerConstantCache <: OrdinaryDiffEqConstantCache
  m::Int #Generalise type of m
  # Can add different order of n_{i} here
  function RichardsonEulerConstantCache()
    m = 2
    new(m)
  end
end

function alg_cache(alg::RichardsonEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = RichardsonEulerConstantCache()
  T = zeros(typeof(u), (tab.m, tab.m))
  RichardsonEulerCache(u,uprev,tmp,k,fsalfirst,tab,T)
end

function alg_cache(alg::RichardsonEuler,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  RichardsonEulerConstantCache()
end
