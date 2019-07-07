@cache struct PDIRK44Cache{uType,rateType,N,JType,WType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k1::Array{rateType}
  k2::Array{rateType}
  nlsolver::N
  J::Array{JType}
  W::Array{WType}
end

struct PDIRK44ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::PDIRK44,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  γ, c = 1.0, 1.0
  J1, W1 = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver1 = iipnlsolve(alg,u,uprev,p,t,dt,f,W1,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  J2, W2 = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver2 = iipnlsolve(alg,u,uprev,p,t,dt,f,W2,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  nlsolver = [nlsolver1, nlsolver2]
  W = [W1, W2]
  J = [J1, J2]
  k1 = [zero(rate_prototype) for i in 1:2 ]
  k2 = [zero(rate_prototype) for i in 1:2 ]
  PDIRK44Cache(u,uprev,k1,k2,nlsolver,J,W)
end

function alg_cache(alg::PDIRK44,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
end
