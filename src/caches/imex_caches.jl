mutable struct CNABConstantCache{uType,tType,F,uEltypeNoUnits,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  tab::Tab
  uprev3::uType
  tprev2::tType
end

function alg_cache(alg::CNAB,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uBottomEltypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  if typeof(f) <: SplitFunction
    uf = DiffEqDiffTools.UDerivativeWrapper(f.f1,t,p)
  else
    uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  end
  ηold = one(uEltypeNoUnits)

  uprev3 = u
  tprev2 = t

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  #tab = CNABTableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  CNABConstantCache(uf,ηold,κ,tol,10000)
end

mutable struct CNABCache{uType,rateType,uNoUnitsType,J,UF,JC,uEltypeNoUnits,Tab,F,kType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  uprev3::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  y1::kType
  y0::kType
  dz::uType
  b::uType
  tmp::uType
  atmp::uNoUnitsType
  J::J
  W::J
  uf::UF
  jac_config::JC
  linsolve::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  tab::Tab
end

u_cache(c::CNABCache)    = (c.z,c.dz)
du_cache(c::CNABCache)   = (c.fsalfirst)

function alg_cache(alg::CNAB,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z = similar(u,indices(u))
  dz = similar(u,indices(u))
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)
  tmp = similar(u); b = similar(u,indices(u));
  atmp = similar(u,uEltypeNoUnits,indices(u))

  if typeof(f) <: SplitFunction
    y0 = similar(u,indices(u)); y1 = similar(u,indices(u))
    uf = DiffEqDiffTools.UJacobianWrapper(f.f1,t,p)
  else
    y0 = nothing; y1 = nothing
    uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  end
  linsolve = alg.linsolve(Val{:init},uf,u)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  #tab = KenCarp3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  ηold = one(uEltypeNoUnits)

  uprev3 = similar(u)
  tprev2 = t

  CNABCache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uEltypeNoUnits,typeof(tab),typeof(linsolve),typeof(k1)}(
              u,uprev,du1,fsalfirst,k,z,y1,y0,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000)
end