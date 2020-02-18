@cache mutable struct DImplicitEulerCache{uType,duType,uNoUnitsType,N} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  du::duType
  atmp::uNoUnitsType
  nlsolver::N
end

mutable struct DImplicitEulerConstantCache{duType,N} <: OrdinaryDiffEqConstantCache
  du::duType
  nlsolver::N
end

function alg_cache(alg::DImplicitEuler,du,u,res_prototype,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  γ, c = 1, 1
  α = 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,res_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,α,Val(false))

  DImplicitEulerConstantCache(du,nlsolver)
end

function alg_cache(alg::DImplicitEuler,du,u,res_prototype,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})

  γ, c = 1, 1
  α = 1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,res_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,α,Val(true))

  atmp = similar(u,uEltypeNoUnits)

  DImplicitEulerCache(u,uprev,uprev2,du,atmp,nlsolver)
end

@cache mutable struct DABDF2ConstantCache{duType,N,dtType,rate_prototype} <: OrdinaryDiffEqConstantCache
  du::duType
  nlsolver::N
  eulercache::DImplicitEulerConstantCache
  dtₙ₋₁::dtType
  fsalfirstprev::rate_prototype
end

function alg_cache(alg::DABDF2,du,u,res_prototype,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  γ, c = 1//1, 1
  α = 1//1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,res_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,α,Val(false))
  eulercache = DImplicitEulerConstantCache(du,nlsolver)

  dtₙ₋₁ = one(dt)
  fsalfirstprev = rate_prototype

  DABDF2ConstantCache(du,nlsolver, eulercache, dtₙ₋₁, fsalfirstprev)
end

@cache mutable struct DABDF2Cache{uType,duType,rateType,uNoUnitsType,N,dtType} <: OrdinaryDiffEqMutableCache
  uₙ::uType
  uₙ₋₁::uType
  uₙ₋₂::uType
  du::duType
  fsalfirst::rateType
  fsalfirstprev::rateType
  atmp::uNoUnitsType
  nlsolver::N
  eulercache::DImplicitEulerCache
  dtₙ₋₁::dtType
end

function alg_cache(alg::DABDF2,du,u,res_prototype,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true})
  γ, c = 1//1, 1
  α = 1//1
  nlsolver = build_nlsolver(alg,u,uprev,p,t,dt,f,res_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,α,Val(true))
  fsalfirst = zero(rate_prototype)

  fsalfirstprev = zero(rate_prototype)
  atmp = similar(u,uEltypeNoUnits)

  eulercache = DImplicitEulerCache(u,uprev,uprev2,du,atmp,nlsolver)

  dtₙ₋₁ = one(dt)

  DABDF2Cache(u,uprev,uprev2,du,fsalfirst,fsalfirstprev,atmp,
              nlsolver,eulercache,dtₙ₋₁)
end
