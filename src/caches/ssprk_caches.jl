struct SSPRK22Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

u_cache(c::SSPRK22Cache) = ()
du_cache(c::SSPRK22Cache) = (c.k,c.fsalfirst)

struct SSPRK22ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK22,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  SSPRK22Cache(u,uprev,k,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK22,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = SSPRK22ConstantCache()


struct SSPRK33Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

u_cache(c::SSPRK33Cache) = ()
du_cache(c::SSPRK33Cache) = (c.k,c.fsalfirst)

struct SSPRK33ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK33,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  SSPRK33Cache(u,uprev,k,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK33,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = SSPRK33ConstantCache()


struct SSPRK432Cache{uType,rateType,uArrayType,uEltypeNoUnits,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  utilde::uArrayType
  atmp::uEltypeNoUnits
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

u_cache(c::SSPRK432Cache) = (c.utilde,c.atmp)
du_cache(c::SSPRK432Cache) = (c.k,c.fsalfirst)

struct SSPRK432ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK432,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u))
  SSPRK432Cache(u,uprev,k,tmp,fsalfirst,utilde,atmp,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK432,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = SSPRK432ConstantCache()


struct SSPRK932Cache{uType,rateType,uArrayType,uEltypeNoUnits,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  utilde::uArrayType
  atmp::uEltypeNoUnits
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

u_cache(c::SSPRK932Cache) = (c.utilde,c.atmp)
du_cache(c::SSPRK932Cache) = (c.k,c.fsalfirst)

struct SSPRK932ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK932,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  utilde = similar(u,indices(u))
  atmp = similar(u,uEltypeNoUnits,indices(u))
  SSPRK932Cache(u,uprev,k,tmp,fsalfirst,utilde,atmp,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK932,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = SSPRK932ConstantCache()


struct SSPRK104Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₄::rateType
  tmp::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

u_cache(c::SSPRK104Cache) = ()
du_cache(c::SSPRK104Cache) = (c.k,c.fsalfirst,c.k₄)

struct SSPRK104ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK104,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tmp = similar(u)
  k = zeros(rate_prototype)
  k₄ = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  SSPRK104Cache(u,uprev,k,k₄,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK104,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}}) = SSPRK104ConstantCache()
