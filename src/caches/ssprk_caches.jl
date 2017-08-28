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


struct SSPRK53Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  u₂::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  α30::Float64
  α32::Float64
  α40::Float64
  α43::Float64
  α52::Float64
  α54::Float64
  β10::Float64
  β21::Float64
  β32::Float64
  β43::Float64
  β54::Float64
  c1::Float64
  c2::Float64
  c3::Float64
  c4::Float64
end

u_cache(c::SSPRK53Cache) = (c.tmp,c.u₂)
du_cache(c::SSPRK53Cache) = (c.k,c.fsalfirst)

struct SSPRK53ConstantCache <: OrdinaryDiffEqConstantCache
  α30::Float64
  α32::Float64
  α40::Float64
  α43::Float64
  α52::Float64
  α54::Float64
  β10::Float64
  β21::Float64
  β32::Float64
  β43::Float64
  β54::Float64
  c1::Float64
  c2::Float64
  c3::Float64
  c4::Float64
end

function alg_cache(alg::SSPRK53,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tmp = similar(u)
  u₂ = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  α30 = 0.355909775063327
  α32 = 0.644090224936674
  α40 = 0.367933791638137
  α43 = 0.632066208361863
  α52 = 0.237593836598569
  α54 = 0.762406163401431
  β10 = 0.377268915331368
  β21 = 0.377268915331368
  β32 = 0.242995220537396
  β43 = 0.238458932846290
  β54 = 0.287632146308408
  c1 = 0.377268915331368
  c2 = 0.754537830662736
  c3 = 0.728985661612188
  c4 = 0.699226135931670
  SSPRK53Cache(u,uprev,k,tmp,u₂,fsalfirst,alg.stage_limiter!,alg.step_limiter!,
                α30,α32,α40,α43,α52,α54,β10,β21,β32,β43,β54,c1,c2,c3,c4)
end

function alg_cache(alg::SSPRK53,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}})
  α30 = 0.355909775063327
  α32 = 0.644090224936674
  α40 = 0.367933791638137
  α43 = 0.632066208361863
  α52 = 0.237593836598569
  α54 = 0.762406163401431
  β10 = 0.377268915331368
  β21 = 0.377268915331368
  β32 = 0.242995220537396
  β43 = 0.238458932846290
  β54 = 0.287632146308408
  c1 = 0.377268915331368
  c2 = 0.754537830662736
  c3 = 0.728985661612188
  c4 = 0.699226135931670
  SSPRK53ConstantCache(α30,α32,α40,α43,α52,α54,β10,β21,β32,β43,β54,c1,c2,c3,c4)
end


struct SSPRK63Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  u₂::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  α40::Float64
  α41::Float64
  α43::Float64
  α62::Float64
  α65::Float64
  β10::Float64
  β21::Float64
  β32::Float64
  β43::Float64
  β54::Float64
  β65::Float64
  c1::Float64
  c2::Float64
  c3::Float64
  c4::Float64
  c5::Float64
end

u_cache(c::SSPRK63Cache) = (c.tmp,c.u₂)
du_cache(c::SSPRK63Cache) = (c.k,c.fsalfirst)

struct SSPRK63ConstantCache <: OrdinaryDiffEqConstantCache
  α40::Float64
  α41::Float64
  α43::Float64
  α62::Float64
  α65::Float64
  β10::Float64
  β21::Float64
  β32::Float64
  β43::Float64
  β54::Float64
  β65::Float64
  c1::Float64
  c2::Float64
  c3::Float64
  c4::Float64
  c5::Float64
end

function alg_cache(alg::SSPRK63,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  tmp = similar(u)
  u₂ = similar(u)
  k = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  α40 = 0.476769811285196
  α41 = 0.098511733286064
  α43 = 0.424718455428740
  α62 = 0.155221702560091
  α65 = 0.844778297439909
  β10 = 0.284220721334261
  β21 = 0.284220721334261
  β32 = 0.284220721334261
  β43 = 0.120713785765930
  β54 = 0.284220721334261
  β65 = 0.240103497065900
  c1 = 0.284220721334261
  c2 = 0.568441442668522
  c3 = 0.852662164002783
  c4 = 0.510854218958172
  c5 = 0.795074940292433
  SSPRK63Cache(u,uprev,k,tmp,u₂,fsalfirst,alg.stage_limiter!,alg.step_limiter!,
                α40,α41,α43,α62,α65,β10,β21,β32,β43,β54,β65,c1,c2,c3,c4,c5)
end

function alg_cache(alg::SSPRK63,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}})
  α40 = 0.476769811285196
  α41 = 0.098511733286064
  α43 = 0.424718455428740
  α62 = 0.155221702560091
  α65 = 0.844778297439909
  β10 = 0.284220721334261
  β21 = 0.284220721334261
  β32 = 0.284220721334261
  β43 = 0.120713785765930
  β54 = 0.284220721334261
  β65 = 0.240103497065900
  c1 = 0.284220721334261
  c2 = 0.568441442668522
  c3 = 0.852662164002783
  c4 = 0.510854218958172
  c5 = 0.795074940292433
  SSPRK63ConstantCache(α40,α41,α43,α62,α65,β10,β21,β32,β43,β54,β65,c1,c2,c3,c4,c5)
end


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


struct SSPRK54Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₃::rateType
  u₂::uType
  u₃::uType
  u₄::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  β10::Float64
  α20::Float64
  α21::Float64
  β21::Float64
  α30::Float64
  α32::Float64
  β32::Float64
  α40::Float64
  α43::Float64
  β43::Float64
  α52::Float64
  α53::Float64
  β53::Float64
  α54::Float64
  β54::Float64
  c1::Float64
  c2::Float64
  c3::Float64
  c4::Float64
end

u_cache(c::SSPRK54Cache) = (c.u₂,c.u₃,c.u₄)
du_cache(c::SSPRK54Cache) = (c.k,c.fsalfirst,c.k₃)

struct SSPRK54ConstantCache <: OrdinaryDiffEqConstantCache
  β10::Float64
  α20::Float64
  α21::Float64
  β21::Float64
  α30::Float64
  α32::Float64
  β32::Float64
  α40::Float64
  α43::Float64
  β43::Float64
  α52::Float64
  α53::Float64
  β53::Float64
  α54::Float64
  β54::Float64
  c1::Float64
  c2::Float64
  c3::Float64
  c4::Float64
end

function alg_cache(alg::SSPRK54,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})
  u₂ = similar(u)
  u₃ = similar(u)
  u₄ = similar(u)
  k = zeros(rate_prototype)
  k₃ = zeros(rate_prototype)
  fsalfirst = zeros(rate_prototype)
  β10 = 0.391752226571890
  α20 = 0.444370493651235
  α21 = 0.555629506348765
  β21 = 0.368410593050371
  α30 = 0.620101851488403
  α32 = 0.379898148511597
  β32 = 0.251891774271694
  α40 = 0.178079954393132
  α43 = 0.821920045606868
  β43 = 0.544974750228521
  α52 = 0.517231671970585
  α53 = 0.096059710526147
  β53 = 0.063692468666290
  α54 = 0.386708617503269
  β54 = 0.226007483236906
  c1 = 0.39175222700391998
  c2 = 0.58607968896779994
  c3 = 0.47454236302687003
  c4 = 0.93501063100923998
  SSPRK54Cache(u,uprev,k,k₃,u₂,u₃,u₄,fsalfirst,alg.stage_limiter!,alg.step_limiter!,
                β10,α20,α21,β21,α30,α32,β32,α40,α43,β43,α52,α53,β53,α54,β54,c1,c2,c3,c4)
end

function alg_cache(alg::SSPRK54,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{false}})
  β10 = 0.391752226571890
  α20 = 0.444370493651235
  α21 = 0.555629506348765
  β21 = 0.368410593050371
  α30 = 0.620101851488403
  α32 = 0.379898148511597
  β32 = 0.251891774271694
  α40 = 0.178079954393132
  α43 = 0.821920045606868
  β43 = 0.544974750228521
  α52 = 0.517231671970585
  α53 = 0.096059710526147
  β53 = 0.063692468666290
  α54 = 0.386708617503269
  β54 = 0.226007483236906
  c1 = 0.39175222700391998
  c2 = 0.58607968896779994
  c3 = 0.47454236302687003
  c4 = 0.93501063100923998
  SSPRK54ConstantCache(β10,α20,α21,β21,α30,α32,β32,α40,α43,β43,α52,α53,β53,α54,β54,c1,c2,c3,c4)
end


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
