@cache struct SSPRK22Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType # superfluous, only needed for callbacks...
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

struct SSPRK22ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK22,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  SSPRK22Cache(u,uprev,k,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK22,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = SSPRK22ConstantCache()

@cache struct SSPRK33Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType # superfluous, only needed for callbacks...
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

struct SSPRK33ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK33,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  SSPRK33Cache(u,uprev,k,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK33,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = SSPRK33ConstantCache()

@cache struct SSPRK53Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK53ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  α30::T
  α32::T
  α40::T
  α43::T
  α52::T
  α54::T
  β10::T
  β21::T
  β32::T
  β43::T
  β54::T
  c1::T2
  c2::T2
  c3::T2
  c4::T2

  function SSPRK53ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    α30 = T(0.355909775063327)
    α32 = T(0.644090224936674)
    α40 = T(0.367933791638137)
    α43 = T(0.632066208361863)
    α52 = T(0.237593836598569)
    α54 = T(0.762406163401431)
    β10 = T(0.377268915331368)
    β21 = T(0.377268915331368)
    β32 = T(0.242995220537396)
    β43 = T(0.238458932846290)
    β54 = T(0.287632146308408)
    c1 = T2(0.377268915331368)
    c2 = T2(0.754537830662736)
    c3 = T2(0.728985661612188)
    c4 = T2(0.699226135931670)

    new{T,T2}(α30, α32, α40, α43, α52, α54, β10, β21, β32, β43, β54, c1, c2, c3, c4)
  end
end

function alg_cache(alg::SSPRK53,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = SSPRK53ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK53Cache(u,uprev,k,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK53,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  SSPRK53ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct SSPRK63Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  u₂::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK63ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  α40::T
  α41::T
  α43::T
  α62::T
  α65::T
  β10::T
  β21::T
  β32::T
  β43::T
  β54::T
  β65::T
  c1::T2
  c2::T2
  c3::T2
  c4::T2
  c5::T2

  function SSPRK63ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    α40 = T(0.476769811285196)
    α41 = T(0.098511733286064)
    α43 = T(0.424718455428740)
    α62 = T(0.155221702560091)
    α65 = T(0.844778297439909)
    β10 = T(0.284220721334261)
    β21 = T(0.284220721334261)
    β32 = T(0.284220721334261)
    β43 = T(0.120713785765930)
    β54 = T(0.284220721334261)
    β65 = T(0.240103497065900)
    c1 = T2(0.284220721334261)
    c2 = T2(0.568441442668522)
    c3 = T2(0.852662164002783)
    c4 = T2(0.510854218958172)
    c5 = T2(0.795074940292433)

    new{T,T2}(α40, α41, α43, α62, α65, β10, β21, β32, β43, β54, β65, c1, c2, c3, c4, c5)
  end
end

function alg_cache(alg::SSPRK63,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  u₂ = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = SSPRK63ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK63Cache(u,uprev,k,tmp,u₂,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK63,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  SSPRK63ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct SSPRK73Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  u₁::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK73ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  α40::T
  α43::T
  α50::T
  α51::T
  α54::T
  α73::T
  α76::T
  β10::T
  β21::T
  β32::T
  β43::T
  β54::T
  β65::T
  β76::T
  c1::T2
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2

  function SSPRK73ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    α40 = T(0.184962588071072)
    α43 = T(0.815037411928928)
    α50 = T(0.180718656570380)
    α51 = T(0.314831034403793)
    α54 = T(0.504450309025826)
    α73 = T(0.120199000000000)
    α76 = T(0.879801000000000)
    β10 = T(0.233213863663009)
    β21 = T(0.233213863663009)
    β32 = T(0.233213863663009)
    β43 = T(0.190078023865845)
    β54 = T(0.117644805593912)
    β65 = T(0.233213863663009)
    β76 = T(0.205181790464579)
    c1 = T2(0.233213863663009)
    c2 = T2(0.466427727326018)
    c3 = T2(0.699641590989027)
    c4 = T2(0.760312095463379)
    c5 = T2(0.574607439040817)
    c6 = T2(0.807821302703826)

    new{T,T2}(α40, α43, α50, α51, α54, α73, α76, β10, β21, β32, β43, β54, β65, β76, c1, c2, c3, c4, c5, c6)
  end
end

function alg_cache(alg::SSPRK73,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  u₁ = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = SSPRK73ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK73Cache(u,uprev,k,tmp,u₁,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK73,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  SSPRK73ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct SSPRK83Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType
  u₂::uType
  u₃::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK83ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  α50::T
  α51::T
  α54::T
  α61::T
  α65::T
  α72::T
  α73::T
  α76::T
  β10::T
  β21::T
  β32::T
  β43::T
  β54::T
  β65::T
  β76::T
  β87::T
  c1::T2
  c2::T2
  c3::T2
  c4::T2
  c5::T2
  c6::T2
  c7::T2

  function SSPRK83ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    α50 = T(0.421366967085359)
    α51 = T(0.005949401107575)
    α54 = T(0.572683631807067)
    α61 = T(0.004254010666365)
    α65 = T(0.995745989333635)
    α72 = T(0.104380143093325)
    α73 = T(0.243265240906726)
    α76 = T(0.652354615999950)
    β10 = T(0.195804015330143)
    β21 = T(0.195804015330143)
    β32 = T(0.195804015330143)
    β43 = T(0.195804015330143)
    β54 = T(0.112133754621673)
    β65 = T(0.194971062960412)
    β76 = T(0.127733653231944)
    β87 = T(0.195804015330143)
    c1 = T2(0.195804015330143)
    c2 = T2(0.391608030660286)
    c3 = T2(0.587412045990429)
    c4 = T2(0.783216061320572)
    c5 = T2(0.561833689734037)
    c6 = T2(0.755247658555329)
    c7 = T2(0.804195984669857)

    new{T,T2}(α50, α51, α54, α61, α65, α72, α73, α76, β10, β21, β32, β43, β54, β65, β76, β87, c1, c2, c3, c4, c5, c6, c7)
  end
end

function alg_cache(alg::SSPRK83,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  u₂ = similar(u)
  u₃ = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = SSPRK83ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK83Cache(u,uprev,k,tmp,u₂,u₃,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK83,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  SSPRK83ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct SSPRK432Cache{uType,rateType,uNoUnitsType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType # superfluous, only needed for callbacks
  fsalfirst::rateType
  utilde::uType
  atmp::uNoUnitsType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

struct SSPRK432ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK432,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  SSPRK432Cache(u,uprev,k,tmp,fsalfirst,utilde,atmp,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK432,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = SSPRK432ConstantCache()

@cache struct SSPRK932Cache{uType,rateType,uNoUnitsType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  tmp::uType # superfluous, only needed for callbacks
  fsalfirst::rateType
  utilde::uType
  atmp::uNoUnitsType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

struct SSPRK932ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK932,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  utilde = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  SSPRK932Cache(u,uprev,k,tmp,fsalfirst,utilde,atmp,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK932,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = SSPRK932ConstantCache()

@cache struct SSPRK54Cache{uType,rateType,StageLimiter,StepLimiter,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₃::rateType
  u₂::uType
  u₃::uType
  tmp::uType # should be u₄, but tmp is needed for callbacks
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  tab::TabType
end

struct SSPRK54ConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  β10::T
  α20::T
  α21::T
  β21::T
  α30::T
  α32::T
  β32::T
  α40::T
  α43::T
  β43::T
  α52::T
  α53::T
  β53::T
  α54::T
  β54::T
  c1::T2
  c2::T2
  c3::T2
  c4::T2

  function SSPRK54ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}
    β10 = T(0.391752226571890)
    α20 = T(0.444370493651235)
    α21 = T(0.555629506348765)
    β21 = T(0.368410593050371)
    α30 = T(0.620101851488403)
    α32 = T(0.379898148511597)
    β32 = T(0.251891774271694)
    α40 = T(0.178079954393132)
    α43 = T(0.821920045606868)
    β43 = T(0.544974750228521)
    α52 = T(0.517231671970585)
    α53 = T(0.096059710526147)
    β53 = T(0.063692468666290)
    α54 = T(0.386708617503269)
    β54 = T(0.226007483236906)
    c1 = T2(0.391752226571890)
    c2 = T2(0.586079689311540)
    c3 = T2(0.474542363121400)
    c4 = T2(0.935010630967653)

    new{T,T2}(β10, α20, α21, β21, α30, α32, β32, α40, α43, β43, α52, α53, β53, α54, β54, c1, c2, c3, c4)
  end
end

function alg_cache(alg::SSPRK54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  u₂ = similar(u)
  u₃ = similar(u)
  tmp = similar(u)
  k = zero(rate_prototype)
  k₃ = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  tab = SSPRK54ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  SSPRK54Cache(u,uprev,k,k₃,u₂,u₃,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!,tab)
end

function alg_cache(alg::SSPRK54,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}})
  SSPRK54ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end

@cache struct SSPRK104Cache{uType,rateType,StageLimiter,StepLimiter} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₄::rateType
  tmp::uType
  fsalfirst::rateType
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end

struct SSPRK104ConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::SSPRK104,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{true}})
  tmp = similar(u)
  k = zero(rate_prototype)
  k₄ = zero(rate_prototype)
  fsalfirst = zero(rate_prototype)
  SSPRK104Cache(u,uprev,k,k₄,tmp,fsalfirst,alg.stage_limiter!,alg.step_limiter!)
end

alg_cache(alg::SSPRK104,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Type{Val{false}}) = SSPRK104ConstantCache()
