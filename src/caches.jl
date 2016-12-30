immutable ODEEmptyCache <: DECache end

alg_cache{F}(alg::OrdinaryDiffEqAlgorithm,prob,callback::F) = ODEEmptyCache()

type ExplicitRKCache{uType,rateType,uEltypeNoUnits,ksEltype} <: DECache
  u::uType
  tmp::uType
  utilde::rateType
  uEEst::rateType
  atmp::uEltypeNoUnits
  uprev::uType
  kprev::ksEltype
  utmp::uType
  kk::Vector{ksEltype}
end

function alg_cache{uType<:AbstractArray}(alg::ExplicitRK,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev)
  kk = Vector{typeof(rate_prototype)}(0)
  for i = 1:tableau.stages
    push!(kk,similar(rate_prototype))
  end
  utilde = similar(rate_prototype)
  tmp = similar(u)
  atmp = similar(u,uEltypeNoUnits)
  utmp = zeros(u)
  uEEst = similar(rate_prototype)
  ExplicitRKCache(u,tmp,utilde,uEEst,atmp,uprev,kprev,utmp,kk)
end

alg_cache{uType}(alg::ExplicitRK,u::uType,rate_prototype,uEltypeNoUnits,tableau,uprev,kprev) = ODEEmptyCache()
