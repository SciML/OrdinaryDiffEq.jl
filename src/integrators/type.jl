type DEOptions{uEltype,uEltypeNoUnits,tTypeNoUnits,tType,F2,F3,F4,F5,tstopsType,ECType}
  maxiters::Int
  timeseries_steps::Int
  save_timeseries::Bool
  adaptive::Bool
  abstol::uEltype
  reltol::uEltypeNoUnits
  gamma::tTypeNoUnits
  qmax::tTypeNoUnits
  qmin::tTypeNoUnits
  dtmax::tType
  dtmin::tType
  internalnorm::F2
  tstops::tstopsType
  saveat::tstopsType
  userdata::ECType
  progress::Bool
  progress_steps::Int
  progress_name::String
  progress_message::F5
  timeseries_errors::Bool
  dense_errors::Bool
  beta1::tTypeNoUnits
  beta2::tTypeNoUnits
  qoldinit::tTypeNoUnits
  dense::Bool
  callback::F3
  isoutofdomain::F4
  calck::Bool
  advance_to_tstop::Bool
  stop_at_next_tstop::Bool
end

type ODEIntegrator{algType<:OrdinaryDiffEqAlgorithm,uType<:Union{AbstractArray,Number},tType,tTypeNoUnits,tdirType,ksEltype,SolType,rateType,F,ProgressType,CacheType,O} <: AbstractODEIntegrator
  sol::SolType
  u::uType
  k::ksEltype
  t::tType
  dt::tType
  f::F
  uprev::uType
  tprev::tType
  alg::algType
  rate_prototype::rateType
  notsaveat_idxs::Vector{Int}
  calcprevs::Bool
  dtcache::tType
  dtchangeable::Bool

  dtpropose::tType
  dt_mod::tTypeNoUnits
  tdir::tdirType
  EEst::tTypeNoUnits
  qold::tTypeNoUnits
  q11::tTypeNoUnits
  iter::Int
  saveiter::Int
  saveiter_dense::Int
  prog::ProgressType
  cache::CacheType
  kshortsize::Int
  just_hit_tstop::Bool
  accept_step::Bool
  reeval_fsal::Bool
  u_modified::Bool
  opts::O
  fsalfirst::rateType
  fsallast::rateType

  ODEIntegrator(sol,u,k,t,dt,f,uprev,tprev,
      alg,rate_prototype,notsaveat_idxs,calcprevs,dtcache,dtchangable,dtpropose,dt_mod,tdir,
      EEst,qold,q11,iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,just_hit_tstop,accept_step,reeval_fsal,u_modified,opts) = new(
      sol,u,k,t,dt,f,uprev,tprev,
      alg,rate_prototype,notsaveat_idxs,calcprevs,dtcache,dtchangable,dtpropose,dt_mod,tdir,
      EEst,qold,q11,iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,just_hit_tstop,accept_step,reeval_fsal,u_modified,opts) # Leave off fsalfirst and last
end
