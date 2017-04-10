type DEOptions{uEltype,uEltypeNoUnits,tTypeNoUnits,tType,F2,F3,F4,F5,F6,tstopsType,ECType,SType}
  maxiters::Int
  timeseries_steps::Int
  save_everystep::Bool
  adaptive::Bool
  abstol::uEltype
  reltol::uEltypeNoUnits
  gamma::tTypeNoUnits
  qmax::tTypeNoUnits
  qmin::tTypeNoUnits
  dtmax::tType
  dtmin::tType
  internalnorm::F2
  save_idxs::SType
  tstops::tstopsType
  saveat::tstopsType
  d_discontinuities::tstopsType
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
  save_start::Bool
  callback::F3
  isoutofdomain::F4
  unstable_check::F6
  verbose::Bool
  calck::Bool
  force_dtmin::Bool
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
  uprev2::uType
  tprev::tType
  alg::algType
  rate_prototype::rateType
  notsaveat_idxs::Vector{Int}
  dtcache::tType
  dtchangeable::Bool
  dtpropose::tType
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
  isout::Bool
  reeval_fsal::Bool
  u_modified::Bool
  opts::O
  fsalfirst::rateType
  fsallast::rateType

  ODEIntegrator(sol,u,k,t,dt,f,uprev,uprev2,tprev,
      alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,dtpropose,tdir,
      EEst,qold,q11,
      iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,just_hit_tstop,accept_step,isout,reeval_fsal,u_modified,opts) = new(
      sol,u,k,t,dt,f,uprev,uprev2,tprev,
      alg,rate_prototype,notsaveat_idxs,dtcache,dtchangeable,dtpropose,tdir,
      EEst,qold,q11,
      iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,just_hit_tstop,accept_step,isout,reeval_fsal,u_modified,opts) # Leave off fsalfirst and last
end

# When this is changed, DelayDiffEq.jl must be changed as well!
