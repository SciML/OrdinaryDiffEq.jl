mutable struct DEOptions{absType,relType,QT,tType,Controller,F1,F2,F3,F4,F5,F6,F7,tstopsType,discType,ECType,SType,MI,tcache,savecache,disccache}
  maxiters::MI
  save_everystep::Bool
  adaptive::Bool
  abstol::absType
  reltol::relType
  gamma::QT
  qmax::QT
  qmin::QT
  qsteady_max::QT
  qsteady_min::QT
  qoldinit::QT
  failfactor::QT
  dtmax::tType
  dtmin::tType
  controller::Controller
  internalnorm::F1
  internalopnorm::F2
  save_idxs::SType
  tstops::tstopsType
  saveat::tstopsType
  d_discontinuities::discType
  tstops_cache::tcache
  saveat_cache::savecache
  d_discontinuities_cache::disccache
  userdata::ECType
  progress::Bool
  progress_steps::Int
  progress_name::String
  progress_message::F6
  timeseries_errors::Bool
  dense_errors::Bool
  dense::Bool
  save_on::Bool
  save_start::Bool
  save_end::Bool
  save_end_user::F3
  callback::F4
  isoutofdomain::F5
  unstable_check::F7
  verbose::Bool
  calck::Bool
  force_dtmin::Bool
  advance_to_tstop::Bool
  stop_at_next_tstop::Bool
end

"""
    ODEIntegrator
Fundamental `struct` allowing interactively stepping through the numerical solving of a differential equation.
The full documentation is hosted here:
[https://diffeq.sciml.ai/latest/basics/integrator/](https://diffeq.sciml.ai/latest/basics/integrator/).
This docstring describes basic functionality only!

Initialize using `integrator = init(prob::ODEProblem, alg; kwargs...)`. The keyword args which are accepted are the same
[common solver options](https://diffeq.sciml.ai/latest/basics/common_solver_opts/)
used by `solve`.


For reference, relevant fields of the `ODEIntegrator` are:

* `t` - time of the proposed step
* `u` - value at the proposed step
* `opts` - common solver options
* `alg` - the algorithm associated with the solution
* `f` - the function being solved
* `sol` - the current state of the solution
* `tprev` - the last timepoint
* `uprev` - the value at the last timepoint

`opts` holds all of the common solver options, and can be mutated to change the solver characteristics.
For example, to modify the absolute tolerance for the future timesteps, one can do:
```julia
integrator.opts.abstol = 1e-9
```
For more info see the linked documentation page.
"""
mutable struct ODEIntegrator{algType<:Union{OrdinaryDiffEqAlgorithm,DAEAlgorithm},IIP,uType,duType,tType,pType,eigenType,EEstT,QT,tdirType,ksEltype,SolType,F,CacheType,O,FSALType,EventErrorType,CallbackCacheType,IA} <: DiffEqBase.AbstractODEIntegrator{algType,IIP,uType,tType}
  sol::SolType
  u::uType
  du::duType
  k::ksEltype
  t::tType
  dt::tType
  f::F
  p::pType
  uprev::uType
  uprev2::uType
  duprev::duType
  tprev::tType
  alg::algType
  dtcache::tType
  dtchangeable::Bool
  dtpropose::tType
  tdir::tdirType
  eigen_est::eigenType
  EEst::EEstT
  qold::QT
  q11::QT
  erracc::QT
  dtacc::tType
  success_iter::Int
  iter::Int
  saveiter::Int
  saveiter_dense::Int
  cache::CacheType
  callback_cache::CallbackCacheType
  kshortsize::Int
  force_stepfail::Bool
  last_stepfail::Bool
  just_hit_tstop::Bool
  do_error_check::Bool
  event_last_time::Int
  vector_event_last_time::Int
  last_event_error::EventErrorType
  accept_step::Bool
  isout::Bool
  reeval_fsal::Bool
  u_modified::Bool
  reinitialize::Bool
  isdae::Bool
  opts::O
  destats::DiffEqBase.DEStats
  initializealg::IA
  fsalfirst::FSALType
  fsallast::FSALType

  function ODEIntegrator{algType,IIP,uType,duType,tType,pType,eigenType,EEstT,tTypeNoUnits,tdirType,ksEltype,SolType,
                F,CacheType,O,FSALType,EventErrorType,CallbackCacheType,InitializeAlgType}(
                sol,u,du,k,t,dt,f,p,uprev,uprev2,duprev,tprev,
      alg,dtcache,dtchangeable,dtpropose,tdir,
      eigen_est,EEst,qold,q11,erracc,dtacc,success_iter,
      iter,saveiter,saveiter_dense,cache,callback_cache,
      kshortsize,force_stepfail,last_stepfail,just_hit_tstop,
      do_error_check,
      event_last_time,vector_event_last_time,last_event_error,
      accept_step,isout,reeval_fsal,u_modified,reinitialize,isdae,
      opts,destats,initializealg) where {algType,IIP,uType,duType,tType,pType,
                                         eigenType,EEstT,tTypeNoUnits,tdirType,
                                         ksEltype,SolType,F,CacheType,O,
                                         FSALType,EventErrorType,
                                         CallbackCacheType,InitializeAlgType}

      new{algType,IIP,uType,duType,tType,pType,eigenType,EEstT,tTypeNoUnits,tdirType,ksEltype,SolType,
                  F,CacheType,O,FSALType,EventErrorType,CallbackCacheType,InitializeAlgType}(
                  sol,u,du,k,t,dt,f,p,uprev,uprev2,duprev,tprev,
      alg,dtcache,dtchangeable,dtpropose,tdir,
      eigen_est,EEst,qold,q11,erracc,dtacc,success_iter,
      iter,saveiter,saveiter_dense,cache,callback_cache,
      kshortsize,force_stepfail,last_stepfail,just_hit_tstop,
      do_error_check,
      event_last_time,vector_event_last_time,last_event_error,
      accept_step,isout,reeval_fsal,u_modified,reinitialize,isdae,
      opts,destats,initializealg) # Leave off fsalfirst and last
  end
end

# When this is changed, DelayDiffEq.jl must be changed as well!
