mutable struct DEOptions{absType,relType,QT,tType,F2,F3,F4,F5,F6,tstopsType,discType,ECType,SType,MI,tcache,savecache,disccache}
  maxiters::MI
  timeseries_steps::Int
  save_everystep::Bool
  adaptive::Bool
  abstol::absType
  reltol::relType
  gamma::QT
  qmax::QT
  qmin::QT
  qsteady_max::QT
  qsteady_min::QT
  failfactor::QT
  dtmax::tType
  dtmin::tType
  internalnorm::F2
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
  progress_message::F5
  timeseries_errors::Bool
  dense_errors::Bool
  beta1::QT
  beta2::QT
  qoldinit::QT
  dense::Bool
  save_start::Bool
  save_end::Bool
  callback::F3
  isoutofdomain::F4
  unstable_check::F6
  verbose::Bool
  calck::Bool
  force_dtmin::Bool
  advance_to_tstop::Bool
  stop_at_next_tstop::Bool
end

"""
    ODEIntegrator
Fundamental `struct` allowing interactively stepping through the numerical solving of a differential equation.
The full documentation is hosted here: [http://docs.juliadiffeq.org/latest/basics/integrator.html](http://docs.juliadiffeq.org/latest/basics/integrator.html). This docstring
describes basic functionality only!

Initialize using `integrator = init(prob::ODEProblem, alg; kwargs...)`. The keyword args which are accepted are the same
[Common Solver Options](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html#Common-Solver-Options-1)
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
mutable struct ODEIntegrator{algType<:OrdinaryDiffEqAlgorithm,uType,tType,pType,eigenType,QT,tdirType,ksEltype,SolType,F,ProgressType,CacheType,O,FSALType} <: AbstractODEIntegrator
  sol::SolType
  u::uType
  k::ksEltype
  t::tType
  dt::tType
  f::F
  p::pType
  uprev::uType
  uprev2::uType
  tprev::tType
  alg::algType
  dtcache::tType
  dtchangeable::Bool
  dtpropose::tType
  tdir::tdirType
  eigen_est::eigenType
  EEst::QT
  qold::QT
  q11::QT
  erracc::QT
  dtacc::tType
  success_iter::Int
  iter::Int
  saveiter::Int
  saveiter_dense::Int
  prog::ProgressType
  cache::CacheType
  kshortsize::Int
  force_stepfail::Bool
  last_stepfail::Bool
  just_hit_tstop::Bool
  event_last_time::Bool
  accept_step::Bool
  isout::Bool
  reeval_fsal::Bool
  u_modified::Bool
  opts::O
  fsalfirst::FSALType
  fsallast::FSALType

  function ODEIntegrator{algType,uType,tType,pType,eigenType,tTypeNoUnits,tdirType,ksEltype,SolType,
                F,ProgressType,CacheType,O,FSALType}(
                sol,u,k,t,dt,f,p,uprev,uprev2,tprev,
      alg,dtcache,dtchangeable,dtpropose,tdir,
      eigen_est,EEst,qold,q11,erracc,dtacc,success_iter,
      iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,force_stepfail,last_stepfail,just_hit_tstop,
      event_last_time,accept_step,isout,reeval_fsal,u_modified,opts) where {algType,uType,tType,pType,eigenType,tTypeNoUnits,tdirType,ksEltype,SolType,
                                     F,ProgressType,CacheType,O,FSALType}

      new{algType,uType,tType,pType,eigenType,tTypeNoUnits,tdirType,ksEltype,SolType,
                  F,ProgressType,CacheType,O,FSALType}(
                  sol,u,k,t,dt,f,p,uprev,uprev2,tprev,
      alg,dtcache,dtchangeable,dtpropose,tdir,
      eigen_est,EEst,qold,q11,erracc,dtacc,success_iter,
      iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,force_stepfail,last_stepfail,just_hit_tstop,
      event_last_time,accept_step,isout,reeval_fsal,u_modified,opts) # Leave off fsalfirst and last
  end
end

# When this is changed, DelayDiffEq.jl must be changed as well!
