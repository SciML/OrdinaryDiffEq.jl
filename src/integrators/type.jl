mutable struct DEOptions{absType,relType,tTypeNoUnits,tType,F2,F3,F4,F5,F6,tstopsType,discType,ECType,SType,MI,tcache,savecache,disccache}
  maxiters::MI
  timeseries_steps::Int
  save_everystep::Bool
  adaptive::Bool
  abstol::absType
  reltol::relType
  gamma::tTypeNoUnits
  qmax::tTypeNoUnits
  qmin::tTypeNoUnits
  qsteady_max::tTypeNoUnits
  qsteady_min::tTypeNoUnits
  failfactor::tTypeNoUnits
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
  beta1::tTypeNoUnits
  beta2::tTypeNoUnits
  qoldinit::tTypeNoUnits
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

## Basic Usage
To initialize an integrator, use the syntax:
```julia
integrator = init(prob::ODEProblem, alg; kwargs...)
```
The keyword args which are accepted are the same 
[Common Solver Options](http://docs.juliadiffeq.org/latest/basics/common_solver_opts.html#Common-Solver-Options-1)
used by `solve`.
 
You can then manually step the integrator using `step!(integrator)`. Other stepping options include:
```julia
for i in take(integrator,n) end
for i in integrator end
for (u,t) in tuples(integrator) end
for (tprev,uprev,u,t) in intervals(integrator) end
```

## Relevant Fields
The `ODEIntegrator` type holds all of the information for the intermediate solution 
of the differential equation. Useful fields are:


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

The `sol` field holds the current solution.
This current solution includes the interpolation function if available, 
and thus `integrator.sol(t)` lets one interpolate efficiently 
over the whole current solution. Additionally, a 
"current interval interpolation function" is provided 
on the integrator type via `integrator(t)`. This uses only the solver 
information from the interval `[tprev,t]` to compute the interpolation, 
and is allowed to extrapolate beyond that interval.

## Basic Stepping Control

* `u_modified!(integrator,bool)`: Bool which states whether a change to u occurred, allowing the solver to handle the discontinuity. By default, this is assumed to be `true` if a callback is used. This will result in the re-calculation of the derivative at `t+dt`, which is not necessary if the algorithm is FSAL and u does not experience a discontinuous change at the end of the interval. Thus if u is unmodified in a callback, a single call to the derivative calculation can be eliminated by `u_modified!(integrator,false)`.
* `get_proposed_dt(integrator)`: Gets the proposed dt for the next timestep.
* `set_proposed_dt!(integrator,dt)`: Sets the proposed dt for the next timestep.
* `proposed_dt(integrator)`: Returns the dt of the proposed step.
* `terminate!(integrator)`: Terminates the integrator by emptying tstops. This can be used in events and callbacks to immediately end the solution process.
* `add_tstop!(integrator,t)`: Adds a tstop at time t.
* `add_saveat!(integrator,t)`: Adds a saveat time point at t.
"""
mutable struct ODEIntegrator{algType<:OrdinaryDiffEqAlgorithm,uType,tType,pType,tTypeNoUnits,tdirType,ksEltype,SolType,F,ProgressType,CacheType,O,FSALType} <: AbstractODEIntegrator
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
  EEst::tTypeNoUnits
  qold::tTypeNoUnits
  q11::tTypeNoUnits
  erracc::tTypeNoUnits
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
  accept_step::Bool
  isout::Bool
  reeval_fsal::Bool
  u_modified::Bool
  opts::O
  fsalfirst::FSALType
  fsallast::FSALType

  function (::Type{ODEIntegrator{algType,uType,tType,pType,tTypeNoUnits,tdirType,ksEltype,SolType,
                F,ProgressType,CacheType,O,FSALType}}){algType,uType,tType,pType,tTypeNoUnits,tdirType,ksEltype,SolType,
                F,ProgressType,CacheType,O,FSALType}(
                sol,u,k,t,dt,f,p,uprev,uprev2,tprev,
      alg,dtcache,dtchangeable,dtpropose,tdir,
      EEst,qold,q11,erracc,dtacc,success_iter,
      iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,force_stepfail,last_stepfail,just_hit_tstop,
      accept_step,isout,reeval_fsal,u_modified,opts)

      new{algType,uType,tType,pType,tTypeNoUnits,tdirType,ksEltype,SolType,
                  F,ProgressType,CacheType,O,FSALType}(
                  sol,u,k,t,dt,f,p,uprev,uprev2,tprev,
      alg,dtcache,dtchangeable,dtpropose,tdir,
      EEst,qold,q11,erracc,dtacc,success_iter,
      iter,saveiter,saveiter_dense,prog,cache,
      kshortsize,force_stepfail,last_stepfail,just_hit_tstop,
      accept_step,isout,reeval_fsal,u_modified,opts) # Leave off fsalfirst and last
  end
end

# When this is changed, DelayDiffEq.jl must be changed as well!
