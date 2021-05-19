save_idxsinitialize(integrator,cache::OrdinaryDiffEqCache,::Type{uType}) where {uType} =
               error("This algorithm does not have an initialization function")

function loopheader!(integrator)
  # Apply right after iterators / callbacks

  # Accept or reject the step
  if integrator.iter > 0
    if ((integrator.opts.adaptive && integrator.accept_step) || !integrator.opts.adaptive) && !integrator.force_stepfail
      integrator.success_iter += 1
      apply_step!(integrator)
    elseif integrator.opts.adaptive && !integrator.accept_step
      if integrator.isout
        integrator.dt = integrator.dt*integrator.opts.qmin
      elseif !integrator.force_stepfail
        step_reject_controller!(integrator,integrator.alg)
      end
    end
  end

  integrator.iter += 1
  choose_algorithm!(integrator,integrator.cache)
  fix_dt_at_bounds!(integrator)
  modify_dt_for_tstops!(integrator)
  integrator.force_stepfail = false
end

last_step_failed(integrator::ODEIntegrator) =
  integrator.last_stepfail && !integrator.opts.adaptive

function modify_dt_for_tstops!(integrator)
  if has_tstop(integrator)
    tdir_t = integrator.tdir * integrator.t
    tdir_tstop = first_tstop(integrator)
    if integrator.opts.adaptive
      integrator.dt = integrator.tdir * min(abs(integrator.dt), abs(tdir_tstop - tdir_t)) # step! to the end
    elseif iszero(integrator.dtcache) && integrator.dtchangeable
      integrator.dt = integrator.tdir * abs(tdir_tstop - tdir_t)
    elseif integrator.dtchangeable && !integrator.force_stepfail
      # always try to step! with dtcache, but lower if a tstop
      # however, if force_stepfail then don't set to dtcache, and no tstop worry
      integrator.dt = integrator.tdir * min(abs(integrator.dtcache), abs(tdir_tstop - tdir_t)) # step! to the end
    end
  end
end

# Want to extend savevalues! for DDEIntegrator
savevalues!(integrator::ODEIntegrator, force_save = false, reduce_size = true) =
  _savevalues!(integrator, force_save, reduce_size)

function _savevalues!(integrator, force_save, reduce_size)::Tuple{Bool,Bool}
  saved, savedexactly = false, false
  !integrator.opts.save_on && return saved, savedexactly
  tdir_t = integrator.tdir * integrator.t
  while !isempty(integrator.opts.saveat) && first(integrator.opts.saveat) <= tdir_t # Perform saveat
    integrator.saveiter += 1; saved = true
    curt = integrator.tdir * pop!(integrator.opts.saveat)
    if curt!=integrator.t # If <t, interpolate
      DiffEqBase.addsteps!(integrator)
      Θ = (curt - integrator.tprev)/integrator.dt
      val = ode_interpolant(Θ,integrator,integrator.opts.save_idxs,Val{0}) # out of place, but no force copy later
      copyat_or_push!(integrator.sol.t,integrator.saveiter,curt)
      save_val = val
      copyat_or_push!(integrator.sol.u,integrator.saveiter,save_val,Val{false})
      if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    else # ==t, just save
      savedexactly = true
      copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
      if integrator.opts.save_idxs === nothing
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
      else
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
      end
      if typeof(integrator.alg) <: FunctionMap || integrator.opts.dense
        integrator.saveiter_dense +=1
        if integrator.opts.dense
          if integrator.opts.save_idxs === nothing
            copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
          else
            copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,[k[integrator.opts.save_idxs] for k in integrator.k],Val{false})
          end
        end
      end
      if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    end
  end
  if force_save || (integrator.opts.save_everystep && (isempty(integrator.sol.t) || (integrator.t !== integrator.sol.t[end])))
    integrator.saveiter += 1; saved, savedexactly = true, true
    if integrator.opts.save_idxs === nothing
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    else
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
    end
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if typeof(integrator.alg) <: FunctionMap || integrator.opts.dense
      integrator.saveiter_dense +=1
      if integrator.opts.dense
        if integrator.opts.save_idxs === nothing
          copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
        else
          copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,[k[integrator.opts.save_idxs] for k in integrator.k],Val{false})
        end
      end
    end
    if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
      copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
    end
  end
  reduce_size && resize!(integrator.k,integrator.kshortsize)
  return saved, savedexactly
end

# Want to extend postamble! for DDEIntegrator
postamble!(integrator::ODEIntegrator) = _postamble!(integrator)

function _postamble!(integrator)
  solution_endpoint_match_cur_integrator!(integrator)
  resize!(integrator.sol.t,integrator.saveiter)
  resize!(integrator.sol.u,integrator.saveiter)
  if !(integrator.sol isa DAESolution)
    resize!(integrator.sol.k,integrator.saveiter_dense)
  end
  if integrator.opts.progress
    @logmsg(LogLevel(-1),
    integrator.opts.progress_name,
    _id = :OrdinaryDiffEq,
    message=integrator.opts.progress_message(integrator.dt,integrator.u,integrator.p,integrator.t),
    progress="done")
  end
end

function solution_endpoint_match_cur_integrator!(integrator)
  if integrator.opts.save_end &&
     (integrator.saveiter == 0 || integrator.sol.t[integrator.saveiter] !=  integrator.t &&
     ((integrator.opts.save_end_user isa Bool && integrator.opts.save_end_user) ||
       integrator.t ∈ integrator.opts.saveat_cache ||
       integrator.t == integrator.sol.prob.tspan[2] ||
       isempty(integrator.opts.saveat_cache)))

    integrator.saveiter += 1
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if integrator.opts.save_idxs === nothing
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    else
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
    end
    if typeof(integrator.alg) <: FunctionMap || integrator.opts.dense
      integrator.saveiter_dense +=1
      if integrator.opts.dense
        if integrator.opts.save_idxs === nothing
          copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
        else
          copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,[k[integrator.opts.save_idxs] for k in integrator.k],Val{false})
        end
      end
    end
    if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
      copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
    end
  end
end

# Want to extend loopfooter! for DDEIntegrator
loopfooter!(integrator::ODEIntegrator) = _loopfooter!(integrator)

function _loopfooter!(integrator)

  # Carry-over from callback
  # This is set to true if u_modified requires callback FSAL reset
  # But not set to false when reset so algorithms can check if reset occurred
  integrator.reeval_fsal = false
  integrator.u_modified = false
  integrator.do_error_check = true
  ttmp = integrator.t + integrator.dt
  if integrator.force_stepfail
      if integrator.opts.adaptive
        integrator.dt = integrator.dt/integrator.opts.failfactor
      elseif integrator.last_stepfail
        return
      end
      integrator.last_stepfail = true
      integrator.accept_step = false
  elseif integrator.opts.adaptive
    q = stepsize_controller!(integrator,integrator.alg)
    integrator.isout = integrator.opts.isoutofdomain(integrator.u,integrator.p,ttmp)
    integrator.accept_step = (!integrator.isout && accept_step_controller(integrator, integrator.opts.controller)) || (integrator.opts.force_dtmin && abs(integrator.dt) <= timedepentdtmin(integrator))
    if integrator.accept_step # Accept
      integrator.destats.naccept += 1
      integrator.last_stepfail = false
      dtnew = DiffEqBase.value(step_accept_controller!(integrator,integrator.alg,q)) * oneunit(integrator.dt)
      integrator.tprev = integrator.t
      if integrator.t isa AbstractFloat && has_tstop(integrator)
        tstop = integrator.tdir * first_tstop(integrator)
        abs(ttmp - tstop) < 100eps(max(integrator.t,tstop)/oneunit(integrator.t))*oneunit(integrator.t) ?
                                  (integrator.t = tstop) : (integrator.t = ttmp)
      else
        integrator.t = ttmp
      end
      calc_dt_propose!(integrator,dtnew)
      handle_callbacks!(integrator)
    else # Reject
      integrator.destats.nreject += 1
    end
  elseif !integrator.opts.adaptive #Not adaptive
    integrator.destats.naccept += 1
    integrator.tprev = integrator.t
    if integrator.t isa AbstractFloat && has_tstop(integrator)
      tstop = integrator.tdir * first_tstop(integrator)
      abs(ttmp - tstop) < 100eps(integrator.t/oneunit(integrator.t))*oneunit(integrator.t) ?
                                  (integrator.t = tstop) : (integrator.t = ttmp)
    else
      integrator.t = ttmp
    end
    integrator.last_stepfail = false
    integrator.accept_step = true
    integrator.dtpropose = integrator.dt
    handle_callbacks!(integrator)
  end
  if integrator.opts.progress && integrator.iter%integrator.opts.progress_steps==0
    t1,t2 = integrator.sol.prob.tspan
    @logmsg(LogLevel(-1),
    integrator.opts.progress_name,
    _id = :OrdinaryDiffEq,
    message=integrator.opts.progress_message(integrator.dt,integrator.u,integrator.p,integrator.t),
    progress=(integrator.t-t1)/(t2-t1))
  end

  # Take value because if t is dual then maxeig can be dual
  if integrator.cache isa CompositeCache
    cur_eigen_est = integrator.opts.internalnorm(DiffEqBase.value(integrator.eigen_est),integrator.t)
    cur_eigen_est > integrator.destats.maxeig && (integrator.destats.maxeig = cur_eigen_est)
  end
  nothing
end

function handle_callbacks!(integrator)
  discrete_callbacks = integrator.opts.callback.discrete_callbacks
  continuous_callbacks = integrator.opts.callback.continuous_callbacks
  atleast_one_callback = false

  continuous_modified = false
  discrete_modified = false
  saved_in_cb = false
  if !(typeof(continuous_callbacks)<:Tuple{})
    time,upcrossing,event_occurred,event_idx,idx,counter =
              DiffEqBase.find_first_continuous_callback(integrator,continuous_callbacks...)
    if event_occurred
      integrator.event_last_time = idx
      integrator.vector_event_last_time = event_idx
      continuous_modified,saved_in_cb = DiffEqBase.apply_callback!(integrator,continuous_callbacks[idx],time,upcrossing,event_idx)
    else
      integrator.event_last_time = 0
      integrator.vector_event_last_time = 1
    end
  end
  if !integrator.force_stepfail && !(typeof(discrete_callbacks)<:Tuple{})
    discrete_modified,saved_in_cb = DiffEqBase.apply_discrete_callback!(integrator,discrete_callbacks...)
  end
  if !saved_in_cb
    savevalues!(integrator)
  end

  integrator.u_modified = continuous_modified || discrete_modified
  if integrator.u_modified
    integrator.do_error_check = false
    handle_callback_modifiers!(integrator)
  end
end

function handle_callback_modifiers!(integrator::ODEIntegrator)
  integrator.reeval_fsal = true
end

function apply_step!(integrator)

  integrator.accept_step = false # yay we got here, don't need this no more

  #Update uprev
  if alg_extrapolates(integrator.alg)
    if isinplace(integrator.sol.prob)
      recursivecopy!(integrator.uprev2,integrator.uprev)
    else
      integrator.uprev2 = integrator.uprev
    end
  end
  if isinplace(integrator.sol.prob)
    recursivecopy!(integrator.uprev,integrator.u)
    if integrator.alg isa DAEAlgorithm
      recursivecopy!(integrator.duprev, integrator.du)
    end
  else
    integrator.uprev = integrator.u
    if integrator.alg isa DAEAlgorithm
      integrator.duprev = integrator.du
    end
  end

  #Update dt if adaptive or if fixed and the dt is allowed to change
  if integrator.opts.adaptive || integrator.dtchangeable
    integrator.dt = integrator.dtpropose
  elseif integrator.dt != integrator.dtpropose && !integrator.dtchangeable
    error("The current setup does not allow for changing dt.")
  end

  # Update fsal if needed
  if has_discontinuity(integrator) && first_discontinuity(integrator) == integrator.tdir * integrator.t
      handle_discontinuities!(integrator)
      get_current_isfsal(integrator.alg, integrator.cache) && reset_fsal!(integrator)
  elseif all_fsal(integrator.alg, integrator.cache) || get_current_isfsal(integrator.alg, integrator.cache)
    if integrator.reeval_fsal || integrator.u_modified || (typeof(integrator.alg)<:DP8 && !integrator.opts.calck) || (typeof(integrator.alg)<:Union{Rosenbrock23,Rosenbrock32} && !integrator.opts.adaptive)
        reset_fsal!(integrator)
    else # Do not reeval_fsal, instead copyto! over
      if isinplace(integrator.sol.prob)
        recursivecopy!(integrator.fsalfirst,integrator.fsallast)
      else
        integrator.fsalfirst = integrator.fsallast
      end
    end
  end
end

handle_discontinuities!(integrator) = pop_discontinuity!(integrator)

function calc_dt_propose!(integrator,dtnew)
  if (typeof(integrator.alg) <: Union{ROCK2,ROCK4,SERK2,ESERK4,ESERK5}) && integrator.opts.adaptive && (integrator.iter >= 1)
    (integrator.alg isa ROCK2) && (dtnew = min(dtnew,typeof(dtnew)((((min(integrator.alg.max_stages,200)^2.0)*.811 - 1.5)/integrator.eigen_est))))
    (integrator.alg isa ROCK4) && (dtnew = min(dtnew,typeof(dtnew)((((min(integrator.alg.max_stages,152)^2.0)*.353 - 3)/integrator.eigen_est))))
    (integrator.alg isa SERK2) && (dtnew = min(dtnew,typeof(dtnew)((0.8*250*250/(integrator.eigen_est+1.0)))))
    (integrator.alg isa ESERK4) && (dtnew = min(dtnew,typeof(dtnew)((0.98*4000*4000/integrator.eigen_est))))
    (integrator.alg isa ESERK5) && (dtnew = min(dtnew,typeof(dtnew)((0.98*2000*2000/integrator.eigen_est))))
  end
  dtpropose = integrator.tdir*min(abs(integrator.opts.dtmax),abs(dtnew))
  dtpropose = integrator.tdir*max(abs(dtpropose),timedepentdtmin(integrator))
  integrator.dtpropose = dtpropose
end

function fix_dt_at_bounds!(integrator)
  if integrator.tdir > 0
    integrator.dt = min(integrator.opts.dtmax,integrator.dt)
  else
    integrator.dt = max(integrator.opts.dtmax,integrator.dt)
  end
  dtmin = timedepentdtmin(integrator)
  if integrator.tdir > 0
    integrator.dt = max(integrator.dt,dtmin) #abs to fix complex sqrt issue at end
  else
    integrator.dt = min(integrator.dt,dtmin) #abs to fix complex sqrt issue at end
  end
end

function handle_tstop!(integrator)
  if has_tstop(integrator)
    tdir_t = integrator.tdir * integrator.t
    tdir_tstop = first_tstop(integrator)
    if tdir_t == tdir_tstop
      while tdir_t == tdir_tstop #remove all redundant copies
        pop_tstop!(integrator)
        has_tstop(integrator) ? (tdir_tstop = first_tstop(integrator)) : break
      end
      integrator.just_hit_tstop = true
    elseif tdir_t > tdir_tstop
      if !integrator.dtchangeable
        DiffEqBase.change_t_via_interpolation!(integrator, integrator.tdir * pop_tstop!(integrator), Val{true})
        integrator.just_hit_tstop = true
      else
        error("Something went wrong. Integrator stepped past tstops but the algorithm was dtchangeable. Please report this error.")
      end
    end
  end
end

function reset_fsal!(integrator)
  # Under these condtions, these algorithms are not FSAL anymore
  integrator.destats.nf += 1
  if typeof(integrator.cache) <: OrdinaryDiffEqMutableCache ||
     (typeof(integrator.cache) <: CompositeCache &&
      typeof(integrator.cache.caches[1]) <: OrdinaryDiffEqMutableCache)
    integrator.f(integrator.fsalfirst,integrator.u,integrator.p,integrator.t)
  else
    integrator.fsalfirst = integrator.f(integrator.u,integrator.p,integrator.t)
  end
  # Do not set false here so it can be checked in the algorithm
  # integrator.reeval_fsal = false
end

DiffEqBase.nlsolve_f(f, alg::OrdinaryDiffEqAlgorithm) = f isa SplitFunction && issplit(alg) ? f.f1 : f
DiffEqBase.nlsolve_f(f, alg::DAEAlgorithm) = f
DiffEqBase.nlsolve_f(integrator::ODEIntegrator) =
  nlsolve_f(integrator.f, unwrap_alg(integrator, true))

function (integrator::ODEIntegrator)(t,deriv::Type=Val{0};idxs=nothing)
  current_interpolant(t,integrator,idxs,deriv)
end

(integrator::ODEIntegrator)(val::AbstractArray,t::Union{Number,AbstractArray},deriv::Type=Val{0};idxs=nothing) = current_interpolant!(val,t,integrator,idxs,deriv)

# Interface used by DelayDiffEq
has_tstop(integrator) = !isempty(integrator.opts.tstops)
first_tstop(integrator) = first(integrator.opts.tstops)
pop_tstop!(integrator) = pop!(integrator.opts.tstops)

has_discontinuity(integrator) = !isempty(integrator.opts.d_discontinuities)
first_discontinuity(integrator) = first(integrator.opts.d_discontinuities)
pop_discontinuity!(integrator) = pop!(integrator.opts.d_discontinuities)
