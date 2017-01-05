initialize{uType}(integrator,cache::OrdinaryDiffEqCache,::Type{uType}) =
                error("This algorithm does not have an initialization function")

@inline function loopheader!(integrator)
  # Apply right after iterators / callbacks

  # Accept or reject the step
  if integrator.iter > 0
    if (integrator.opts.adaptive && integrator.accept_step) || !integrator.opts.adaptive
      apply_step!(integrator)
    elseif integrator.opts.adaptive && !integrator.accept_step
      integrator.dt = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
    end
  end

  integrator.iter += 1
  modify_dt_for_tstops!(integrator)
  choose_algorithm!(integrator,integrator.cache)
end

@inline function modify_dt_for_tstops!(integrator)
  tstops = integrator.opts.tstops
  if !isempty(tstops)
    if integrator.opts.adaptive
      if integrator.tdir > 0
        integrator.dt = min(abs(integrator.dt),abs(top(tstops)-integrator.t)) # step! to the end
      else
        integrator.dt = -min(abs(integrator.dt),abs(top(tstops)-integrator.t))
      end
    elseif integrator.dtcache == zero(integrator.t) && integrator.dtchangeable # Use integrator.opts.tstops
      integrator.dt = integrator.tdir*abs(top(tstops)-integrator.t)
    elseif integrator.dtchangeable # always try to step! with dtcache, but lower if a tstops
      integrator.dt = integrator.tdir*min(abs(integrator.dtcache),abs(top(tstops)-integrator.t)) # step! to the end
    end
  end
end

@def ode_exit_conditions begin
  if integrator.iter > integrator.opts.maxiters
    warn("Interrupted. Larger maxiters is needed.")
    postamble!(integrator)
    return nothing
  end
  if integrator.dt == zero(integrator.t)
    warn("dt == 0. Aborting")
    postamble!(integrator)
    return nothing
  end
  if any(isnan,integrator.uprev)
    warn("NaNs detected. Aborting")
    postamble!(integrator)
    return nothing
  end
end

@inline function savevalues!(integrator::ODEIntegrator)
  while !isempty(integrator.opts.saveat) && integrator.tdir*top(integrator.opts.saveat) <= integrator.tdir*integrator.t # Perform saveat
    integrator.saveiter += 1
    curt = pop!(integrator.opts.saveat)
    if integrator.opts.saveat!=integrator.t # If <t, interpolate
      ode_addsteps!(integrator)
      Θ = (curt - integrator.tprev)/integrator.dt
      val = ode_interpolant(Θ,integrator)
      copyat_or_push!(integrator.sol.t,integrator.saveiter,curt)
      copyat_or_push!(integrator.sol.u,integrator.saveiter,val)
      if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    else # ==t, just save
      copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
      if integrator.opts.dense
        integrator.saveiter_dense += 1
        copyat_or_push!(integrator.notsaveat_idxs,integrator.saveiter_dense,integrator.saveiter)
        copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
      end
      if typeof(alg) <: OrdinaryDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    end
  end
  if integrator.opts.save_timeseries && integrator.iter%integrator.opts.timeseries_steps==0
    integrator.saveiter += 1
    copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if integrator.opts.dense
      integrator.saveiter_dense += 1
      copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
      copyat_or_push!(integrator.notsaveat_idxs,integrator.saveiter_dense,integrator.saveiter)
    end
    if typeof(integrator.alg) <: OrdinaryDiffEqCompositeAlgorithm
      copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
    end
  end
  if isspecialdense(integrator.alg)
    resize!(integrator.k,integrator.kshortsize)
  end
end

@inline function postamble!(integrator)

  if integrator.opts.save_timeseries==false
    solution_endpoint_match_cur_integrator!(integrator)
  end

  resize!(integrator.sol.t,integrator.saveiter)
  resize!(integrator.sol.u,integrator.saveiter)
  resize!(integrator.sol.k,integrator.saveiter_dense)

  if integrator.sol.t[end] !=  integrator.t
    error("Solution endpoint doesn't match the current time in the postamble. This should never happen.")
  end
  !(typeof(integrator.prog)<:Void) && Juno.done(integrator.prog)
end

@inline function solution_endpoint_match_cur_integrator!(integrator)
  if integrator.sol.t[integrator.saveiter] !=  integrator.t
    integrator.saveiter += 1
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    if integrator.opts.dense
      integrator.saveiter_dense +=1
      copyat_or_push!(integrator.sol.k,integrator.saveiter_dense,integrator.k)
      copyat_or_push!(integrator.notsaveat_idxs,integrator.saveiter_dense,integrator.saveiter)
    end
  end
end

@inline function loopfooter!(integrator)
  if integrator.opts.adaptive
    integrator.q11 = integrator.EEst^integrator.opts.beta1
    q = integrator.q11/(integrator.qold^integrator.opts.beta2)
    q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),q/integrator.opts.gamma))
    dtnew = integrator.dt/q
    ttmp = integrator.t + integrator.dt
    integrator.accept_step = (!integrator.opts.isoutofdomain(ttmp,integrator.u) && integrator.EEst <= 1.0)
    if integrator.accept_step # Accept
      integrator.t = ttmp
      calc_dt_propose!(integrator,dtnew)
      if !(typeof(integrator.opts.callback)<:Void)
        for c in integrator.opts.callback
          apply_callback!(integrator,c)
        end
      else
        savevalues!(integrator)
      end
    end
  else #Not adaptive
    integrator.t += integrator.dt
    integrator.accept_step = true
    integrator.dtpropose = integrator.dt
    if !(typeof(integrator.opts.callback)<:Void)
      for c in integrator.opts.callback
        apply_callback!(integrator,c)
      end
    else
      savevalues!(integrator)
    end
  end
  if !(typeof(integrator.prog)<:Void) && integrator.opts.progress && integrator.iter%integrator.opts.progress_steps==0
    Juno.msg(integrator.prog,integrator.opts.progress_message(integrator.dt,integrator.t,integrator.u))
    Juno.progress(integrator.prog,integrator.t/integrator.sol.prob.tspan[2])
  end
end

@inline function apply_step!(integrator)

  integrator.accept_step = false # yay we got here, don't need this no more

  #Update uprev
  if typeof(integrator.u) <: AbstractArray
    recursivecopy!(integrator.uprev,integrator.u)
  else
    integrator.uprev = integrator.u
  end

  #Update dt if adaptive or if fixed and the dt is allowed to change
  if integrator.opts.adaptive || integrator.dtchangeable
    integrator.dt = integrator.dt_mod*integrator.dtpropose
  elseif integrator.dt != integrator.dt_mod*integrator.dtpropose && !integrator.dtchangeable
    error("The current setup does not allow for changing dt.")
  end

  # Update fsal if needed
  if isfsal(integrator.alg)
    if integrator.reeval_fsal || (typeof(integrator.alg)<:DP8 && !integrator.opts.calck) || (typeof(integrator.alg)<:Union{Rosenbrock23,Rosenbrock32} && !integrator.opts.adaptive)
      # Under these condtions, these algorithms are not FSAL anymore
      if typeof(integrator.cache) <: OrdinaryDiffEqMutableCache
        integrator.f(integrator.t,integrator.u,integrator.fsalfirst)
      else
        integrator.fsalfirst = integrator.f(integrator.t,integrator.u)
      end
      integrator.reeval_fsal = false
    else
      if typeof(integrator.fsalfirst) <: AbstractArray
        recursivecopy!(integrator.fsalfirst,integrator.fsallast)
      else
        integrator.fsalfirst = integrator.fsallast
      end
    end
  end

  # Update kprev
  integrator.tprev = integrator.t
  #=
  if integrator.calcprevs # Is this a micro-optimization that can be removed?
    if !isspecialdense(integrator.alg) && integrator.opts.calck
      if typeof(integrator.k) <: AbstractArray
        recursivecopy!(integrator.kprev,integrator.k)
      else
        integrator.kprev = integrator.k
      end
    end
  end
  =#

  integrator.dt_mod = typeof(integrator.t)(1)
end



@inline function calc_dt_propose!(integrator,dtnew)
  integrator.qold = max(integrator.EEst,integrator.opts.qoldinit)
  if integrator.tdir > 0
    integrator.dtpropose = min(integrator.opts.dtmax,dtnew)
  else
    integrator.dtpropose = max(integrator.opts.dtmax,dtnew)
  end
  if integrator.tdir > 0
    integrator.dtpropose = max(integrator.dtpropose,integrator.opts.dtmin) #abs to fix complex sqrt issue at end
  else
    integrator.dtpropose = min(integrator.dtpropose,integrator.opts.dtmin) #abs to fix complex sqrt issue at end
  end
end

@inline function handle_tstop!(integrator)
  tstops = integrator.opts.tstops
  if !isempty(tstops)
    t = integrator.t
    ts_top = top(tstops)
    if t == ts_top
      pop!(tstops)
      integrator.just_hit_tstop = true
    elseif t > ts_top
      if !integrator.dtchangeable
        change_t_via_interpolation!(integrator, pop!(tstops), Val{true})
        integrator.just_hit_tstop = true
      else
        error("Something went wrong. Integrator stepped past tstops but the algorithm was dtchangeable. Please report this error.")
      end
    end
  end
end

(integrator::ODEIntegrator)(t) = current_interpolant(t,integrator)
