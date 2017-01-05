function start(integrator::ODEIntegrator)
  0
end

function next(integrator::ODEIntegrator,state)
  state += 1
  step!(integrator) # Iter updated in the step! header
  # Next is callbacks -> iterator  -> top
  integrator,state
end

done(integrator::ODEIntegrator) = done(integrator,integrator.iter)

function done(integrator::ODEIntegrator,state)
  if integrator.iter > integrator.opts.maxiters
    warn("Interrupted. Larger maxiters is needed.")
    ode_postamble!(integrator)
    return true
  end
  if any(isnan,integrator.uprev)
    warn("NaNs detected. Aborting")
    ode_postamble!(integrator)
    return true
  end
  if isempty(integrator.opts.tstops)
    ode_postamble!(integrator)
    return true
  elseif integrator.just_hit_tstop
    integrator.just_hit_tstop = false
    if integrator.opts.stop_at_next_tstop
      ode_postamble!(integrator)
      return true
    end
  end
  false
end

function step!(integrator::ODEIntegrator)
  if integrator.opts.advance_to_tstop
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.opts.tstops)
      loopheader!(integrator)
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
    end
  else
    loopheader!(integrator)
    perform_step!(integrator,integrator.cache)
    loopfooter!(integrator)
    while !integrator.accept_step
      loopheader!(integrator)
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
    end
  end
  if !isempty(integrator.opts.tstops) && integrator.t == top(integrator.opts.tstops)
   pop!(integrator.opts.tstops)
   integrator.just_hit_tstop = true
  elseif !isempty(integrator.opts.tstops) && integrator.t > top(integrator.opts.tstops) && !integrator.dtchangeable
    change_t_via_interpolation!(integrator,pop!(integrator.opts.tstops),Val{false})
    integrator.just_hit_tstop = true
  elseif !isempty(integrator.opts.tstops) && integrator.t > top(integrator.opts.tstops)
    error("Something went wrong. Integrator stepped past tstops but the algorithm was dtchangeable. Please report this error.")
  end
end

eltype(integrator::ODEIntegrator) = typeof(integrator)
