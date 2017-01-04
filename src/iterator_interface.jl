function start(integrator::ODEIntegrator)
  initialize!(integrator,integrator.cache)
  integrator.iter
end

function next(integrator::ODEIntegrator,state)
  integrator.iter += 1
  step(integrator)
  integrator,integrator.iter
end

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

function step(integrator::ODEIntegrator)
  if integrator.opts.advance_to_tstop
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.opts.tstops)
      ode_loopheader!(integrator)
      perform_step!(integrator,integrator.cache)
      ode_loopfooter!(integrator)
    end
  else
    ode_loopheader!(integrator)
    perform_step!(integrator,integrator.cache)
    ode_loopfooter!(integrator)
  end
  if !isempty(integrator.opts.tstops) && integrator.t == top(integrator.opts.tstops)
   pop!(integrator.opts.tstops)
   integrator.just_hit_tstop = true
  end
end
eltype(integrator::ODEIntegrator) = typeof(integrator)
