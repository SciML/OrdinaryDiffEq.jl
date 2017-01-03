function start(integrator::ODEIntegrator)
  initialize!(integrator,integrator.cache)
  integrator.iter
end

function next(integrator::ODEIntegrator,state)
  state += 1
  if integrator.advance_to_tstop
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      perform_step!(integrator,integrator.cache)
      ode_loopfooter!(integrator)
    end
  else
    ode_loopheader!(integrator)
    perform_step!(integrator,integrator.cache)
    ode_loopfooter!(integrator)
  end
  !isempty(integrator.tstops) && integrator.t == top(integrator.tstops) && pop!(integrator.tstops)
  integrator,state
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
  if isempty(integrator.tstops)
    ode_postamble!(integrator)
    return true
  end
  false
end

eltype(integrator::ODEIntegrator) = typeof(integrator)
