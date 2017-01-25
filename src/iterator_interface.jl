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
    if integrator.opts.verbose
      warn("Interrupted. Larger maxiters is needed.")
    end
    postamble!(integrator)
    return true
  end
  if integrator.opts.unstable_check(integrator.dt,integrator.t,integrator.u)
    if integrator.opts.verbose
      warn("Instability detected. Aborting")
    end
    postamble!(integrator)
    return nothing
  end
  if integrator.dt == zero(integrator.t)
    if integrator.opts.verbose
      warn("dt == 0. Aborting")
    end
    postamble!(integrator)
    return nothing
  end
  if isempty(integrator.opts.tstops)
    postamble!(integrator)
    return true
  elseif integrator.just_hit_tstop
    integrator.just_hit_tstop = false
    if integrator.opts.stop_at_next_tstop
      postamble!(integrator)
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
  handle_tstop!(integrator)
end

eltype(integrator::ODEIntegrator) = typeof(integrator)
