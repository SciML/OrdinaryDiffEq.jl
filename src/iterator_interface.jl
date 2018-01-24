function start(integrator::ODEIntegrator)
  0
end

@inline function next(integrator::ODEIntegrator,state)
  state += 1
  step!(integrator) # Iter updated in the step! header
  # Next is callbacks -> iterator  -> top
  integrator,state
end

done(integrator::ODEIntegrator) = done(integrator,integrator.iter)

@inline function done(integrator::ODEIntegrator,state)
  if integrator.iter > integrator.opts.maxiters
    if integrator.opts.verbose
      warn("Interrupted. Larger maxiters is needed.")
    end
    postamble!(integrator)
    return true
  end
  if integrator.opts.unstable_check(integrator.dt,integrator.u,integrator.p,integrator.t)
    if integrator.opts.verbose
      warn("Instability detected. Aborting")
    end
    postamble!(integrator)
    return true
  end
  if !integrator.opts.force_dtmin && integrator.opts.adaptive && abs(integrator.dt) <= abs(integrator.opts.dtmin)
    if integrator.opts.verbose
      warn("dt <= dtmin. Aborting. If you would like to force continuation with dt=dtmin, set force_dtmin=true")
    end
    postamble!(integrator)
    return true
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

@inline function step!(integrator::ODEIntegrator)
  if integrator.opts.advance_to_tstop
    @inbounds while integrator.tdir*integrator.t < integrator.tdir*top(integrator.opts.tstops)
      loopheader!(integrator)
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
    end
  else
    @inbounds loopheader!(integrator)
    @inbounds perform_step!(integrator,integrator.cache)
    @inbounds loopfooter!(integrator)
    @inbounds while !integrator.accept_step
      loopheader!(integrator)
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
    end
  end
  @inbounds handle_tstop!(integrator)
end

eltype(integrator::ODEIntegrator) = typeof(integrator)
