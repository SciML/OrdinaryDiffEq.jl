@inline function step!(integrator::ODEIntegrator)
  if integrator.opts.advance_to_tstop
    @inbounds while integrator.tdir*integrator.t < integrator.tdir*top(integrator.opts.tstops)
      loopheader!(integrator)
      check_error!(integrator) != :Success && return
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
    end
  else
    @inbounds loopheader!(integrator)
    check_error!(integrator) != :Success && return
    @inbounds perform_step!(integrator,integrator.cache)
    @inbounds loopfooter!(integrator)
    @inbounds while !integrator.accept_step
      loopheader!(integrator)
      perform_step!(integrator,integrator.cache)
      loopfooter!(integrator)
    end
  end
  @inbounds handle_tstop!(integrator)
  nothing
end
