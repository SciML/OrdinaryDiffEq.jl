function step!(integrator::ODEIntegrator)
    if integrator.opts.advance_to_tstop
        while integrator.tdir * integrator.t < first(integrator.opts.tstops)
            loopheader!(integrator)
            (integrator.do_error_check && check_error!(integrator) != ReturnCode.Success) &&
                return integrator.sol.retcode
            perform_step!(integrator, integrator.cache)
            loopfooter!(integrator)
        end
    else
        loopheader!(integrator)
        (integrator.do_error_check && check_error!(integrator) != ReturnCode.Success) &&
            return integrator.sol.retcode
        perform_step!(integrator, integrator.cache)
        loopfooter!(integrator)
        while !integrator.accept_step
            loopheader!(integrator)
            (integrator.do_error_check && check_error!(integrator) != ReturnCode.Success) &&
                return integrator.sol.retcode
            perform_step!(integrator, integrator.cache)
            loopfooter!(integrator)
        end
    end
    handle_tstop!(integrator)
    return integrator.sol.retcode
end
