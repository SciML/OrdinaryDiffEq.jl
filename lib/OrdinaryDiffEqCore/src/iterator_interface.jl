# Hook so SDE can redirect perform_step! to its own module's function.
@inline _perform_step!(integrator) = perform_step!(integrator, integrator.cache)

# Shared implementation — no type annotation so SDE can call it too.
function _step!(integrator)
    if integrator.opts.advance_to_tstop
        while integrator.tdir * integrator.t < first(integrator.opts.tstops)
            loopheader!(integrator)
            (integrator.do_error_check && check_error!(integrator) != ReturnCode.Success) &&
                return integrator.sol.retcode
            _perform_step!(integrator)
            _loopfooter!(integrator)
        end
    else
        loopheader!(integrator)
        (integrator.do_error_check && check_error!(integrator) != ReturnCode.Success) &&
            return integrator.sol.retcode
        _perform_step!(integrator)
        _loopfooter!(integrator)
        while !integrator.accept_step
            loopheader!(integrator)
            (integrator.do_error_check && check_error!(integrator) != ReturnCode.Success) &&
                return integrator.sol.retcode
            _perform_step!(integrator)
            _loopfooter!(integrator)
        end
    end
    handle_tstop!(integrator)
    return integrator.sol.retcode
end

step!(integrator::ODEIntegrator) = _step!(integrator)
