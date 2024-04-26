# New structure

abstract type PerformStepUpdate end

struct Relaxation{OPT, INV} <: PerformStepUpdate
    opt::OPT
    inv::INV
end

function (r::Relaxtion)(integrator, integrator.cache)

    γmin = integrator.dtmin / integrator.dt  
    γmax = integrator.dtmax / integrator.dt 

    ## need to define or use a method find root with three main arguments : 
        # - a method to minimize
        # - the function to minimize
        # - interval where we look for
    γ_opt = minimize(r.opt, γ -> norm(γ*integrator.dt*(integrator.u_propose-integrator.uprev)) + integrator.uprev, [γmin, γmax])

    integrator.dt_changed = integrator.dt * γ_opt
    integrator_dt_has_changed_in_performstep = true
end



# In the solve function
while integrator.tdir * integrator.t < first(integrator.opts.tstops)
    loopheader!(integrator)
    if integrator.do_error_check && check_error!(integrator) != ReturnCode.Success
        return integrator.sol
    end
    if isInTheNewSyntax(alg)
        perform_step2!(integrator, integrator.cache)
    else
        perform_step!(integrator, integrator.cache)
    end
    loopfooter!(integrator)
    if isempty(integrator.opts.tstops)
        break
    end
end

# Change the perform_step function
function perform_step2!(integrator, integrator.cache, method)
    integrator_dt_has_changed_in_performstep = false
    
    # Computations will really only contain the mathematical scheme
    computations!(integrator, integrator.cache)

    user_update!(integrator, integrator.cache, method)

    if integrator_dt_has_changed_in_performstep
        #check dt in [dtmin, dtmax] and no conflic with tstop
    end

    # Recording will contain all update we need to do into the integrator 
    recording!(integrator)
end

function user_update!(integrator, integrator.cache, method::PerformStepUpdate)


end


