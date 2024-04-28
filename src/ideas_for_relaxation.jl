#= Modification of the perform_step function into three part:
    - computations
    - modif_step
    - finalize_step
    Their role is described below
=#
function perform_stepNEW!(integrator, cache)

    # Variable to know if dt has changed during perform_step
    integrator.dt_has_changed_in_performstep = false

    # computations! will only contain the mathematical scheme
    # i.e the computations of the u(t+dt)
    # the result is store not in integrator.u but integrator.u_propose
    computations!(integrator, cache)

    # modif_step! enables to modify the step like when we want to perform a relaxation
    # for this we give a new struture that can be defined either by us for already known
    # modification we want to do or by a user (see below)
    modif_step!(integrator, cache, integrator.modif)

    # finalize_step! will do staff related to the solver like integrator.stats, register integrator.fsal
    # and register integrator.u
    finalize_step!(integrator, cache)
end


# To not break all the code we can just add a function to know if the algorithms are 
# written with the new form of perform_step! or not.
while integrator.tdir * integrator.t < first(integrator.opts.tstops)
    loopheader!(integrator)
    if integrator.do_error_check && check_error!(integrator) != ReturnCode.Success
        return integrator.sol
    end

    # This is just to not break all the package
    if isInTheNewFrameWork(integrator.alg)
        # The new performstep!
        perform_stepNEW!(integrator, integrator.cache)
    else
        # The already existing perform_step!
        perform_step!(integrator, integrator.cache)
    end

    loopfooter!(integrator)
    if isempty(integrator.opts.tstops)
        break
    end
end


# modif_step! will look like :

modif_step!(integrator, cache, ::Nothing) = nothing

function modif_step!(integrator, cache, modif)
    
    # Perform the modifications
    modif(integrator, cache)

    # Here we check the validity of chaging dt if it has changed
    # if it is valid integrator.changed_valid will be true, if not it will be false
    integrator.changed_valid = true
    if integrator_dt_has_changed_in_performstep
        # check dt in [dtmin, dtmax]
        # things related to tstops
        # surely other things
        if integrator.changed_valid
            integrator.u_propose = integrator.u_changed
            integrator.dt = integrator.dt_changed
        else
            # print error or warning
        end
    else if integrator.dt_changed != zero(integrator.dt_changed)
        # print error
    end
end

######
# The example of a modification will be like the relaxation :
struct Relaxation{OPT, INV}
    opt::OPT
    inv::INV
end

# Here is the modif function for relaxation
function (r::Relaxtion)(integrator)
    # We fix here the bounds of interval where we are going to look for the relaxation
    # and taking accound the bounds [dtmin, dtmax] and the presence of tstops
    γmin = integrator.dtmin / integrator.dt  
    γmax = min(integrator.dtmax / first(integrator.opts.tstops)) / integrator.dt

    S_u = integrator.dt*(integrator.u_propose-integrator.uprev) 
    ## need to define or use a method "minimize" with foure main arguments : 
        # - a method of minimization
        # - the function to minimize
        # - interval where we look for
        # - starting point (fix to 1) for the minimizing method
    γ_opt = minimize(   r.opt, 
                        γ -> norm(r.inv(γ*S_u .+ integrator.uprev) .- r.inv(integrator.uprev)), 
                        (γmin, γmax), 
                        1.0)
    # new dt
    integrator.dt_changed = integrator.dt * γ_opt
    integrator.dt_has_changed = true

    # update u
    integrator.u_changed = integrator.uprev + γ_opt*S_u
end


## example of a code 






