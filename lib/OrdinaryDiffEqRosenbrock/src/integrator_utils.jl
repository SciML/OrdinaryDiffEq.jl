function apply_step!(integrator)
    update_uprev!(integrator)

    #Update dt if adaptive or if fixed and the dt is allowed to change
    if integrator.opts.adaptive || integrator.dtchangeable
        integrator.dt = integrator.dtpropose
    elseif integrator.dt != integrator.dtpropose && !integrator.dtchangeable
        error("The current setup does not allow for changing dt.")
    end

    # Update fsal if needed
    if has_discontinuity(integrator) &&
       first_discontinuity(integrator) == integrator.tdir * integrator.t
        handle_discontinuities!(integrator)
        get_current_isfsal(integrator.alg, integrator.cache) && reset_fsal!(integrator)
    elseif all_fsal(integrator.alg, integrator.cache) ||
           get_current_isfsal(integrator.alg, integrator.cache)
        if integrator.reeval_fsal || integrator.u_modified ||
           (integrator.alg isa DP8 && !integrator.opts.calck) ||
           (integrator.alg isa Union{Rosenbrock23, Rosenbrock32} &&
            !integrator.opts.adaptive)
            reset_fsal!(integrator)
        else # Do not reeval_fsal, instead copyto! over
            if isinplace(integrator.sol.prob)
                recursivecopy!(integrator.fsalfirst, integrator.fsallast)
            else
                integrator.fsalfirst = integrator.fsallast
            end
        end
    end
    return nothing
end
