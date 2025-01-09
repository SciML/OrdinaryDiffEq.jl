### AutoSwitch
### Designed to switch between two solvers, stiff and non-stiff

function AutoSwitch(nonstiffalg, stiffalg;
        maxstiffstep = 10, maxnonstiffstep = 3,
        nonstifftol = 9 // 10, stifftol = 9 // 10, dtfac = 2,
        stiffalgfirst = false,
        switch_max = 5)
    AutoSwitch(nonstiffalg, stiffalg, maxstiffstep, maxnonstiffstep,
        promote(nonstifftol, stifftol)..., dtfac, stiffalgfirst, switch_max)
end

#= TODO: possible improvements to stiffness detection:
1. Cost of implicit solver goes up with length(u)^3 (or ~length(u)^2 with sparsity), factor this in
2. Take into account order difference between nonstiff and stiff solver.
    If nonstiff solver is higher order (which is common), the nonstiff solver can be stability constrained,
    but switching to a stiff solver won't help since it will be order constrained to a smaller timestep
3. Track how close dt is to dtmin to avoid aborting even if stiffness not found
4. Track dt for stiff and nonstiff alg separetly and only use the stiff alg if it is maintaining higher dt
=#
function is_stiff(integrator, alg, ntol, stol, is_stiffalg)
    eigen_est, dt = integrator.eigen_est, integrator.dt
    stiffness = abs(eigen_est * dt / alg_stability_size(alg)) # `abs` here is just for safety
    tol = is_stiffalg ? stol : ntol
    os = oneunit(stiffness)
    bool = !(stiffness <= os * tol)

    if !bool
        integrator.alg.choice_function.successive_switches += 1
    else
        integrator.alg.choice_function.successive_switches = 0
    end

    integrator.do_error_check = (integrator.alg.choice_function.successive_switches >
                                 integrator.alg.choice_function.switch_max || !bool) ||
                                is_stiffalg
    bool
end

function default_autoswitch end

function (AS::AutoSwitchCache)(integrator)
    #horrible awful hack
    isdefault = isdefaultalg(integrator.alg)
    if isdefault
        return default_autoswitch(AS, integrator)
    end
    if AS.current == 0
        AS.current = Int(AS.stiffalgfirst) + 1
        return AS.current
    end

    dt = integrator.dt
    # Successive stiffness test positives are counted by a positive integer,
    # and successive stiffness test negatives are counted by a negative integer
    AS.count = is_stiff(integrator, AS.nonstiffalg, AS.nonstifftol, AS.stifftol,
        AS.is_stiffalg) ?
               AS.count < 0 ? 1 : AS.count + 1 :
               AS.count > 0 ? -1 : AS.count - 1
    if (!AS.is_stiffalg && AS.count > AS.maxstiffstep)
        integrator.dt = dt * AS.dtfac
        AS.is_stiffalg = true
    elseif (AS.is_stiffalg && AS.count < -AS.maxnonstiffstep)
        integrator.dt = dt / AS.dtfac
        AS.is_stiffalg = false
    end
    AS.current = Int(AS.is_stiffalg) + 1
    return AS.current
end

function AutoAlgSwitch(
        nonstiffalg::OrdinaryDiffEqAlgorithm, stiffalg::OrdinaryDiffEqAlgorithm; kwargs...)
    AS = AutoSwitch(nonstiffalg, stiffalg; kwargs...)
    CompositeAlgorithm((nonstiffalg, stiffalg), AS)
end

function AutoAlgSwitch(nonstiffalg::Tuple, stiffalg::Tuple; kwargs...)
    AS = AutoSwitch(nonstiffalg, stiffalg; kwargs...)
    CompositeAlgorithm((nonstiffalg..., stiffalg...), AS)
end
