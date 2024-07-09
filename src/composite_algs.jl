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
    isdefault = integrator.alg isa CompositeAlgorithm{
        <:Any, <:Tuple{Tsit5, Vern7, Rosenbrock23, Rodas5P, FBDF, FBDF}}
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

AutoTsit5(alg; kwargs...) = AutoAlgSwitch(Tsit5(), alg; kwargs...)
AutoDP5(alg; kwargs...) = AutoAlgSwitch(DP5(), alg; kwargs...)
