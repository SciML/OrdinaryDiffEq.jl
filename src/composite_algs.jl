### AutoSwitch
### Designed to switch between two solvers, stiff and non-stiff

function AutoSwitch(nonstiffalg, stiffalg, algtrait = nothing;
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
    bool = stiffness > os * tol

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

function (AS::AutoSwitchCache)(integrator)
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

function AutoAlgSwitch(nonstiffalg::OrdinaryDiffEqAlgorithm, stiffalg::OrdinaryDiffEqAlgorithm, algtrait = nothing; kwargs...)
    AS = AutoSwitch(nonstiffalg, stiffalg, algtrait; kwargs...)
    CompositeAlgorithm((nonstiffalg, stiffalg), AS)
end

function AutoAlgSwitch(nonstiffalg::Tuple, stiffalg::Tuple, algtrait; kwargs...)
    AS = AutoSwitch(nonstiffalg, stiffalg, algtrait; kwargs...)
    CompositeAlgorithm((nonstiffalg..., stiffalg...), AS)
end

AutoTsit5(alg; kwargs...) = AutoAlgSwitch(Tsit5(), alg; kwargs...)
AutoDP5(alg; kwargs...) = AutoAlgSwitch(DP5(), alg; kwargs...)
AutoVern6(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern6(lazy = lazy), alg; kwargs...)
AutoVern7(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern7(lazy = lazy), alg; kwargs...)
AutoVern8(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern8(lazy = lazy), alg; kwargs...)
AutoVern9(alg; lazy = true, kwargs...) = AutoAlgSwitch(Vern9(lazy = lazy), alg; kwargs...)

### Default ODE Solver

EnumX.@enumx DefaultSolverChoice begin
    Tsit5 = 1
    Vern7 = 2
    Rosenbrock23 = 3
    Rodas5P = 4
    FBDF = 5
    KrylovFBDF = 6
end

const NUM_NONSTIFF = 2
const NUM_STIFF = 4
const LOW_TOL = 1e-6
const MED_TOL = 1e-2
const EXTREME_TOL = 1e-9
const SMALLSIZE = 50
const MEDIUMSIZE = 500
const STABILITY_SIZES = (alg_stability_size(Tsit5()), alg_stability_size(Vern7()))
const DEFAULTBETA2S = (beta2_default(Tsit5()), beta2_default(Vern7()), beta2_default(Rosenbrock23()), beta2_default(Rodas5P()), beta2_default(FBDF()), beta2_default(FBDF()))
const DEFAULTBETA1S = (beta1_default(Tsit5(),DEFAULTBETA2S[1]), beta1_default(Vern7(),DEFAULTBETA2S[2]),
                      beta1_default(Rosenbrock23(), DEFAULTBETA2S[3]), beta1_default(Rodas5P(), DEFAULTBETA2S[4]),
                      beta1_default(FBDF(), DEFAULTBETA2S[5]), beta1_default(FBDF(), DEFAULTBETA2S[6]))

callbacks_exists(integrator) = !isempty(integrator.opts.callbacks)
current_nonstiff(current) = ifelse(current <= NUM_NONSTIFF,current,current-NUM_STIFF)

function DefaultODEAlgorithm(; lazy = true, stiffalgfirst = false, kwargs...)
    nonstiff = (Tsit5(), Vern7(lazy = lazy))
    stiff = (Rosenbrock23(;kwargs...), Rodas5P(;kwargs...), FBDF(;kwargs...), FBDF(;linsolve = LinearSolve.KrylovJL_GMRES()))
    AutoAlgSwitch(nonstiff, stiff, DefaultODESolver(); stiffalgfirst)
end

function is_stiff(integrator, alg, ntol, stol, is_stiffalg, current)
    eigen_est, dt = integrator.eigen_est, integrator.dt
    stiffness = abs(eigen_est * dt / STABILITY_SIZES[nonstiffchoice(integrator.opts.reltol)]) # `abs` here is just for safety
    tol = is_stiffalg ? stol : ntol
    os = oneunit(stiffness)
    bool = stiffness > os * tol

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

function nonstiffchoice(reltol)
    x = if reltol < LOW_TOL
        DefaultSolverChoice.Vern7
    else
        DefaultSolverChoice.Tsit5
    end
    Int(x)
end

function stiffchoice(reltol, len)
    x = if len > MEDIUMSIZE
        DefaultSolverChoice.KrylovFBDF
    elseif len > SMALLSIZE
        DefaultSolverChoice.FBDF
    else
        if reltol < LOW_TOL
            DefaultSolverChoice.Rodas5P
        else
            DefaultSolverChoice.Rosenbrock23
        end
    end
    Int(x)
end

function (AS::AutoSwitchCache{DefaultODESolver})(integrator)

    len = length(integrator.u)
    reltol = integrator.opts.reltol

    # Chooose the starting method
    if AS.current == 0
        choice = if AS.stiffalgfirst || integrator.f.mass_matrix != I
            stiffchoice(reltol, len)
        else
            nonstiffchoice(reltol)
        end
        AS.current = choice
        return AS.current
    end

    dt = integrator.dt
    # Successive stiffness test positives are counted by a positive integer,
    # and successive stiffness test negatives are counted by a negative integer
    AS.count = is_stiff(integrator, AS.nonstiffalg, AS.nonstifftol, AS.stifftol,
        AS.is_stiffalg, AS.current) ?
               AS.count < 0 ? 1 : AS.count + 1 :
               AS.count > 0 ? -1 : AS.count - 1
    if (!AS.is_stiffalg && AS.count > AS.maxstiffstep)
        integrator.dt = dt * AS.dtfac
        AS.is_stiffalg = true
        AS.current = stiffchoice(reltol, len)
    elseif (AS.is_stiffalg && AS.count < -AS.maxnonstiffstep)
        integrator.dt = dt / AS.dtfac
        AS.is_stiffalg = false
        AS.current = nonstiffchoice(reltol)
    end
    return AS.current
end
