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
const DEFAULTBETA2S = (
    beta2_default(Tsit5()), beta2_default(Vern7()), beta2_default(Rosenbrock23()),
    beta2_default(Rodas5P()), beta2_default(FBDF()), beta2_default(FBDF()))
const DEFAULTBETA1S = (
    beta1_default(Tsit5(), DEFAULTBETA2S[1]), beta1_default(Vern7(), DEFAULTBETA2S[2]),
    beta1_default(Rosenbrock23(), DEFAULTBETA2S[3]), beta1_default(
        Rodas5P(), DEFAULTBETA2S[4]),
    beta1_default(FBDF(), DEFAULTBETA2S[5]), beta1_default(FBDF(), DEFAULTBETA2S[6]))

callbacks_exists(integrator) = !isempty(integrator.opts.callbacks)
current_nonstiff(current) = ifelse(current <= NUM_NONSTIFF, current, current - NUM_STIFF)

function DefaultODEAlgorithm(; lazy = true, stiffalgfirst = false, kwargs...)
    nonstiff = (Tsit5(), Vern7(lazy = lazy))
    stiff = (Rosenbrock23(; kwargs...), Rodas5P(; kwargs...), FBDF(; kwargs...),
        FBDF(; linsolve = LinearSolve.KrylovJL_GMRES(), kwargs...))
    AutoAlgSwitch(nonstiff, stiff; stiffalgfirst)
end

function isdefaultalg(alg::CompositeAlgorithm{
        <:Any, <:Tuple{Tsit5, Vern7, Rosenbrock23, Rodas5P, FBDF, FBDF}})
    true
end

function DiffEqBase.__init(prob::ODEProblem, ::Nothing, args...; kwargs...)
    DiffEqBase.init(
        prob, DefaultODEAlgorithm(autodiff = AutoFiniteDiff()),
        args...; wrap = Val(false), kwargs...)
end
function DiffEqBase.__solve(prob::ODEProblem, ::Nothing, args...; kwargs...)
    DiffEqBase.solve(
        prob, DefaultODEAlgorithm(autodiff = AutoFiniteDiff()),
        args...; wrap = Val(false), kwargs...)
end

function is_stiff(integrator, alg, ntol, stol, is_stiffalg, current)
    eigen_est, dt = integrator.eigen_est, integrator.dt
    stiffness = abs(eigen_est * dt /
                    STABILITY_SIZES[nonstiffchoice(integrator.opts.reltol)]) # `abs` here is just for safety
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

function stiffchoice(reltol, len, mass_matrix)
    x = if len > MEDIUMSIZE
        DefaultSolverChoice.KrylovFBDF
    elseif len > SMALLSIZE
        DefaultSolverChoice.FBDF
    else
        if reltol < LOW_TOL || mass_matrix != I
            DefaultSolverChoice.Rodas5P
        else
            DefaultSolverChoice.Rosenbrock23
        end
    end
    Int(x)
end

function default_autoswitch(AS::AutoSwitchCache, integrator)
    len = length(integrator.u)
    reltol = integrator.opts.reltol

    # Choose the starting method
    if AS.current == 0
        choice = if AS.stiffalgfirst || integrator.f.mass_matrix != I
            stiffchoice(reltol, len, integrator.f.mass_matrix)
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
    if integrator.f.mass_matrix != I
        #don't change anything
    elseif (!AS.is_stiffalg && AS.count > AS.maxstiffstep)
        integrator.dt = dt * AS.dtfac
        AS.is_stiffalg = true
        AS.current = stiffchoice(reltol, len, integrator.f.mass_matrix)
    elseif (AS.is_stiffalg && AS.count < -AS.maxnonstiffstep)
        integrator.dt = dt / AS.dtfac
        AS.is_stiffalg = false
        AS.current = nonstiffchoice(reltol)
    end
    return AS.current
end

# hack for the default alg
function is_mass_matrix_alg(alg::CompositeAlgorithm{
        <:Any, <:Tuple{Tsit5, Vern7, Rosenbrock23, Rodas5P, FBDF, FBDF}})
    true
end
