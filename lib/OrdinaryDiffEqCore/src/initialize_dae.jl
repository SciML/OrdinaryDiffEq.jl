struct DefaultInit <: DiffEqBase.DAEInitializationAlgorithm end

struct ShampineCollocationInit{T, F} <: DiffEqBase.DAEInitializationAlgorithm
    initdt::T
    nlsolve::F
end
function ShampineCollocationInit(; initdt = nothing, nlsolve = nothing)
    ShampineCollocationInit(initdt, nlsolve)
end
function ShampineCollocationInit(initdt)
    ShampineCollocationInit(; initdt = initdt, nlsolve = nothing)
end

struct BrownFullBasicInit{T, F} <: DiffEqBase.DAEInitializationAlgorithm
    abstol::T
    nlsolve::F
end
function BrownFullBasicInit(; abstol = 1e-10, nlsolve = nothing)
    BrownFullBasicInit(abstol, nlsolve)
end
BrownFullBasicInit(abstol) = BrownFullBasicInit(; abstol = abstol, nlsolve = nothing)

## Notes

#=
differential_vars = [any(!iszero,x) for x in eachcol(M)]

A column should be zero for an algebraic variable, since that means that the
derivative term doesn't show up in any equations (i.e. is an algebraic variable).
The rows are not necessarily non-zero, for example a flux condition between two
differential variables. But if it's a condition that doesn't involve the algebraic
variable, then the system is not Index 1!

=#

## Expansion

function DiffEqBase.initialize_dae!(integrator::ODEIntegrator,
        initializealg = integrator.initializealg)
    if !hasmethod(_initialize_dae!, (typeof(integrator),
        typeof(integrator.sol.prob), BrownFullBasicInit, 
        Val{DiffEqBase.isinplace(integrator.sol.prob)}))
            error("This ODE requires a DAE initialization and thus a nonlinear solve but no nonlinear solve has been loaded. To solve this problem, do `using OrdinaryDiffEqNonlinearSolve` or pass a custom `nlsolve` choice into the `initializealg`.")
    end
    _initialize_dae!(integrator, integrator.sol.prob,
        initializealg,
        Val(DiffEqBase.isinplace(integrator.sol.prob)))
end

## Default algorithms

function _initialize_dae!(integrator, prob::ODEProblem,
        alg::DefaultInit, x::Union{Val{true}, Val{false}})
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(integrator.opts.abstol), x)
    else
        _initialize_dae!(integrator, prob,
            CheckInit(), x)
    end
end

function _initialize_dae!(integrator, prob::ODEProblem,
        alg::DefaultInit, x::Val{false})
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(integrator.opts.abstol), x)
    else
        _initialize_dae!(integrator, prob,
            CheckInit(), x)
    end
end

function _initialize_dae!(integrator, prob::DAEProblem,
        alg::DefaultInit, x::Union{Val{true}, Val{false}})
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(integrator.opts.abstol), x)
    elseif !applicable(_initialize_dae!, integrator, prob,
        BrownFullBasicInit(), x) &&
           !applicable(_initialize_dae!,
        integrator, prob, ShampineCollocationInit(), x)
        error("`OrdinaryDiffEqNonlinearSolve` is not loaded, which is required for the default initialization algorithm (`BrownFullBasicInit` or `ShampineCollocationInit`). To solve this problem, either do `using OrdinaryDiffEqNonlinearSolve` or pass `initializealg = CheckInit()` to the `solve` function. This second option requires consistent `u0`.")
    elseif prob.differential_vars === nothing
        _initialize_dae!(integrator, prob,
            ShampineCollocationInit(), x)
    else
        _initialize_dae!(integrator, prob,
            CheckInit(), x)
    end
end

function _initialize_dae!(integrator, prob::DiscreteProblem,
        alg::DefaultInit, x::Union{Val{true}, Val{false}})
    if SciMLBase.has_initializeprob(prob.f)
        # integrator.opts.abstol is `false` for `DiscreteProblem`.
        _initialize_dae!(integrator, prob, OverrideInit(one(eltype(prob.u0)) * 1e-12), x)
    end
end

## Nonlinear Solver Defaulting

## If an alg is given use it
default_nlsolve(alg, isinplace, u, initprob, autodiff = false) = alg

## If the initialization is trivial just use nothing alg
function default_nlsolve(
        ::Nothing, isinplace::Val{true}, u::Nothing, ::AbstractNonlinearProblem, autodiff = false)
    nothing
end

function default_nlsolve(
        ::Nothing, isinplace::Val{true}, u::Nothing, ::NonlinearLeastSquaresProblem, autodiff = false)
    nothing
end

function default_nlsolve(
        ::Nothing, isinplace::Val{false}, u::Nothing, ::AbstractNonlinearProblem, autodiff = false)
    nothing
end

function default_nlsolve(
        ::Nothing, isinplace::Val{false}, u::Nothing,
        ::NonlinearLeastSquaresProblem, autodiff = false)
    nothing
end

function OrdinaryDiffEqCore.default_nlsolve(
        ::Nothing, isinplace, u, ::AbstractNonlinearProblem, autodiff = false)
    error("This ODE requires a DAE initialization and thus a nonlinear solve but no nonlinear solve has been loaded. To solve this problem, do `using OrdinaryDiffEqNonlinearSolve` or pass a custom `nlsolve` choice into the `initializealg`.")
end

function OrdinaryDiffEqCore.default_nlsolve(
        ::Nothing, isinplace, u, ::NonlinearLeastSquaresProblem, autodiff = false)
    error("This ODE requires a DAE initialization and thus a nonlinear solve but no nonlinear solve has been loaded. To solve this problem, do `using OrdinaryDiffEqNonlinearSolve` or pass a custom `nlsolve` choice into the `initializealg`.")
end

## NoInit

function _initialize_dae!(integrator, prob::AbstractDEProblem,
        alg::NoInit, x::Union{Val{true}, Val{false}})
end

## OverrideInit

function _initialize_dae!(integrator, prob::AbstractDEProblem,
        alg::OverrideInit, isinplace::Union{Val{true}, Val{false}})
    initializeprob = prob.f.initialization_data.initializeprob

    # If it doesn't have autodiff, assume it comes from symbolic system like ModelingToolkit
    # Since then it's the case of not a DAE but has initializeprob
    # In which case, it should be differentiable
    iu0 = state_values(initializeprob)
    isAD = if iu0 === nothing
        AutoForwardDiff
    elseif has_autodiff(integrator.alg)
        alg_autodiff(integrator.alg) isa AutoForwardDiff
    else
        true
    end

    nlsolve_alg = default_nlsolve(alg.nlsolve, isinplace, iu0, initializeprob, isAD)

    u0, p, success = SciMLBase.get_initial_values(
        prob, integrator, prob.f, alg, isinplace; nlsolve_alg,
        abstol = integrator.opts.abstol, reltol = integrator.opts.reltol)

    if isinplace === Val{true}()
        integrator.u .= u0
    elseif isinplace === Val{false}()
        integrator.u = u0
    else
        error("Unreachable reached. Report this error.")
    end
    integrator.p = p
    sol = integrator.sol
    @reset sol.prob.p = integrator.p
    integrator.sol = sol

    if !success
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol,
            ReturnCode.InitialFailure)
    end
end

## CheckInit
struct CheckInitFailureError <: Exception
    normresid::Any
    abstol::Any
end

function Base.showerror(io::IO, e::CheckInitFailureError)
    print(io,
        "DAE initialization failed: your u0 did not satisfy the initialization requirements, 
        normresid = $(e.normresid) > abstol = $(e.abstol). If you wish for the system to 
        automatically change the algebraic variables to satisfy the algebraic constraints, 
        please pass `initializealg = BrownBasicInit()` to solve (this option will require 
        `using OrdinaryDiffEqNonlinearSolve`). If you wish to perform an initialization on the
        complete u0, please pass initializealg = ShampineCollocationInit() to solve. Note that 
        initialization can be a very difficult process for DAEs and in many cases can be 
        numerically intractable without symbolic manipulation of the system. For an automated 
        system that will generate numerically stable initializations, see ModelingToolkit.jl 
        structural simplification for more details."
    )
end

function _initialize_dae!(integrator, prob::ODEProblem, alg::CheckInit,
        isinplace::Val{true})
    @unpack p, t, f = integrator
    M = integrator.f.mass_matrix
    tmp = first(get_tmp_cache(integrator))
    u0 = integrator.u

    algebraic_vars = [all(iszero, x) for x in eachcol(M)]
    algebraic_eqs = [all(iszero, x) for x in eachrow(M)]
    (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
    update_coefficients!(M, u0, p, t)
    f(tmp, u0, p, t)
    tmp .= ArrayInterface.restructure(tmp, algebraic_eqs .* _vec(tmp))

    normresid = integrator.opts.internalnorm(tmp, t)
    if normresid > integrator.opts.abstol
        throw(CheckInitFailureError(normresid, integrator.opts.abstol))
    end
end

function _initialize_dae!(integrator, prob::ODEProblem, alg::CheckInit,
        isinplace::Val{false})
    @unpack p, t, f = integrator
    u0 = integrator.u
    M = integrator.f.mass_matrix

    algebraic_vars = [all(iszero, x) for x in eachcol(M)]
    algebraic_eqs = [all(iszero, x) for x in eachrow(M)]
    (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
    update_coefficients!(M, u0, p, t)
    du = f(u0, p, t)
    resid = _vec(du)[algebraic_eqs]

    normresid = integrator.opts.internalnorm(resid, t)
    if normresid > integrator.opts.abstol
        throw(CheckInitFailureError(normresid, integrator.opts.abstol))
    end
end

function _initialize_dae!(integrator, prob::DAEProblem,
        alg::CheckInit, isinplace::Val{true})
    @unpack p, t, f = integrator
    u0 = integrator.u
    resid = get_tmp_cache(integrator)[2]

    f(resid, integrator.du, u0, p, t)
    normresid = integrator.opts.internalnorm(resid, t)
    if normresid > integrator.opts.abstol
        throw(CheckInitFailureError(normresid, integrator.opts.abstol))
    end
end

function _initialize_dae!(integrator, prob::DAEProblem,
        alg::CheckInit, isinplace::Val{false})
    @unpack p, t, f = integrator
    u0 = integrator.u

    nlequation_oop = u -> begin
        f((u - u0) / dt, u, p, t)
    end

    nlequation = (u, _) -> nlequation_oop(u)

    resid = f(integrator.du, u0, p, t)
    normresid = integrator.opts.internalnorm(resid, t)
    if normresid > integrator.opts.abstol
        throw(CheckInitFailureError(normresid, integrator.opts.abstol))
    end
function _initialize_dae!(integrator, prob::AbstractDEProblem, alg::CheckInit,
        isinplace::Union{Val{true}, Val{false}})
    SciMLBase.get_initial_values(
        prob, integrator, prob.f, alg, isinplace; abstol = integrator.opts.abstol)
end
