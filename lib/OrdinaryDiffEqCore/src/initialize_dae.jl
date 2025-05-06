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
    elseif !applicable(_initialize_dae!, integrator, prob,
        BrownFullBasicInit(integrator.opts.abstol), x)
        error("`OrdinaryDiffEqNonlinearSolve` is not loaded, which is required for the default initialization algorithm (`BrownFullBasicInit` or `ShampineCollocationInit`). To solve this problem, either do `using OrdinaryDiffEqNonlinearSolve` or pass `initializealg = CheckInit()` to the `solve` function. This second option requires consistent `u0`.")
    else
        _initialize_dae!(integrator, prob,
            BrownFullBasicInit(integrator.opts.abstol), x)
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
            BrownFullBasicInit(integrator.opts.abstol), x)
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
function _initialize_dae!(integrator, prob::AbstractDEProblem, alg::CheckInit,
        isinplace::Union{Val{true}, Val{false}})
    SciMLBase.get_initial_values(
        prob, integrator, prob.f, alg, isinplace; abstol = integrator.opts.abstol)
end
