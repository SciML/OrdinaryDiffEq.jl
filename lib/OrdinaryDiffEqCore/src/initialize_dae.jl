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

struct OverrideInit{T, F} <: DiffEqBase.DAEInitializationAlgorithm
    abstol::T
    nlsolve::F
end

function OverrideInit(; abstol = 1e-10, nlsolve = nothing)
    OverrideInit(abstol, nlsolve)
end
OverrideInit(abstol) = OverrideInit(; abstol = abstol, nlsolve = nothing)

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
        alg::DefaultInit, x::Val{true})
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(integrator.opts.abstol), x)
    else
        _initialize_dae!(integrator, prob,
            BrownFullBasicInit(integrator.opts.abstol), x)
    end
end

function _initialize_dae!(integrator, prob::ODEProblem,
        alg::DefaultInit, x::Val{false})
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(integrator.opts.abstol), x)
    else
        _initialize_dae!(integrator, prob,
            BrownFullBasicInit(integrator.opts.abstol), x)
    end
end

function _initialize_dae!(integrator, prob::DAEProblem,
        alg::DefaultInit, x::Val{false})
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(integrator.opts.abstol), x)
    elseif prob.differential_vars === nothing
        _initialize_dae!(integrator, prob,
            ShampineCollocationInit(), x)
    else
        _initialize_dae!(integrator, prob,
            BrownFullBasicInit(integrator.opts.abstol), x)
    end
end

function _initialize_dae!(integrator, prob::DAEProblem,
        alg::DefaultInit, x::Val{true})
    if SciMLBase.has_initializeprob(prob.f)
        _initialize_dae!(integrator, prob,
            OverrideInit(integrator.opts.abstol), x)
    elseif prob.differential_vars === nothing
        _initialize_dae!(integrator, prob,
            ShampineCollocationInit(), x)
    else
        _initialize_dae!(integrator, prob,
            BrownFullBasicInit(integrator.opts.abstol), x)
    end
end

## Nonlinear Solver Defaulting

## If an alg is given use it
default_nlsolve(alg, isinplace, u, initprob, autodiff = false) = alg

## If the initialization is trivial just use nothing alg
function default_nlsolve(
        ::Nothing, isinplace, u::Nothing, ::NonlinearProblem, autodiff = false)
    nothing
end

function default_nlsolve(
        ::Nothing, isinplace, u::Nothing, ::NonlinearLeastSquaresProblem, autodiff = false)
    nothing
end

function OrdinaryDiffEqCore.default_nlsolve(::Nothing, isinplace, u, ::NonlinearProblem, autodiff = false)
    error("This ODE requires a DAE initialization and thus a nonlinear solve but no nonlinear solve has been loaded. To solve this problem, do `using OrdinaryDiffEqNonlinearSolve` or pass a custom `nlsolve` choice into the `initializealg`.")
end

function OrdinaryDiffEqCore.default_nlsolve(::Nothing, isinplace, u, ::NonlinearLeastSquaresProblem, autodiff = false)
    error("This ODE requires a DAE initialization and thus a nonlinear solve but no nonlinear solve has been loaded. To solve this problem, do `using OrdinaryDiffEqNonlinearSolve` or pass a custom `nlsolve` choice into the `initializealg`.")
end

## NoInit

function _initialize_dae!(integrator, prob::Union{ODEProblem, DAEProblem},
        alg::NoInit, x::Union{Val{true}, Val{false}})
end

## OverrideInit

function _initialize_dae!(integrator, prob::Union{ODEProblem, DAEProblem},
        alg::OverrideInit, isinplace::Union{Val{true}, Val{false}})
    initializeprob = prob.f.initializeprob

    # If it doesn't have autodiff, assume it comes from symbolic system like ModelingToolkit
    # Since then it's the case of not a DAE but has initializeprob
    # In which case, it should be differentiable
    isAD = if initializeprob.u0 === nothing
        AutoForwardDiff
    elseif has_autodiff(integrator.alg)
        alg_autodiff(integrator.alg) isa AutoForwardDiff
    else
        true
    end

    alg = default_nlsolve(alg.nlsolve, isinplace, initializeprob.u0, initializeprob, isAD)
    nlsol = solve(initializeprob, alg)
    if isinplace === Val{true}()
        integrator.u .= prob.f.initializeprobmap(nlsol)
    elseif isinplace === Val{false}()
        integrator.u = prob.f.initializeprobmap(nlsol)
    else
        error("Unreachable reached. Report this error.")
    end

    if nlsol.retcode != ReturnCode.Success
        integrator.sol = SciMLBase.solution_new_retcode(integrator.sol,
            ReturnCode.InitialFailure)
    end
end
