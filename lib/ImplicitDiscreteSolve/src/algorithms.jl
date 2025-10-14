"""
    IDSolve()

First order solver for `ImplicitDiscreteSystems`.
"""
struct IDSolve{NLS} <:
       OrdinaryDiffEqAlgorithm
    nlsolve::NLS
    extrapolant::Symbol
    controller::Symbol
end

function IDSolve(;
        nlsolve = NewtonRaphson(),
        extrapolant = :constant,
        controller = :PI
)
    IDSolve{typeof(nlsolve)}(nlsolve, extrapolant, controller)
end
