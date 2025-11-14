"""
    IDSolve()

First order solver for `ImplicitDiscreteSystems`.
"""
struct IDSolve{NLS} <:
       OrdinaryDiffEqAlgorithm
    nlsolve::NLS
    extrapolant::Symbol
end

function IDSolve(;
        nlsolve = NewtonRaphson(),
        extrapolant = :constant,
)
    IDSolve{typeof(nlsolve)}(nlsolve, extrapolant)
end
