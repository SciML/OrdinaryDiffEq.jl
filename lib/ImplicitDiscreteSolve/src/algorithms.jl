"""
    IDSolve()

First order solver for `ImplicitDiscreteSystems`.
"""
# struct IDSolve{CS, AD, NLS, FDT, ST, CJ} <:
struct IDSolve{NLS} <:
       OrdinaryDiffEqAlgorithm
    nlsolve::NLS
    extrapolant::Symbol
    controller::Symbol
end

function IDSolve(;
        nlsolve = NewtonRaphson(), #NLNewton(),
        extrapolant = :constant,
        controller = :PI,
    )

    IDSolve{typeof(nlsolve)}(nlsolve, extrapolant, controller)
end
