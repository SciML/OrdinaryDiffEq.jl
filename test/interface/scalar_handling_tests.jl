using OrdinaryDiffEq

# https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/390
solve(ODEProblem((x, p, t) -> -x, 1.0, (0.0, 50.0)), Rosenbrock23(autodiff = false))
