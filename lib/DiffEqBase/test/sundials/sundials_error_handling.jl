using OrdinaryDiffEq, Sundials, Test

f_oop(u, p, t) = 2u
u0 = 0.5
tspan = (0.0, 1.0)

# Sundials CVODE_Adams does not support complex initial conditions
prob = ODEProblem{false}(f_oop, 1.0 + im, tspan)
@test_throws SciMLBase.ComplexSupportError solve(prob, CVODE_Adams())

# DFBDF is a DAE solver, not compatible with standard ODEProblem
@test_throws SciMLBase.ProblemSolverPairingError solve(prob, DFBDF())
