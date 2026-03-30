# Developing A New Problem

New problems should be defined for new types of differential equations, new partial
differential equations, and special subclasses of differential equations for which
solvers can dispatch on for better performance.

To develop a new problem, you need to make a new `DEProblem` and a new `DESolution`.
These types belong in DiffEqBase and should be exported.
The `DEProblem` type should hold all of the mathematical information about the
problem (including all of the meshing information in both space and time),
and the `DESolution` should hold all of the information for the solution.
Then all that is required is to define a `__solve(::DEProblem,alg;kwargs)`
which takes in the problem and returns a solution.

Then to check that the algorithm works, add a dispatch for `test_convergence`
which makes a `ConvergenceSimulation` type. This type already has a plot recipe, so
plotting functionality will already be embedded. This requires that your
problem can take in a true solution, and has a field `errors` which is a
dictionary of symbols for the different error estimates (L2,L infinity, etc.)

After these steps, update the documentation to include the new problem types and
the new associated solvers.
