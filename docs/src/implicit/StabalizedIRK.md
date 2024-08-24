```@meta
CollapsedDocStrings = true
```
# OrdinaryDiffEqStabalizedIRK

Stabilized Explicit Runge-Kutta methods,
like Runge-Kutta-Chebyshev methods and ROCK methods
are explicit methods which chain together many stages in a specific way to get large stability regions.
they are made in such a way to converge to a large stability region,
and thus suitable to stiff equations.
However, they converge to having a large stability region in the direction of the negative real axis,
and thus are only stable on a subset of stiff equations which are not dominated by large complex eigenvalues in the Jacobian.

Stabilized implicit methods try to mitigate this problem by being an IMEX type scheme,
requiring a SplitODEProblem where the splitting is designed to treat the large complex eigenvalues implicitly
while treating the large real eigenvalues using a fast explicit stabilized RK type of method.

These methods utilize an upper bound on the spectral radius of the Jacobian.
Users can supply an upper bound by specifying the keyword argument `eigen_est`, for example

```julia
`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`
```

## Installation

To be able to access the solvers in `OrdinaryDiffEqStabalizedIRK`, you must first install them use the Julia package manager:

```julia
using Pkg
Pkg.add("OrdinaryDiffEqStabalizedIRK")
```

This will only install the solvers listed at the bottom of this page.
If you want to explore other solvers for your problem,
you will need to install some of the other libraries listed in the navigation bar on the left.

## Example usage

```julia
using OrdinaryDiffEqStabilizedIRK
A = randn(20, 20)
B = randn(20, 20)
f1 = (u, p, t) -> A * u
f2 = (u, p, t) -> B * u
u0 = randn(20, 1)
tspan = (0.0, 1.0)
prob = SplitODEProblem(f1, f2, u0, tspan)
sol = solve(prob, IRKC())
```

## Full list of solvers

```@docs
IRKC
```
