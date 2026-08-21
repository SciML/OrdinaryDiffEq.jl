# GlobalDiffEq: Global Error Estimation and Control

Standard adaptive ODE solvers control the *local* error of each step. The
local tolerances only indirectly control the *global* (accumulated) error of
the solution, which can grow arbitrarily large over long integrations or on
unstable problems even when every step satisfies its local tolerance. The
GlobalDiffEq sublibrary provides solvers and solver wrappers that estimate
the global error, and in several cases control it to a requested global
tolerance.

To use these methods:

```julia
using GlobalDiffEq
```

The global-error-estimating solvers [`GLEE24`](@ref), [`GLEE35`](@ref)
(Constantinescu 2016), and the Dormand-Prince-based [`MM5GEE`](@ref) (Makazaga
and Murua 2003) carry a running, asymptotically correct estimate of the
solution's global error at every time point, at the cost of a few extra stages
per step. Solving with them produces a solution whose states are
`ArrayPartition`s (`sol.u[i].x[1]` the solution, `sol.u[i].x[2]` the global
error estimate); [`global_error_estimate`](@ref) extracts the estimates.

[`GlobalErrorTransport`](@ref) wraps any adaptive solver with global error
estimation and `gtol`-based control by integrating the linearized
error-transport equation driven by the dense-output defect;
[`global_error_estimate`](@ref)`(prob, alg)` exposes the standalone estimator.

[`GlobalRichardson`](@ref) wraps any fixed-step method in global Richardson
extrapolation over whole solves, interpreting `abstol` and `reltol` as global
tolerances. It is the most robust and most expensive option.

For example, solving while tracking the global error along the trajectory:

```julia
using GlobalDiffEq

function lorenz!(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    return
end
prob = ODEProblem(lorenz!, [1.0; 0.0; 0.0], (0.0, 10.0))
sol = solve(prob, GLEE35(); abstol = 1.0e-8, reltol = 1.0e-8)
errs = global_error_estimate(sol)  # global error estimate at every sol.t
```

## Global-error-estimating solvers

```@docs
GLEE23
GLEE24
GLEE35
MM5GEE
global_error_estimate
```

## Global error controlling wrappers

```@docs
GlobalRichardson
GlobalErrorTransport
```
