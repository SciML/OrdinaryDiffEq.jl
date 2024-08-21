# OrdinaryDiffEqRKN

Second order solvers.

To be able to access the solvers in `OrdinaryDiffEqRKN`, you must first install them use the Julia package manager:

```julia
using Pkg
Pkg.add("OrdinaryDiffEqRKN")
```
This will only install the solvers listed at the bottom of this page.
If you want to explore other solvers for your problem,
you will need to install some of the other libraries listed in the navigation bar on the left.

## Example usage
```julia
using OrdinaryDiffEqOrdinaryDiffEqRKN
function HH_acceleration!(dv, v, u, p, t)
    x, y = u
    dx, dy = dv
    dv[1] = -x - 2x * y
    dv[2] = y^2 - y - x^2
end
initial_positions = [0.0, 0.1]
initial_velocities = [0.5, 0.0]
tspan = (0.0, 1.0)
prob = SecondOrderODEProblem(HH_acceleration!, initial_velocities, initial_positions, tspan)
sol = solve(prob, Nystrom4(), dt = 1 / 10)
```

## Full list of solvers

```@docs
IRKN3
IRKN4
Nystrom4
Nystrom4VelocityIndependent
Nystrom5VelocityIndependent
FineRKN4
FineRKN5
DPRKN4
DPRKN5
DPRKN6
DPRKN6FM
DPRKN8
DPRKN12
ERKN4
ERKN5
ERKN7
RKN4
```
