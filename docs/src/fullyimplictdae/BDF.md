```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqBDF

Multistep BDF methods, good for large stiff systems.

```@eval
using OrdinaryDiffEqBDF

function f2(out, du, u, p, t)
    out[1] = -0.04u[1] + 1e4 * u[2] * u[3] - du[1]
    out[2] = +0.04u[1] - 3e7 * u[2]^2 - 1e4 * u[2] * u[3] - du[2]
    out[3] = u[1] + u[2] + u[3] - 1.0
end
u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]
tspan = (0.0, 100000.0)
differential_vars = [true, true, false]
prob = DAEProblem(f2, du₀, u₀, tspan, differential_vars = differential_vars)
sol = solve(prob, DFBDF())
```

## Full list of solvers

### DAE

```@docs
DImplicitEuler
DABDF2
DFBDF
```
