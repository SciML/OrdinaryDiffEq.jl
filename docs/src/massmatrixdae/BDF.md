```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqBDF

Multistep BDF methods, good for large stiff systems.

```julia
function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    nothing
end
M = [1.0 0 0
     0 1.0 0
     0 0 0]
f = ODEFunction(rober, mass_matrix = M)
prob_mm = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1e5), (0.04, 3e7, 1e4))
sol = solve(prob_mm, FBDF(), reltol = 1e-8, abstol = 1e-8)
```

## Full list of solvers

```@docs; canonical=false
ABDF2
QNDF
QNDF1
QNDF2
QBDF
QBDF1
QBDF2
MEBDF2
FBDF
```
