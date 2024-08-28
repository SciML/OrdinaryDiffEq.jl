```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqRosenbrock

Methods for small and medium sized stiff systems of differential equations.
At high tolerances, `>1e-2`, `Rosenbrock23` is a good choice.
At medium tolerances `>1e-8` it is recommended you use `Rodas5P` or `Rodas4P`,
the former is more efficient, but the latter is more reliable.
For larger systems look at multistep methods.

## Example usage

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
sol = solve(prob_mm, Rodas5(), reltol = 1e-8, abstol = 1e-8)
```

## Full list of solvers

```@docs; canonical=false
Rosenbrock23
Rosenbrock32
ROS3P
Rodas3
Rodas23W
Rodas3P
Rodas4
Rodas42
Rodas4P
Rodas4P2
Rodas5
Rodas5P
Rodas5Pe
Rodas5Pr
RosenbrockW6S4OS
ROS2
ROS2PR
ROS2S
ROS3
ROS3PR
Scholz4_7
ROS34PW1a
ROS34PW1b
ROS34PW2
ROS34PW3
ROS34PRw
ROS3PRL
ROS3PRL2
ROK4a
RosShamp4
Veldd4
Velds4
GRK4T
GRK4A
Ros4LStab
```
