
```@meta
CollapsedDocStrings = true
```
# OrdinaryDiffEqRosenbrock

Methods for small and medium sized stiff systems of differential equations.
At high tolerances, `>1e-2`, `Rosenbrock23` is a good choice.
At medium tolerances `>1e-8` it is recommended you use `Rodas5P` or `Rodas4P`,
the former is more efficient, but the latter is more reliable.
For larger systems look at multistep methods.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqRosenbrock", "Rodas5P")
```

## Full list of solvers

```@docs
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
