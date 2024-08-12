# OrdinaryDiffEqSymplecticRK

Solvers for Hamiltonian systems.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqSymplecticRK", "Nystrom4")
```

## Full list of solvers

```@docs
SymplecticEuler
VelocityVerlet
VerletLeapfrog
PseudoVerletLeapfrog
McAte2
Ruth3
McAte3
CandyRoz4
McAte4
CalvoSanz4
McAte42
McAte5
Yoshida6
KahanLi6
McAte8
KahanLi8
SofSpa10
```
