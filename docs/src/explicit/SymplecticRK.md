# OrdinaryDiffEqSymplecticRK

A symplectic integrator is an integrator whose solution resides on a symplectic manifold.
Because of discretization error, when it is solving a Hamiltonian system it doesn't get exactly the correct trajectory on the manifold.
Instead, that trajectory itself is perturbed `O(Δtn)` for the order n from the true trajectory.
Then there's a linear drift due to numerical error of this trajectory over time
Normal integrators tend to have a quadratic (or more) drift, and do not have any good global guarantees about this phase space path (just local).
What means is that symplectic integrators tend to capture the long-time patterns better than normal integrators because of this lack of drift and this almost guarantee of periodicity. 

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqSymplecticRK", "McAte2")
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
