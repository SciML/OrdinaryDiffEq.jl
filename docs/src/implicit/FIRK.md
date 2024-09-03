```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqFIRK

FIRK methods are fully implicit Runge-Kutta methods.
They can have special properties, like be symplectic integrators, and can achieve higher order for the same number of stage in comparison to diagonal methods.
However, the fully implicit methods have a larger implicit system to solve and thus have a higher linear algebra cost.
This can be useful in some contexts to promote more parallelism,
but also since the size of the factorization is cubic and the dominant cost for large equations,
multiplying `O(n^3)` operations to `O((sn)^3)` can be a considerable cost increase for FIRK tableaus,
where `s`, the number of stages, is particularly large.
That said, the restriction to diagonal implicitness imposes order restrictions,
such as SDIRK methods having a maximum order of 5, which can restrict the problems best suited for SDIRK methods.

The most common FIRK method in production are those based on RadauIIA tableaus,
which is an ODE representation of Gaussian collocation.
Like Gaussian collocation, it achieves higher order convergence than its stages, namely order 2s+1 for s stages.
Thus RadauIIA FIRK methods tend to be some of the highest order methods (excluding extrapolation methods).
This means that high order RadauIIA methods are recommended in the same scenarios that high-order explicit Runge-Kutta methods are recommended simply with the restriction of being a stiff equation.
Such scenarios include cases like very low tolerances: RadauIIA methods can be the best performing methods for scenarios where tolerances are `1e-9` and below.
Additionally, for ODE systems of size less than 200, the increased size of the Jacobian is mitigated by improved multithreading,
since BLAS implementations are only good at multithreading LU factorizations after a certain matrix size.
For this reason, RadauIIA methods tend to be recommended in cases where ODE size is small to intermediate and very accurate solutions are required.

They should be tested against the parallel implicit extrapolation which also specialize in this regime.

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqFIRK", "RadauIIA5")
```

## Full list of solvers

```@docs
RadauIIA3
RadauIIA5
RadauIIA9
```
