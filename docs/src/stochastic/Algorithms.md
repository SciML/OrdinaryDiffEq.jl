# Stochastic solver API

This page lists the stochastic algorithms available through the
OrdinaryDiffEq umbrella package. Choose an algorithm according to the SDE
interpretation, noise structure, stiffness, and required strong or weak order.

## High-order stochastic Runge-Kutta methods

```@docs
StochasticDiffEqHighOrder.RosslerSRA
StochasticDiffEqHighOrder.RosslerSRI
StochasticDiffEqHighOrder.SOSRA
StochasticDiffEqHighOrder.SOSRA2
StochasticDiffEqHighOrder.SOSRI
StochasticDiffEqHighOrder.SOSRI2
StochasticDiffEqHighOrder.SRA
StochasticDiffEqHighOrder.SRA1
StochasticDiffEqHighOrder.SRA2
StochasticDiffEqHighOrder.SRA3
StochasticDiffEqHighOrder.SRI
StochasticDiffEqHighOrder.SRIW1
StochasticDiffEqHighOrder.SRIW2
```

## Iterated-integral-free methods

```@docs
StochasticDiffEqIIF.IIF1M
StochasticDiffEqIIF.IIF1Mil
StochasticDiffEqIIF.IIF2M
```

## Implicit stochastic methods

```@docs
StochasticDiffEqImplicit.ISSEM
StochasticDiffEqImplicit.ISSEulerHeun
StochasticDiffEqImplicit.ImplicitEM
StochasticDiffEqImplicit.ImplicitEulerHeun
StochasticDiffEqImplicit.ImplicitRKMil
StochasticDiffEqImplicit.SImplicitMidpoint
StochasticDiffEqImplicit.SKenCarp
StochasticDiffEqImplicit.STrapezoid
```

## Tau-leaping methods

```@docs
StochasticDiffEqLeaping.CaoTauLeaping
StochasticDiffEqLeaping.ImplicitTauLeaping
StochasticDiffEqLeaping.TauLeaping
StochasticDiffEqLeaping.ThetaTrapezoidalTauLeaping
```

## Low-order stochastic methods

```@docs
StochasticDiffEqLowOrder.EM
StochasticDiffEqLowOrder.EulerHeun
StochasticDiffEqLowOrder.LambaEM
StochasticDiffEqLowOrder.LambaEulerHeun
StochasticDiffEqLowOrder.PCEuler
StochasticDiffEqLowOrder.RKMil
StochasticDiffEqLowOrder.RKMilCommute
StochasticDiffEqLowOrder.SimplifiedEM
StochasticDiffEqLowOrder.SplitEM
```

## Milstein methods

```@docs
StochasticDiffEqMilstein.RKMilGeneral
StochasticDiffEqMilstein.WangLi3SMil_A
StochasticDiffEqMilstein.WangLi3SMil_B
StochasticDiffEqMilstein.WangLi3SMil_C
StochasticDiffEqMilstein.WangLi3SMil_D
StochasticDiffEqMilstein.WangLi3SMil_E
StochasticDiffEqMilstein.WangLi3SMil_F
```

## Stabilized stochastic methods

```@docs
StochasticDiffEqROCK.KomBurSROCK2
StochasticDiffEqROCK.SKSROCK
StochasticDiffEqROCK.SROCK1
StochasticDiffEqROCK.SROCK2
StochasticDiffEqROCK.SROCKC2
StochasticDiffEqROCK.SROCKEM
StochasticDiffEqROCK.TangXiaoSROCK2
```

## Random ordinary differential equation methods

```@docs
StochasticDiffEqRODE.BAOAB
StochasticDiffEqRODE.RandomEM
StochasticDiffEqRODE.RandomHeun
StochasticDiffEqRODE.RandomTamedEM
```

## Mass-action tau leaping

`TauLeaping`, `CaoTauLeaping`, `ImplicitTauLeaping`, and
`ThetaTrapezoidalTauLeaping` accept a `JumpProblem` built from a `DiscreteProblem`,
`PureLeaping()`, and a `MassActionJump`. These solvers retain the mass-action
representation, use its stored rate constants for error control, and apply its
stoichiometry directly. The implicit methods evaluate the mass-action drift
without an intermediate propensity vector during nonlinear iteration.

```@example massaction_leaping
using JumpProcesses, StochasticDiffEqLeaping

jump = MassActionJump([0.1], [[1 => 1]], [[1 => -1, 2 => 1]])
prob = JumpProblem(DiscreteProblem([1000.0, 0.0], (0.0, 1.0)), PureLeaping(), jump)
sol = solve(prob, ImplicitTauLeaping(); dt = 0.01, adaptive = false)
```

For more general propensity functions and count-based updates, these methods
also accept `RegularJump`. Combining it with a `MassActionJump` in the same
leaping problem is not supported. These StochasticDiffEq solvers do not provide
`EnsembleGPUKernel` implementations; see JumpProcesses for its supported GPU
leaping algorithms.

`CaoTauLeaping` currently lacks its adaptive step-selection calculation. Use
`dt` with `adaptive = false` for that method; use `TauLeaping` when adaptive
steps are needed.
