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
