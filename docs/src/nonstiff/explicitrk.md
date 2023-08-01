# Explicit Runge-Kutta Methods

With the help of [FastBroadcast.jl](https://github.com/YingboMa/FastBroadcast.jl), 
we can use threaded parallelism to reduce compute time for all of the explicit Runge-Kutta methods!
The `thread` option determines whether internal broadcasting on appropriate CPU arrays should be serial
(`thread = OrdinaryDiffEq.False()`, default) or use multiple threads 
(`thread = OrdinaryDiffEq.True()`) when Julia is started with multiple threads. 
When we call `solve(prob, alg(thread=OrdinaryDiffEq.True()))`,
we can turn on the multithreading option to achieve acceleration
(for sufficiently large problems).


## Standard Explicit Runge-Kutta Methods

```@docs
Heun
Ralston
Midpoint
RK4
RKM
MSRK5
MSRK6
Anas5
RKO65
OwrenZen3
OwrenZen4
OwrenZen5
BS3
DP5
Tsit5
DP8
TanYam7
TsitPap8
Feagin10
Feagin12
Feagin14
FRK65
PFRK87
Stepanov5
SIR54
Alshina2
Alshina3
Alshina6
```

## Lazy Interpolation Explicit Runge-Kutta Methods

```@docs
BS5
Vern6
Vern7
Vern8
Vern9
```

## Fixed Timestep Only Explicit Runge-Kutta Methods

```@docs
Euler
RK46NL
ORK256
```

## Parallel Explicit Runge-Kutta Methods

```@docs
KuttaPRK2p5
```
