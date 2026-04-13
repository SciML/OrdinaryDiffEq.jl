# StochasticDiffEqWeak.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

StochasticDiffEqWeak.jl is a component of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo providing weak-order stochastic Runge-Kutta methods for computing expectations of SDE solutions. While completely independent and usable on its own, users wanting all SDE solvers should use [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

## Solvers

### Roessler SRK Methods (weak order 2)
- `DRI1`, `DRI1NM` - Debrabant-Roessler methods for Ito SDEs
- `RI1`, `RI3`, `RI5`, `RI6` - Roessler weak order 2 methods for Ito SDEs
- `RDI1WM`, `RDI2WM`, `RDI3WM`, `RDI4WM` - Roessler-Debrabant methods

### Stratonovich Weak Methods
- `RS1`, `RS2` - Roessler weak order 2 methods for Stratonovich SDEs

### Platen Methods
- `PL1WM`, `PL1WMA` - Platen's weak order 2 method

### Komori-Burrage Methods
- `NON`, `NON2` - Komori non-commutative noise methods
- `COM` - Komori commutative noise method

### Semi-Implicit/Semi-Explicit
- `SIEA`, `SIEB` - Semi-implicit Euler type A/B
- `SMEA`, `SMEB` - Semi-implicit Milstein type A/B

### Tang-Xiao Methods
- `W2Ito1` - Weak order 2 Ito method

### Implicit Weak Methods
- `IRI1` - Implicit Roessler weak order 2 (drift-implicit)

## Example

```julia
using StochasticDiffEqWeak
f(u, p, t) = 1.01 * u
g(u, p, t) = 0.87 * u
u0 = 0.5
tspan = (0.0, 1.0)
prob = SDEProblem(f, g, u0, tspan)
sol = solve(prob, DRI1())
```

## Available Solvers

For the full list of available SDE solvers, refer to the [DifferentialEquations.jl SDE Solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/sde_solve/) documentation.
