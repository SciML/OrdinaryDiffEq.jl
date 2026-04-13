# OrdinaryDiffEqMultirate

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/OrdinaryDiffEq/stable/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

OrdinaryDiffEqMultirate is a subpackage of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo that provides **multirate explicit integrators** for split ODEs of the form

```
du/dt = f1(u, t) + f2(u, t)
```

where `f1` is the *fast* component and `f2` is the *slow* component (the standard SciML split-ODE convention). The slow rate is frozen over each macro interval while the fast rate is integrated with many substeps, so you pay the cost of the slow term only once per macro step.

## Installation

```julia
import Pkg
Pkg.add("OrdinaryDiffEqMultirate")
```

## Exports

- `MREEF` — **M**ultirate **R**ichardson **E**xtrapolation with **E**uler as the base **F**ast method. An adaptive explicit integrator that stacks multiple Euler-based fast-slow integrations into an Aitken–Neville Richardson-extrapolation table to boost accuracy.

`MREEF` keyword arguments:

| kwarg | default | meaning |
| --- | --- | --- |
| `m` | `4` | number of fast substeps per macro interval |
| `order` | `4` | extrapolation order (number of base solutions in the Richardson table) |
| `seq` | `:harmonic` | subdivision sequence used for extrapolation; `:harmonic` or `:romberg` |

## Usage

```julia
using OrdinaryDiffEqMultirate
using SciMLBase

# A split ODE: f1 is the fast rate, f2 is the slow rate.
f1!(du, u, p, t) = (@. du = -50 * u)    # fast, stiff-ish
f2!(du, u, p, t) = (@. du = -u)         # slow

u0    = [1.0, 2.0, 3.0]
tspan = (0.0, 1.0)

prob = SplitODEProblem(f1!, f2!, u0, tspan)
sol  = solve(prob, MREEF())               # defaults: m=4, order=4, :harmonic

# Tune the substep count / extrapolation order explicitly:
sol = solve(prob, MREEF(m = 8, order = 6, seq = :romberg))
```

## When to reach for it

Multirate methods pay off when the two halves of your split have very different characteristic timescales and the slow term is expensive to evaluate (so you really want to evaluate it as few times as possible while still resolving the fast dynamics). For problems where both halves cost roughly the same per call, a standard explicit RK on the combined right-hand side is usually simpler and faster.

## Documentation

See the main [OrdinaryDiffEq.jl documentation](https://docs.sciml.ai/OrdinaryDiffEq/stable/) for background on split-ODE solvers and the broader SciML ecosystem.
