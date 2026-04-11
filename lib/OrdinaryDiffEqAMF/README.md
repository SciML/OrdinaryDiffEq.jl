# OrdinaryDiffEqAMF

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/OrdinaryDiffEq/stable/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

OrdinaryDiffEqAMF is a subpackage of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) monorepo that provides **Approximate Matrix Factorization (AMF)** wrappers for Rosenbrock-W methods. AMF lets you replace the dense `W = I - γJ` solve at each Rosenbrock-W stage with a product of simpler factors, `W ≈ ∏ᵢ (I - γJᵢ)`, which is often much cheaper for problems whose Jacobian has exploitable structure (e.g. dimensionally-split discretizations of PDEs).

## Installation

```julia
import Pkg
Pkg.add("OrdinaryDiffEqAMF")
```

For targeted dependencies, pair it with a Rosenbrock-W subpackage such as `OrdinaryDiffEqRosenbrock`:

```julia
Pkg.add(["OrdinaryDiffEqAMF", "OrdinaryDiffEqRosenbrock"])
```

## How it works

`AMF` is a thin algorithm wrapper around any Rosenbrock-W method (`ROS34PW1a`, `Rosenbrock23`, `Rodas4`, …). When you call `solve(prob, AMF(ROS34PW1a))` it:

1. Configures the inner Rosenbrock-W algorithm with `linsolve = SciMLOpFactorization()`, a linear solver that keeps the `W_prototype` as a `SciMLOperator` product instead of concretizing it to a dense matrix.
2. Expects the `ODEFunction` to carry a structured `jac_prototype` and `W_prototype` made of `SciMLOperators`. Build them with `build_amf_function`.

The result: each stage's `W` solve is done one factor at a time, typically with structured solves (triangular, banded, FFT-diagonalizable, etc.) that are asymptotically cheaper than a dense factorization.

## Usage

```julia
using OrdinaryDiffEqAMF
using OrdinaryDiffEqRosenbrock
using SciMLOperators
using LinearAlgebra

# A 20×20 problem whose Jacobian splits into upper + lower triangular parts.
N = 20

function f!(du, u, p, t)
    # ... fill du ...
end

# Per-part Jacobians, expressed as SciMLOperators.
J1_op = MatrixOperator(UpperTriangular(zeros(N, N)); update_func! = fjac_upper)
J2_op = MatrixOperator(LowerTriangular(zeros(N, N)); update_func! = fjac_lower)
J_op  = cache_operator(J1_op + J2_op, zeros(N))

# Build the ODEFunction with AMF-aware jac_prototype / W_prototype.
func = build_amf_function(f!; jac = J_op, split = (J1_op, J2_op))

u0    = zeros(N)
tspan = (0.0, 1.0)
prob  = ODEProblem(func, u0, tspan)

sol = solve(prob, AMF(ROS34PW1a))
```

If you want to supply custom structured AMF factors directly (rather than letting the wrapper build `I - γJᵢ` for you), pass `amf_factors = (F1, F2, …)` to `build_amf_function`.

Extra keyword arguments given to `AMF(alg; kwargs...)` are forwarded to the inner algorithm constructor, so you can still tune autodiff, step control, etc.:

```julia
sol = solve(prob, AMF(ROS34PW1a; autodiff = AutoFiniteDiff()))
```

## Exports

- `AMF` — algorithm wrapper
- `AMFOperator` — builds the `W_prototype` operator from a Jacobian, a split, or explicit factors
- `build_amf_function` — helper that assembles the `ODEFunction` with the right `jac_prototype`, `W_prototype`, and sparsity pattern
- `SciMLOpFactorization` — the `SciMLOperator`-aware linear solver used internally

## Documentation

See the main [OrdinaryDiffEq.jl documentation](https://docs.sciml.ai/OrdinaryDiffEq/stable/) for background on Rosenbrock-W methods and the broader SciML ecosystem.
