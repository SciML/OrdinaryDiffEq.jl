```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqAMF

**Approximate Matrix Factorization (AMF)** is a wrapper around Rosenbrock-W methods that replaces the dense `W = I - γJ` solve at each stage with a product of simpler factors:

```
W ≈ -∏ᵢ (I - γJᵢ) / γ
```

where the `Jᵢ` are individually-structured pieces of the Jacobian (e.g. upper and lower triangular parts, directional operators in dimensional splittings, mode blocks in spectral discretizations). Each factor can be solved independently with a cheap structured solver, so the overall per-step work drops from "one dense solve" to "a handful of structured solves".

## Key Properties

AMF methods provide:

  - **Structured linear solves** — each factor is handled by its own `SciMLOperator`-aware solver, so you never materialize a dense `W`.
  - **Any Rosenbrock-W base method** — `AMF` is an algorithm wrapper, not a fixed integrator; plug in `Rosenbrock23`, `ROS34PW1a`, `Rodas4`, etc.
  - **Works with the existing SciML split/operator infrastructure** — you describe the Jacobian split as `SciMLOperators` and pass them through `build_amf_function`.
  - **Preserves the Rosenbrock-W order** under the usual AMF consistency assumptions.

## When to Use AMF

AMF is recommended when:

  - Your Jacobian has **exploitable structure** that a dense LU solve would throw away — e.g. ADI-style dimensional splits, upper+lower triangular decompositions, block-diagonal plus a low-rank coupling, FFT/diagonalizable pieces.
  - You're using a **Rosenbrock-W method** and the `W` solve is the dominant per-step cost.
  - You can write the structured factors as `SciMLOperators`.

If your Jacobian is genuinely dense and unstructured, plain Rosenbrock with a standard `LinearSolve` algorithm will be simpler and faster — AMF only helps when the product-of-factors approximation is significantly cheaper than the true `W` solve.

## How It Works

`AMF(AlgType)` wraps a Rosenbrock-W algorithm type. When you call `solve(prob, AMF(ROS34PW1a))`, the wrapper:

 1. Instantiates the inner algorithm with `linsolve = SciMLOpFactorization()`, a `LinearSolve` backend that respects `SciMLOperator` products and solves them factor-by-factor rather than concretizing them.
 2. Delegates the rest of the integration to the inner Rosenbrock-W method unchanged.

To use it, the `ODEFunction` must carry a `jac_prototype` (the true Jacobian, as a `SciMLOperator`) *and* a `W_prototype` (the AMF approximation of `-(I - γJ)/γ`, built from the factors). `build_amf_function` assembles both for you.

## Usage

First install the package alongside a Rosenbrock-W subpackage:

```julia
using Pkg
Pkg.add(["OrdinaryDiffEqAMF", "OrdinaryDiffEqRosenbrock"])
```

Minimal example with a 2-way split Jacobian:

```julia
using OrdinaryDiffEqAMF
using OrdinaryDiffEqRosenbrock
using SciMLOperators
using LinearAlgebra

# Problem: du/dt = f(u, t), with Jacobian J = J_upper + J_lower.
function f!(du, u, p, t)
    # ... fill du ...
end

# Per-part update functions populate the structured operators in place.
fjac_upper(J, u, p, t) = (# fill upper-triangular J)
fjac_lower(J, u, p, t) = (# fill lower-triangular J)

N = 20
J1_op = MatrixOperator(UpperTriangular(zeros(N, N)); update_func! = fjac_upper)
J2_op = MatrixOperator(LowerTriangular(zeros(N, N)); update_func! = fjac_lower)
J_op  = cache_operator(J1_op + J2_op, zeros(N))

# Hand the split to build_amf_function so it can build the W_prototype for you.
func = build_amf_function(f!; jac = J_op, split = (J1_op, J2_op))

u0    = zeros(N)
tspan = (0.0, 1.0)
prob  = ODEProblem(func, u0, tspan)

sol = solve(prob, AMF(ROS34PW1a))
```

### Supplying AMF factors directly

Sometimes you want to control the factors yourself — for example when the factors aren't just `I - γJᵢ` but have their own structure (pre-factorized banded matrices, FFT plans, etc.). Pass them as `amf_factors`:

```julia
func = build_amf_function(f!;
    jac = J_op,
    split = (J1_op, J2_op),
    amf_factors = (F1_op, F2_op),
)

sol = solve(prob, AMF(Rosenbrock23))
```

`split` and `amf_factors` must contain the same number of operators; the split is used to propagate Jacobian updates while the factors determine what `W_prototype` actually looks like.

### Forwarding keyword arguments

Any extra kwargs passed to `AMF(alg; kwargs...)` are forwarded to the inner algorithm constructor:

```julia
sol = solve(prob, AMF(ROS34PW1a; autodiff = AutoFiniteDiff()))
```

## Full list of exports

```@docs
OrdinaryDiffEqAMF.AMF
OrdinaryDiffEqAMF.AMFOperator
OrdinaryDiffEqAMF.build_amf_function
```
