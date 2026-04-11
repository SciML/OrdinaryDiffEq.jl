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

## Usage: 2D reaction-diffusion with ADI-style splitting

The canonical AMF use case is a PDE whose discrete Laplacian separates along each spatial dimension — the split reduces the dense Rosenbrock-W `W` solve to a pair of Kronecker-structured solves that can be done much faster. This example is a 2D reaction-diffusion problem,

```
uₜ = A uₓₓ + B u_yy + g(u, t),     g(u, t) = u²(1 − u) + eᵗ,
```

on the unit square with homogeneous Dirichlet boundary conditions. After second-order finite differencing on an `N × N` grid, the Jacobian of the linear (diffusion) part splits as

```
J = A (Lx ⊗ I) + B (I ⊗ Lx) = Jx + Jy
```

where `Lx` is the tridiagonal 1D Laplacian. AMF approximates the stage operator as

```
W ≈ −(I − γ Jx)(I − γ Jy) / γ,
```

so each stage does two independent block-diagonal tridiagonal solves instead of one dense solve — the standard ADI preconditioner, applied automatically by `AMF(...)` whenever you hand it the split.

Install the package alongside a Rosenbrock-W subpackage:

```julia
using Pkg
Pkg.add(["OrdinaryDiffEqAMF", "OrdinaryDiffEqRosenbrock"])
```

Set up the problem:

```julia
using OrdinaryDiffEqAMF
using OrdinaryDiffEqRosenbrock
using SciMLOperators
using LinearAlgebra

function setup_fd2d_problem(; A = 0.1, B = 0.1, N = 40, final_t = 1.0)
    h = 1 / (N + 1)

    # 1D Laplacian with Dirichlet zero BCs, as a SciMLOperator.
    D    = (1 / h^2) * Tridiagonal(ones(N - 1), -2 * ones(N), ones(N - 1))
    D_op = MatrixOperator(D)

    # 2D Jacobian pieces: diffusion in x and in y, via Kronecker products.
    Jx_op = A * Base.kron(D_op, IdentityOperator(N))
    Jy_op = B * Base.kron(IdentityOperator(N), D_op)
    J_op  = cache_operator(Jx_op + Jy_op, zeros(N^2))

    # Nonlinear reaction term.
    g!(du, u, p, t) = (@. du += u^2 * (1 - u) + exp(t); return nothing)

    # Full RHS: linear diffusion + nonlinear reaction.
    function f!(du, u, p, t)
        mul!(du, J_op, u)
        g!(du, u, p, t)
        return nothing
    end

    # Initial condition: a bump that vanishes on the boundary.
    u0 = [16 * (h * i) * (h * j) * (1 - h * i) * (1 - h * j) for i in 1:N for j in 1:N]

    # Hand the split to build_amf_function so Rosenbrock-W sees the
    # structured W_prototype instead of a dense one.
    amf_func = build_amf_function(f!; jac = J_op, split = (Jx_op, Jy_op))
    return ODEProblem(amf_func, u0, (0.0, final_t))
end

prob = setup_fd2d_problem()
sol  = solve(prob, AMF(ROS34PW1a); abstol = 1e-8, reltol = 1e-8)
```

For this 1600-unknown problem the AMF stage solve is substantially faster than the dense equivalent, and the speedup grows with `N` because the dense `W` solve is `O(N⁶)` while the two Kronecker-structured AMF solves are `O(N³)` total. A regression test inside the package (`lib/OrdinaryDiffEqAMF/test/test_fd2d.jl`) measures the speedup on the same problem against a baseline `solve(..., ROS34PW1a())` without AMF.

### Supplying AMF factors directly

Sometimes you want to control the factors yourself — for instance to pre-factorize each 1D operator with something cheaper than a generic structured solver, or to bake in an FFT plan. Pass the custom factors as `amf_factors`:

```julia
using SciMLOperators: ScalarOperator

gamma_op = ScalarOperator(
    1.0;
    update_func       = (old, u, p, t; gamma = 1.0) -> gamma,
    accepted_kwargs   = Val((:gamma,)),
)

# Each factor is (I − γ Jᵢ), pre-assembled as a Kronecker-structured operator.
custom_factors = (
    Base.kron(IdentityOperator(N) - gamma_op * A * D_op, IdentityOperator(N)),
    Base.kron(IdentityOperator(N), IdentityOperator(N) - gamma_op * B * D_op),
)

amf_func = build_amf_function(
    f!;
    jac         = J_op,
    split       = (Jx_op, Jy_op),
    amf_factors = custom_factors,
)

sol = solve(ODEProblem(amf_func, u0, (0.0, 1.0)), AMF(ROS34PW1a); reltol = 1e-8)
```

`split` and `amf_factors` must contain the same number of operators: the split is used to propagate Jacobian updates, and the factors determine the actual `W_prototype` that the linear solver sees.

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
