# OrdinaryDiffEqSDC

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/OrdinaryDiffEq/stable/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

OrdinaryDiffEqSDC is a subpackage of the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)
monorepo providing **Spectral Deferred Correction** (SDC) methods.

SDC solves the collocation problem on a step

```
u_m = u_n + Δt Σ_j Q[m,j] f(u_j, t_n + τ_j Δt),    m = 1 … M
```

not directly, but by `K` sweeps of the preconditioned iteration

```
u^{k+1} − Δt QΔ f^{k+1} = u_n + Δt (Q − QΔ) f^k
```

where `QΔ` is a cheap lower-triangular approximation of the dense quadrature
matrix `Q`. Because `QΔ` is lower triangular each sweep decomposes into `M`
implicit solves of exactly the shape a DIRK stage solves, so one sweep costs
about one `M`-stage DIRK step. Each sweep adds (at least) one order until the
order of the underlying collocation method is reached, so accuracy is set by the
runtime parameters `num_nodes` and `num_sweeps` rather than by choosing a
different tableau.

## Exports

- `SDC` — serial spectral deferred correction, fixed step size.

| keyword | default | meaning |
|---|---|---|
| `num_nodes` | `3` | number of collocation nodes `M` |
| `node_type` | `SDCNodes.Legendre` | node distribution: `Legendre`, `Equidistant` |
| `quad_type` | `SDCQuadrature.RadauRight` | which endpoints are nodes: `Gauss` (neither), `RadauLeft`, `RadauRight`, `Lobatto` (both) |
| `num_sweeps` | `3` | number of sweeps `K` |
| `sweeper` | `SDCSweeper.BE` | the preconditioner `QΔ`, see below |
| `step_update` | `SDCStepUpdate.Quadrature` | `Quadrature` (`u_n + Δt Σ_m w_m f_m`) or `LastNode` (`u_M`) |

Standard implicit-solver keywords (`autodiff`, `concrete_jac`, `linsolve`,
`nlsolve`) are accepted and behave as elsewhere in the library.

Order: `min(K + s, p*)`, where `s = 1` for the quadrature step update (`0` for
`:lastnode`), plus one extra on the first sweep for the second-order
`:Trapezoid` sweeper, and `p*` is the order of the collocation method —
`2M` for Legendre–Gauss, `2M−1` for Legendre–Radau, `2M−2` for Legendre–Lobatto,
`M` (or `M+1` for even `M` on the symmetric rules) for equidistant nodes.

## Usage

```julia
using OrdinaryDiffEqSDC

prob = ODEProblem((du, u, p, t) -> (du[1] = -u[2]; du[2] = u[1]), [1.0, 0.0], (0.0, 2π))

# 4th order: 3 Radau-right nodes, 3 backward-Euler sweeps, quadrature update
sol = solve(prob, SDC(num_nodes = 3, num_sweeps = 3); abstol = 1e-8, reltol = 1e-8)

# 7th order, and a much better preconditioner for stiff problems
sol = solve(
    prob, SDC(num_nodes = 4, num_sweeps = 6, sweeper = SDCSweeper.LU);
    abstol = 1e-10, reltol = 1e-10
)
```

Adaptive. The embedded solution is the step update formed from the previous
sweep, so the estimate costs a few `axpy`s. Because it measures how far the
iteration is from the collocation solution rather than how far that solution is
from the truth, a stiff problem needs enough sweeps for the iteration to
converge before the step size becomes accuracy limited; too few sweeps and the
estimate sets the step size instead. Pass `adaptive = false` with a `dt` for
fixed steps.

## Sweepers

With `δ₁ = τ₁` and `δ_m = τ_m − τ_{m−1}`:

| `sweeper` | `QΔ` | notes |
|---|---|---|
| `BE` | `QΔ[i,j] = δ_j` for `j ≤ i` | implicit Euler between the nodes; the original Dutt–Greengard–Rokhlin sweep |
| `FE` | `QΔ[i,j] = δ_{j+1}` for `j < i` | explicit Euler between the nodes; strictly lower triangular, so the sweep needs no solve at all |
| `Trapezoid` | `(QΔ_BE + QΔ_FE)/2` | second order, so it gains two orders on the first sweep |
| `LU` | `Uᵀ` where `Qᵀ = LU` | Weiser's LU trick; nilpotent stiff-limit iteration matrix, and the best serial choice |
| `Picard` | `0` | unpreconditioned Picard iteration; diverges for stiff problems, useful for teaching and testing |
| `BEpar` | `diag(τ)` | **diagonal**: implicit Euler from the step start to each node |
| `MIN_SR_NS` | `diag(τ)/M` | **diagonal**: nilpotent non-stiff iteration matrix; optimal as `Δt → 0` and not A-stable |
| `MIN_SR_S` | tabulated | **diagonal**: nilpotent stiff-limit iteration matrix; matches `LU` once the iteration has converged |
| `MIN_SR_FLEX` | `diag(τ)/k` on sweep `k` | **diagonal**, changes every sweep, falls back to `MIN_SR_S` past sweep `M` |

The last five are diagonal, so their sweeps decouple across the nodes. Pass
`threading = true` to run the `M` node solves of each sweep concurrently
(`OrdinaryDiffEqCore`'s `BaseThreads()` and `PolyesterThreads()` also work, for
control over the backend); a non-diagonal sweeper with threading is rejected,
since it couples the nodes within a sweep.

```julia
solve(
    prob,
    SDC(num_nodes = 4, num_sweeps = 6, sweeper = SDCSweeper.MIN_SR_S,
        threading = true);
    abstol = 1e-9, reltol = 1e-9
)
```

Threaded and sequential runs produce bitwise identical solutions. Measured on a
stiff 1-D reaction-diffusion problem, `M = 4` nodes on 4 threads:

| unknowns | sequential | threaded | speedup |
|---|---|---|---|
| 64 | 0.0093 s | 0.0072 s | 1.3x |
| 256 | 0.0905 s | 0.0460 s | 2.0x |
| 1024 | 2.42 s | 1.42 s | 1.7x |

The ceiling is `M`, and the fall-off at 1024 unknowns is the memory cost: each
node keeps its own Jacobian and `W`, so the working set grows as `M x (J + W)`
and leaves cache. Two caveats. Solver statistics (`nf`, `nsolve`, `nw`)
undercount when threaded, because the counters they accumulate into are shared
and updated without synchronisation — this affects every threaded solver in the
library, not just this one. And `f` must be safe to call concurrently.

The three `MIN_SR_*` families are from Čaklović, Lunet, Götschel and Ruprecht,
SIAM J. Sci. Comput. 47 (2025) A430-A453. `MIN_SR_S` has no closed form — its entries solve a nilpotency
condition numerically — so they are tabulated for 2 to 10 nodes, and asking for
more raises an error naming the limit.

On Prothero-Robinson with `lambda = -1e4`, `M = 3` and `dt = 0.1`, the diagonal
sweepers separate clearly:

| sweeper | `K = 3` | `K = 4` | `K = 6` | `K = 8` |
|---|---|---|---|---|
| `LU` (serial) | 3.6e-3 | 8.2e-4 | 2.9e-7 | 1.2e-9 |
| `MIN_SR_S` | 2.3e+0 | 3.1e-3 | 2.0e-6 | 1.3e-9 |
| `MIN_SR_FLEX` | 7.9e-2 | 5.7e-2 | 1.8e-4 | 2.6e-7 |
| `BEpar` | 1.4e+1 | 9.9e+0 | 1.5e+0 | 2.7e-1 |

## References

- A. Dutt, L. Greengard and V. Rokhlin, *Spectral deferred correction methods for
  ordinary differential equations*, BIT Numerical Mathematics 40 (2000) 241–266.
- M. Weiser, *Faster SDC convergence on non-equidistant grids by DIRK sweeps*,
  BIT Numerical Mathematics 55 (2015) 1219–1241.
- R. Speck, *Parallelizing spectral deferred corrections across the method*,
  Computing and Visualization in Science 19 (2018) 75–83.
- G. Čaklović, T. Lunet, S. Götschel and D. Ruprecht, *Improving efficiency of
  parallel across the method spectral deferred corrections*, SIAM J. Sci. Comput. 47 (2025) A430-A453.

The coefficient conventions and the convergence-order test follow
[qmat](https://github.com/Parallel-in-Time/qmat).
