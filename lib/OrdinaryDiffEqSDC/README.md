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
| `node_type` | `:Legendre` | node distribution: `:Legendre`, `:Equidistant` |
| `quad_type` | `:RadauRight` | which endpoints are nodes: `:Gauss` (neither), `:RadauLeft`, `:RadauRight`, `:Lobatto` (both) |
| `num_sweeps` | `3` | number of sweeps `K` |
| `sweeper` | `:BE` | the preconditioner `QΔ`: `:BE`, `:FE`, `:Trapezoid`, `:LU`, `:Picard`, `:BEpar`, `:MIN_SR_NS` |
| `step_update` | `:quadrature` | `:quadrature` (`u_n + Δt Σ_m w_m f_m`) or `:lastnode` (`u_M`) |

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
sol = solve(prob, SDC(num_nodes = 3, num_sweeps = 3); dt = 0.1, adaptive = false)

# 7th order, and a much better preconditioner for stiff problems
sol = solve(
    prob, SDC(num_nodes = 4, num_sweeps = 6, sweeper = :LU);
    dt = 0.1, adaptive = false
)
```

Fixed step size only — pass `adaptive = false` and a `dt`.

## Sweepers

With `δ₁ = τ₁` and `δ_m = τ_m − τ_{m−1}`:

| `sweeper` | `QΔ` | notes |
|---|---|---|
| `:BE` | `QΔ[i,j] = δ_j` for `j ≤ i` | implicit Euler between the nodes; the original Dutt–Greengard–Rokhlin sweep |
| `:FE` | `QΔ[i,j] = δ_{j+1}` for `j < i` | explicit Euler between the nodes; strictly lower triangular, so the sweep needs no solve at all |
| `:Trapezoid` | `(QΔ_BE + QΔ_FE)/2` | second order, so it gains two orders on the first sweep |
| `:LU` | `Uᵀ` where `Qᵀ = LU` | Weiser's LU trick; makes the stiff-limit iteration matrix nilpotent, and is usually the best serial choice |
| `:Picard` | `0` | unpreconditioned Picard iteration; diverges for stiff problems, useful for teaching and testing |
| `:BEpar` | `diag(τ)` | **diagonal**: implicit Euler from the step start to each node |
| `:MIN_SR_NS` | `diag(τ)/M` | **diagonal**: Čaklović et al. 2024, nilpotent non-stiff iteration matrix |

The two diagonal sweepers decouple the sweep across the nodes. They run serially
here; running them in parallel across `M` threads is the next step.

## References

- A. Dutt, L. Greengard and V. Rokhlin, *Spectral deferred correction methods for
  ordinary differential equations*, BIT Numerical Mathematics 40 (2000) 241–266.
- M. Weiser, *Faster SDC convergence on non-equidistant grids by DIRK sweeps*,
  BIT Numerical Mathematics 55 (2015) 1219–1241.
- R. Speck, *Parallelizing spectral deferred corrections across the method*,
  Computing and Visualization in Science 19 (2018) 75–83.
- G. Čaklović, T. Lunet, S. Götschel and D. Ruprecht, *Improving efficiency of
  parallel across the method spectral deferred corrections*, arXiv:2403.18641.

The coefficient conventions and the convergence-order test follow
[qmat](https://github.com/Parallel-in-Time/qmat).
