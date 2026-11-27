```@meta
CollapsedDocStrings = true
```

# OrdinaryDiffEqBDF

Backward Differentiation Formula (BDF) methods are multistep implicit methods specifically designed for solving large stiff systems of differential equations. They are the preferred choice for very large systems (>1000 equations) where other implicit methods become computationally expensive.

## Key Properties

BDF methods offer:

  - **Excellent efficiency for large systems** (>1000 ODEs)
  - **L-stable behavior** for orders 1 and 2 only
  - **Adaptive order and stepsize** control for optimal performance
  - **Alpha-stability** for higher orders (but less stable than L-stable methods for problems with large complex eigenvalues)

## When to Use BDF Methods

BDF methods are recommended for:

  - **Large stiff systems** with more than 1000 equations
  - **Very stiff problems** where other implicit methods struggle
  - **Long-time integration** of stiff systems
  - **Parabolic PDEs** after spatial discretization
  - **Reaction-diffusion systems** and chemical kinetics
  - **Circuit simulation** and other engineering applications with large stiff systems

## Solver Selection Guide

### Recommended methods

  - **`QNDF`**: Adaptive order quasi-constant timestep BDF, best general choice for large systems
  - **`FBDF`**: Fixed-leading coefficient BDF, often more efficient than QNDF

## Matrix-free Krylov solves of DAEs

For an in-place `DAEProblem`, `DFBDF` with its default Newton solver can apply the Newton
matrix without forming it by passing a supported Krylov `linsolve` and leaving `concrete_jac`
unset or false, for example `DFBDF(; linsolve = KrylovJL_GMRES())`. The Newton matrix
`∂F/∂u + (α/(γ h)) ∂F/∂du` is then applied as a Jacobian-vector product, and
`stats.njacs` stays zero. The `cj = α/(γ h)` of the current step is used at every stage. The
matrix-free path uses the selected `autodiff` backend for Jacobian-vector products, including
when `DAEFunction.jvp` is supplied.

The supported combinations of Jacobian-vector product and Krylov solver are:

| `autodiff`                                  | `KrylovJL_GMRES` | `KrylovJL_FGMRES` | `KrylovJL_BICGSTAB` |
|:------------------------------------------- |:----------------:|:-----------------:|:-------------------:|
| `AutoForwardDiff()` (the default)           | yes              | yes               | yes                 |
| `AutoFiniteDiff(; fdtype = Val(:central))`  | yes              | yes               | no                  |
| any other finite-difference backend         | no               | no                | no                  |

Forward finite differences destabilize stiff DAEs such as Robertson's through the matrix-free
operator, and BiCGStab fails with the centered stencil on it, so those combinations throw an
`ArgumentError` when the integrator is initialized; use forward AD or a concrete `linsolve`.
The Newton iteration sets the linear tolerances and measures the residual in the left
preconditioner of the solve, so these are also rejected with an `ArgumentError`: other
Krylov.jl methods, `KrylovJL` with `ldiv = false` or the Krylov.jl keywords `M`, `N`, `W`,
`atol` and `rtol`, and any other `linsolve` that does not need a concrete matrix, such as
`SimpleGMRES` and the `IterativeSolversJL*` wrappers. Pass preconditioners through `precs` and tolerances to `solve`.

A Krylov solve takes one scalar relative tolerance. It is the requested `reltol` (also for a
fixed step); a vector `reltol` is reduced to its smallest component, and the result is floored
at working precision (`eps`), so a zero component cannot make the stopping tolerances zero.
For a finite-difference Jacobian product, which is linear only to about `sqrt(eps)`, the
tolerance `r` used below is further floored at `cbrt(eps)`; an AD product keeps `reltol`.
The integration tolerances stay componentwise. A correction is accepted only if its true
left-preconditioned residual `‖Pl \ (b - A u)‖` is at most `r` times `‖Pl \ b‖` for the
current Newton iteration, plus the roundoff of forming it: the right-hand side and product of the
current iterate alone. The allowance taken from the step's first residual only lets the Krylov
solve stop early; it never excuses an unresolved correction. A zero correction is never accepted for
a nonzero right-hand side, nor is a correction whose residual or bound is not finite, or for
which the relative bound `r ‖Pl \ b‖` falls below `floatmin`. Otherwise one cold retry without the allowance decides,
and the step is rejected unless the retry supplies a nonzero correction that clearly reduced the
residual (the same bound with `sqrt(r)` in place of `r`). Every linear solve, a rejected retry included, is counted in
`stats.nsolve`. Because the stopping test is on a residual norm, rows of the residual
scaled orders of magnitude below the others are invisible to it: give such a problem a
preconditioner that equilibrates the rows. On the matrix-free path `resize!` of a
`DFBDF` integrator throws an `ArgumentError`: the state dimension is fixed at initialization.

## Optional FBDF time filtering

`FBDF(time_filter = true)` adds embedded candidates to the fixed-coefficient
nonlinear solve. Once sufficient history is available, base orders 1–4 can
produce an order-`k+1` candidate. Order 3 also supplies a second-order BDF3-Stab
candidate, including when `max_order = Val(3)`. The maximum order bounds the
accepted candidate; an order-raising candidate is eligible during adaptive
stepping only when the order-change waiting period has elapsed.

The filters follow [DeCaria, Guzel, Layton, and Li](https://arxiv.org/abs/1810.06670),
with an adjusted order-raising weight for FBDF's fixed-coefficient equation on
unequal steps. If `e` is that equation's solution error on the monic polynomial
of degree `k+1` with exact history, and `c` is the newest coefficient in its
`(k+1)`st divided difference, the weight is `e / (1 + c*e)`. This cancels the
leading local error and reduces to the paper's weight on equal steps. The
paper's constant-step stability result should not be read as a stability
guarantee for arbitrary sequences of unequal steps.

Adaptive stepping keeps each candidate's state, derivative, and error estimate
together and ranks eligible candidates by the permitted step size. Step-size
control uses the selected candidate's order. An order-raising selection promotes the
base order for the next step. Selecting BDF3-Stab retains the BDF3 solve so its
embedded candidates remain available; repeated failures can reduce the base
order. STALD's unfiltered derivative
estimates are not applied to filtered candidates. With `adaptive = false`, the
highest available candidate within `max_order` is used; the startup order still
requires enough history. Dense output includes the additional history point
needed by an order-raising candidate.

Time filtering currently requires the identity mass matrix. Use
`FBDF(time_filter = false)` for mass-matrix ODEs and DAEs. Filtering is opt-in
and adds right-hand-side evaluations; compare work at matched achieved accuracy
before choosing it for performance. The default remains `time_filter = false`.

## Performance Characteristics

  - **Most efficient for systems with >1000 equations**
  - **Outperform Runge-Kutta methods** on very large stiff systems
  - **Memory efficient** due to multistep structure
  - **Excel at very low tolerances** (1e-9 and below)
  - **Particularly effective** for problems arising from PDE discretizations

## Comparison with Other Methods

Choose BDF methods over:

  - **Rosenbrock methods**: When system size > 1000 equations
  - **SDIRK methods**: For very large stiff systems where RK methods become expensive
  - **Explicit methods**: For any stiff problem

Choose other methods over BDF when:

  - **System size < 100**: Rosenbrock or SDIRK methods often more efficient
  - **Problems with large complex eigenvalues**: Rosenbrock and L-stable SDIRK methods are more stable due to BDF methods only being alpha-stable
  - **Moderate stiffness**: SDIRK methods may be more robust
  - **Non-stiff problems**: Use explicit methods like Tsit5

```@eval
first_steps = evalfile("./common_first_steps.jl")
first_steps("OrdinaryDiffEqBDF", "QNDF")
```

## Full list of solvers

```@docs
ABDF2
QNDF
QNDF1
QNDF2
QBDF
QBDF1
QBDF2
MEBDF2
FBDF
NordsieckBDF
DNordsieckBDF
```
