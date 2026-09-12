@doc generic_solver_docstring(
    "Spectral Deferred Correction.

Solves the collocation problem on `num_nodes` quadrature nodes by preconditioned
iteration. On each step the collocation system

```
u_m = u_n + Δt Σ_j Q[m,j] f(u_j, t_n + τ_j Δt),   m = 1, …, M
```

is not solved directly. Instead `num_sweeps` sweeps of

```
u^{k+1} - Δt QΔ f^{k+1} = u_n + Δt (Q - QΔ) f^k
```

are taken, where `QΔ` is a cheap lower-triangular approximation of the dense
quadrature matrix `Q` chosen by `sweeper`. Because `QΔ` is lower triangular the
sweep decomposes into `M` successive implicit solves of exactly the shape a DIRK
stage solves, so each sweep costs about as much as one `M`-stage DIRK step.

Each sweep raises the order of the method by (at least) one until the order of
the underlying collocation method is reached, so the accuracy is tuned by
`num_sweeps` and `num_nodes` rather than by picking a different tableau.

Adaptive. The embedded estimate is the difference between the step updates
formed from the last two sweeps, which costs a handful of `axpy`s because both
iterates are already in the cache.",
    "SDC",
    "Spectral Deferred Correction method.",
    """@article{dutt2000spectral,
    title={Spectral deferred correction methods for ordinary differential equations},
    author={Dutt, Alok and Greengard, Leslie and Rokhlin, Vladimir},
    journal={BIT Numerical Mathematics},
    volume={40},
    number={2},
    pages={241--266},
    year={2000},
    publisher={Springer}}
    @article{weiser2015faster,
    title={Faster {SDC} convergence on non-equidistant grids by {DIRK} sweeps},
    author={Weiser, Martin},
    journal={BIT Numerical Mathematics},
    volume={55},
    number={4},
    pages={1219--1241},
    year={2015},
    publisher={Springer}}""",
    """
    - `num_nodes`: number of collocation nodes `M` per step.
    - `node_type`: node distribution, `SDCNodes.Legendre` or `.Equidistant`.
    - `quad_type`: which interval endpoints are nodes, `SDCQuadrature.Gauss`
        (neither), `.RadauLeft` (left), `.RadauRight` (right) or `.Lobatto` (both).
    - `num_sweeps`: number of sweeps `K` per step.
    - `sweeper`: the preconditioner `QΔ`. `SDCSweeper.BE` (implicit Euler between
        the nodes), `.FE` (explicit Euler between the nodes), `.Trapezoid`, `.LU`
        (Weiser's LU trick), `.Picard` (`QΔ = 0`), `.BEpar` (diagonal, implicit
        Euler from the step start to each node) or `.MIN_SR_NS` (diagonal,
        `diag(τ)/M`).
    - `step_update`: how the step solution is formed from the node values,
        `SDCStepUpdate.Quadrature` (`u_n + Δt Σ_m w_m f_m`) or `.LastNode`
        (`u_M`, which requires the last node to be the right endpoint).
    """,
    """
    num_nodes = 3,
    node_type = SDCNodes.Legendre,
    quad_type = SDCQuadrature.RadauRight,
    num_sweeps = 3,
    sweeper = SDCSweeper.BE,
    step_update = SDCStepUpdate.Quadrature,
    """
)
struct SDC{AD, F, F2, CJ} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm
    num_nodes::Int
    node_type::SDCNodes.T
    quad_type::SDCQuadrature.T
    num_sweeps::Int
    sweeper::SDCSweeper.T
    step_update::SDCStepUpdate.T
    linsolve::F
    nlsolve::F2
    autodiff::AD
    concrete_jac::CJ
end

function SDC(;
        num_nodes::Int = 3,
        node_type::SDCNodes.T = SDCNodes.Legendre,
        quad_type::SDCQuadrature.T = SDCQuadrature.RadauRight,
        num_sweeps::Int = 3,
        sweeper::SDCSweeper.T = SDCSweeper.BE,
        step_update::SDCStepUpdate.T = SDCStepUpdate.Quadrature,
        autodiff = AutoForwardDiff(), concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton()
    )
    sdc_validate(num_nodes, quad_type, num_sweeps, step_update)
    autodiff = _fixup_ad(autodiff)
    return SDC(
        num_nodes, node_type, quad_type, num_sweeps, sweeper, step_update,
        linsolve, nlsolve, autodiff, _unwrap_val(concrete_jac)
    )
end
