@doc generic_solver_docstring(
    "Implicit Runge-Kutta-Chebyshev method.",
    "IRKC",
    "Stabilized Implicit Runge Kutta method.",
    "REF TBD",
    "- `eigen_est`: function of the form
    `(integrator) -> integrator.eigen_est = upper_bound`,
    where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
    If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.",
    "eigen_est = nothing,"
)
struct IRKC{AD, F, F2, P, K, T, E} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    controller::Symbol
    eigen_est::E
    autodiff::AD
    concrete_jac::Union{Nothing, Bool}
end

function IRKC(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, controller = :Standard, eigen_est = nothing
    )
    autodiff = _fixup_ad(autodiff)

    return IRKC(
        linsolve, nlsolve, precs, κ, tol,
        extrapolant, controller, eigen_est, autodiff,
        _unwrap_val(concrete_jac)

    )
end
