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
struct IRKC{AD, F, F2, K, T, E, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    κ::K
    tol::T
    extrapolant::Symbol
    eigen_est::E
    autodiff::AD
    concrete_jac::CJ
end

function IRKC(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, eigen_est = nothing
    )
    autodiff = _fixup_ad(autodiff)

    return IRKC(
        linsolve, nlsolve, κ, tol,
        extrapolant, eigen_est, autodiff,
        _unwrap_val(concrete_jac)
    )
end
