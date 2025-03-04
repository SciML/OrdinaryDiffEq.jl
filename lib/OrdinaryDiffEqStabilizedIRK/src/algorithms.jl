@doc generic_solver_docstring("Implicit Runge-Kutta-Chebyshev method.",
    "IRKC",
    "Stabilized Implicit Runge Kutta method.",
    "REF TBD",
    "- `eigen_est`: function of the form
    `(integrator) -> integrator.eigen_est = upper_bound`,
    where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
    If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.",
    "eigen_est = nothing,")
struct IRKC{CS, AD, F, F2, P, FDT, ST, CJ, K, T, E} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    controller::Symbol
    eigen_est::E
    autodiff::AD
end

function IRKC(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, controller = :Standard, eigen_est = nothing)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    IRKC{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol), typeof(eigen_est)}(linsolve, nlsolve, precs, κ, tol,
        extrapolant, controller, eigen_est, AD_choice)
end
