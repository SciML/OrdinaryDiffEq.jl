@doc generic_solver_docstring(
    "Implicit Taylor and Taylor-Padé solvers. Each step solves one nonlinear
    system with a simplified Newton iteration whose matrix is assembled from the
    Jacobians of the Taylor coefficients; its factorization is reused while the
    matrix changes little. The local error is estimated from an embedded
    Taylor-Padé formula (order ≥ 3) or from the next Taylor term (μ-Taylor),
    filtered through `(I - hJ/(P+Q))^Q` with `J` the Jacobian of `f`; the first-
    and second-order methods use divided differences of the last states. The
    Taylor-Padé method with degrees (P, Q) has the [P/Q] Padé stability function:
    the Gauss methods (P = Q) are A-stable with |R(∞)| = 1 and do not damp the
    error of stiff components; the Radau and Lobatto methods (P < Q) are L-stable.", # description
    "ImplicitTaylor", # name
    "Implicit Taylor Series Method", # solver_class
    "@techreport{scott_solving_2000,
	title = {Solving ODE Initial Value Problems With Implicit Taylor Series Methods},
	url = {https://ntrs.nasa.gov/citations/20000034027},
	number = {NASA/TM-2000-209400},
	author = {Scott, James R.},
	month = mar,
	year = {2000},
    }", # references
    """
    - `extrapolant`: initial guess of the Newton iteration. `:extrapolate`
      (default) extrapolates the last three states, which stays on the slow
      manifold of a stiff problem; `:constant` starts at the previous state;
      `:linear` uses `dt * f(uprev)` and `:polynomial` a second-order Taylor
      predictor (P, Q > 0, order > 2 only). The last two are built from
      derivatives and can start far off the manifold on stiff components.
    - `controller`: `:Predictive` (default, Gustafsson) or `:PI`. dt does not grow
      after a step whose Newton iteration failed, needed recovery or converged
      slowly, and after a failure it grows back only gradually.
    - `κ`: Newton stopping tolerance relative to the integration tolerances
      (default `0.01`): the iteration stops when `η * ‖Δz‖ < κ`. Adaptive steps
      also validate the residual through the error filter before accepting.
    - `maxiters`: maximum Newton corrections per attempt (default `10`).
    - `real_function`: project μ-Taylor steps onto the real axis (default `true`);
      set to `false` for complex-valued solutions.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """, # keyword_description
    """
    extrapolant = :extrapolate,
    controller = :Predictive,
    step_limiter! = trivial_limiter!,
    """ # keyword_default
)
struct ImplicitTaylor{P, Q, T, F, Tol, StepLimiter} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    order::Val{P}
    order_q::Val{Q}
    μ::T
    real_function::Bool
    linsolve::F
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    κ::Tol
    maxiters::Int
end

# order_q = Val(0): μ-Taylor method; otherwise Taylor-Padé with degrees (P, Q).
function ImplicitTaylor(;
        order = Val(1), order_q = Val(0), μ = 1.0, real_function = true,
        linsolve = nothing,
        extrapolant = :extrapolate,
        controller = :Predictive, step_limiter! = trivial_limiter!, κ = nothing, maxiters = 10
    )
    P, Q = _unwrap_val(order), _unwrap_val(order_q)
    P >= 0 && Q >= 0 && P + Q > 0 || throw(ArgumentError("Taylor degrees must be nonnegative with P + Q > 0."))
    maxiters > 0 || throw(ArgumentError("maxiters must be positive."))
    κ === nothing || (isfinite(κ) && κ > 0) || throw(ArgumentError("κ must be finite and positive."))
    extrapolant in (:extrapolate, :constant, :linear, :polynomial) ||
        throw(ArgumentError("extrapolant must be :extrapolate, :constant, :linear or :polynomial."))
    controller in (:PI, :Predictive) || throw(ArgumentError("controller must be :PI or :Predictive."))
    return ImplicitTaylor{
        P, Q, typeof(μ), typeof(linsolve),
        typeof(κ), typeof(step_limiter!),
    }(
        order, order_q, μ, real_function, linsolve,
        extrapolant, controller, step_limiter!, κ, maxiters
    )
end

# Polynomial degree 1, L-stable, first-order (same as implicit Euler)
ImplicitTaylor1L1(; kwargs...) = ImplicitTaylor(; order = Val(1), μ = 1.0, kwargs...)
# Polynomial degree 1, A-stable, second-order
ImplicitTaylor1A2(; kwargs...) = ImplicitTaylor(; order = Val(1), μ = 0.5, kwargs...)
# Polynomial degree 2, A-stable, second-order
ImplicitTaylor2A2(; kwargs...) = ImplicitTaylor(; order = Val(2), μ = 1.0, kwargs...)
# Polynomial degree 2, A-stable, second-order
ImplicitTaylor2A2Alt(; kwargs...) = ImplicitTaylor(; order = Val(2), μ = 0.5, kwargs...)
# Polynomial degree 2, A-stable, fourth-order
ImplicitTaylor2A4(; kwargs...) = ImplicitTaylor(; order = Val(2), μ = complex(0.5, √3 / 6), kwargs...)
# Polynomial degree 3, A-stable, fourth-order
ImplicitTaylor3A4(; kwargs...) = ImplicitTaylor(; order = Val(3), μ = 0.5, kwargs...)
# Polynomial degree 4, A-stable, fourth-order
ImplicitTaylor4A4(; kwargs...) = ImplicitTaylor(; order = Val(4), μ = 0.5, kwargs...)

export ImplicitTaylor, ImplicitTaylor1L1, ImplicitTaylor1A2, ImplicitTaylor2A2, ImplicitTaylor2A2Alt, ImplicitTaylor2A4, ImplicitTaylor3A4, ImplicitTaylor4A4

# Taylor-Padé methods: [P/Q] Padé stability function, order P + Q. Gauss
# (P = Q): A-stable, |R(∞)| = 1. Radau and Lobatto (P < Q): L-stable.
TaylorGaussA2(; kwargs...) = ImplicitTaylor(; order = Val(1), order_q = Val(1), kwargs...)
TaylorGaussA4(; kwargs...) = ImplicitTaylor(; order = Val(2), order_q = Val(2), kwargs...)
TaylorGaussA6(; kwargs...) = ImplicitTaylor(; order = Val(3), order_q = Val(3), kwargs...)

TaylorRadauLA1(; kwargs...) = ImplicitTaylor(; order = Val(0), order_q = Val(1), kwargs...)
TaylorRadauLA3(; kwargs...) = ImplicitTaylor(; order = Val(1), order_q = Val(2), kwargs...)
TaylorRadauLA5(; kwargs...) = ImplicitTaylor(; order = Val(2), order_q = Val(3), kwargs...)

TaylorLobattoLA2(; kwargs...) = ImplicitTaylor(; order = Val(0), order_q = Val(2), kwargs...)
TaylorLobattoLA4(; kwargs...) = ImplicitTaylor(; order = Val(1), order_q = Val(3), kwargs...)
TaylorLobattoLA6(; kwargs...) = ImplicitTaylor(; order = Val(2), order_q = Val(4), kwargs...)

export TaylorGaussA2, TaylorGaussA4, TaylorGaussA6, TaylorRadauLA1, TaylorRadauLA3, TaylorRadauLA5, TaylorLobattoLA2, TaylorLobattoLA4, TaylorLobattoLA6
