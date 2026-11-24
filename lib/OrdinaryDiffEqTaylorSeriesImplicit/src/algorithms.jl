@doc generic_solver_docstring(
    "Implicit μ-Taylor and Taylor-Padé methods with simplified Newton solves.
    The adaptive error estimate uses an embedded Padé formula, the highest Taylor
    term, or divided differences, depending on the method. A Taylor-Padé method
    of degrees (P, Q) has order P + Q and the [P/Q] Padé stability function.
    P = Q gives A-stability and P < Q gives L-stability.", # description
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
    - `extrapolant`: Newton predictor. `:extrapolate` (default) extrapolates the
      state history, `:constant` uses the last state, `:linear` adds `dt * f(uprev)`,
      and `:polynomial` uses a second-order Taylor predictor (P, Q > 0 and order > 2).
    - `controller`: `:Predictive` (default, Gustafsson) or `:PI`. Slow and
      recovered nonlinear solves limit step growth.
    - `κ`: Newton tolerance relative to the integration tolerances (default
      `0.01`). Adaptive steps also check the filtered residual.
    - `maxiters`: maximum Newton corrections per attempt (default `10`).
    - `newton_matrix`: `:assembled` (default), `:factored` or `:auto`. The
      factored increment form uses the Padé denominator factors `I - hJ/r_i`,
      with `r_i` its roots. `:auto` can switch based on conditioning and convergence.
      These alternatives apply to adaptive Taylor-Padé methods with Q ≥ 3
      and real dense matrices.
    - `stages`: `:reduced` solves for the state increment. `:derivative` also
      solves for the scaled Taylor coefficients, with one Newton factor per
      real Padé root or conjugate pair. `:auto` (default) starts with the reduced
      form and switches based on Newton convergence. Derivative stages require
      real Taylor-Padé equations with Q ≥ 2 and `newton_matrix = :assembled`.
    - `real_function`: project μ-Taylor steps onto the real axis (default `true`).
      Set to `false` for complex-valued solutions.
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
    newton_matrix::Symbol
    stages::Symbol
end

# order_q = Val(0): μ-Taylor method; otherwise Taylor-Padé with degrees (P, Q).
function ImplicitTaylor(;
        order = Val(1), order_q = Val(0), μ = 1.0, real_function = true,
        linsolve = nothing,
        extrapolant = :extrapolate,
        controller = :Predictive, step_limiter! = trivial_limiter!, κ = nothing, maxiters = 10,
        newton_matrix = :assembled, stages = :auto
    )
    P, Q = _unwrap_val(order), _unwrap_val(order_q)
    P >= 0 && Q >= 0 && P + Q > 0 || throw(ArgumentError("Taylor degrees must be nonnegative with P + Q > 0."))
    maxiters > 0 || throw(ArgumentError("maxiters must be positive."))
    κ === nothing || (isfinite(κ) && κ > 0) || throw(ArgumentError("κ must be finite and positive."))
    extrapolant in (:extrapolate, :constant, :linear, :polynomial) ||
        throw(ArgumentError("extrapolant must be :extrapolate, :constant, :linear or :polynomial."))
    controller in (:PI, :Predictive) || throw(ArgumentError("controller must be :PI or :Predictive."))
    newton_matrix in (:auto, :assembled, :factored) ||
        throw(ArgumentError("newton_matrix must be :auto, :assembled or :factored."))
    stages in (:auto, :derivative, :reduced) ||
        throw(ArgumentError("stages must be :auto, :derivative or :reduced."))
    return ImplicitTaylor{
        P, Q, typeof(μ), typeof(linsolve),
        typeof(κ), typeof(step_limiter!),
    }(
        order, order_q, μ, real_function, linsolve,
        extrapolant, controller, step_limiter!, κ, maxiters, newton_matrix, stages
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
