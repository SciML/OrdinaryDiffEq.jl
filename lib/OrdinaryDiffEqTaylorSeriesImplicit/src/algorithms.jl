@doc generic_solver_docstring(
    "Taylor implicit solvers. A-stable.
    Adaptive timestepping through a divided differences estimate via memory.", # description
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
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """, # keyword_description
    """
    extrapolant = :constant,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """ # keyword_default
)
struct ImplicitTaylor{P, Q, T, F, Prec, Tol, StepLimiter} <:
    OrdinaryDiffEqAdaptiveAlgorithm
    order::Val{P}
    order_q::Val{Q}
    μ::T
    real_function::Bool
    linsolve::F
    precs::Prec
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    κ::Tol
    maxiters::Int
end

# This simultaneously represents Taylor-Padé methods and μ-Taylor methods.
# When order_q = Val(0), it's a μ-Taylor method; otherwise, it's a Taylor-Padé method.
function ImplicitTaylor(;
        order = Val(1), order_q = Val(0), μ = 1.0, real_function = true,
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :constant,
        controller = :PI, step_limiter! = trivial_limiter!, κ = nothing, maxiters = 10
    )

    return ImplicitTaylor{
        _unwrap_val(order), _unwrap_val(order_q), typeof(μ), typeof(linsolve),
        typeof(precs), typeof(κ), typeof(step_limiter!),
    }(
        order, order_q, μ, real_function, linsolve,
        precs, extrapolant, controller, step_limiter!, κ, maxiters
    )
end

# Polynomial degree 1, L-stable, first-order (same as implicit Euler)
ImplicitTaylor1L1(; kwargs...) = ImplicitTaylor(order = Val(1), μ = 1.0, kwargs...)
# Polynomial degree 1, A-stable, second-order
ImplicitTaylor1A2(; kwargs...) = ImplicitTaylor(order = Val(1), μ = 0.5, kwargs...)
# Polynomial degree 2, A-stable, second-order
ImplicitTaylor2A2(; kwargs...) = ImplicitTaylor(order = Val(2), μ = 1.0, kwargs...)
# Polynomial degree 2, A-stable, second-order
ImplicitTaylor2A2Alt(; kwargs...) = ImplicitTaylor(order = Val(2), μ = 0.5, kwargs...)
# Polynomial degree 2, A-stable, fourth-order
ImplicitTaylor2A4(; kwargs...) = ImplicitTaylor(order = Val(2), μ = complex(0.5, √3 / 6), kwargs...)
# Polynomial degree 3, A-stable, fourth-order
ImplicitTaylor3A4(; kwargs...) = ImplicitTaylor(order = Val(3), μ = 0.5, kwargs...)
# Polynomial degree 4, A-stable, fourth-order
ImplicitTaylor4A4(; kwargs...) = ImplicitTaylor(order = Val(4), μ = 0.5, kwargs...)

export ImplicitTaylor, ImplicitTaylor1L1, ImplicitTaylor1A2, ImplicitTaylor2A2, ImplicitTaylor2A2Alt, ImplicitTaylor2A4, ImplicitTaylor3A4, ImplicitTaylor4A4

TaylorGaussA2(; kwargs...) = ImplicitTaylor(order = Val(1), order_q = Val(1), kwargs...)
TaylorGaussA4(; kwargs...) = ImplicitTaylor(order = Val(2), order_q = Val(2), kwargs...)
TaylorGaussA6(; kwargs...) = ImplicitTaylor(order = Val(3), order_q = Val(3), kwargs...)

TaylorRadauLA1(; kwargs...) = ImplicitTaylor(order = Val(0), order_q = Val(1), kwargs...)
TaylorRadauLA3(; kwargs...) = ImplicitTaylor(order = Val(1), order_q = Val(2), kwargs...)
TaylorRadauLA5(; kwargs...) = ImplicitTaylor(order = Val(2), order_q = Val(3), kwargs...)

TaylorLobattoLA2(; kwargs...) = ImplicitTaylor(order = Val(0), order_q = Val(2), kwargs...)
TaylorLobattoLA4(; kwargs...) = ImplicitTaylor(order = Val(1), order_q = Val(3), kwargs...)
TaylorLobattoLA6(; kwargs...) = ImplicitTaylor(order = Val(2), order_q = Val(4), kwargs...)

export TaylorGaussA2, TaylorGaussA4, TaylorGaussA6, TaylorRadauLA1, TaylorRadauLA3, TaylorRadauLA5, TaylorLobattoLA2, TaylorLobattoLA4, TaylorLobattoLA6
