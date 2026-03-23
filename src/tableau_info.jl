"""
`Base.length(tab::ODERKTableau)`

Defines the length of a Runge-Kutta method to be the number of stages.
"""
Base.length(tab::ODERKTableau) = tab.stages

"""
`stability_region(z, tab::ODERKTableau; embedded = false)`

Calculates the stability function from the tableau at `z`. Stable if <1.
If `embedded = true`, the stability function is calculated for the embedded
method. Otherwise, the stability function is calculated for the main method
(default).

```math
r(z) = 1 + z bᵀ(I - zA)⁻¹ e
```
where e denotes a vector of ones.
"""
function stability_region(z::Number, tab::ODERKTableau; embedded::Bool = false)
    A, c = tab.A, tab.c
    b = embedded ? tab.αEEst : tab.α
    e = ones(eltype(A), length(b))
    stages = (I - z * A) \ e
    return 1 + z * (transpose(b) * stages)
end

"""
    stability_region(z, alg::AbstractODEAlgorithm)

Calculates the stability function from the algorithm `alg` at `z`.
The stability region of a possible embedded method cannot be calculated
using this method.

If you use an implicit method, you may run into convergence issues when
the value of `z` is outside of the stability region, e.g.,

```julia-repl
julia> typemin(Float64)
-Inf

julia> stability_region(typemin(Float64), ImplicitEuler())
┌ Warning: Newton steps could not converge and algorithm is not adaptive. Use a lower dt.

julia> nextfloat(typemin(Float64))
-1.7976931348623157e308

julia> stability_region(nextfloat(typemin(Float64)), ImplicitEuler())
0.0
```
"""
function stability_region(z::Number, alg::AbstractODEAlgorithm)
    u0 = one(z)
    tspan = (zero(real(z)), one(real(z)))
    ode = ODEProblem{false}((u, p, t) -> z * u, u0, tspan)
    integrator = init(ode, alg; adaptive = false, dt = 1)
    step!(integrator)
    return integrator.u[end]
end

"""
`stability_region(tab_or_alg::Union{ODERKTableau, AbstractODEAlgorithm};
                  initial_guess=-3.0)`

Calculates the length of the stability region in the real axis.
See also [`imaginary_stability_interval`](@ref).
"""
function stability_region(
        tab_or_alg::Union{ODERKTableau, AbstractODEAlgorithm};
        initial_guess = -3.0, kw...
    )
    f = (x, p) -> abs(stability_region(x, tab_or_alg)) - 1
    T = typeof(initial_guess)
    # Bracket between a point near 0 (inside stability region) and a point
    # far from 0 (outside stability region). Solvers are 0-stable (|r(0)|=1)
    # and not ∞-stable, so a sign change must exist between these points.
    near_zero = sign(initial_guess) * one(T) / T(10)
    far = initial_guess * T(10)
    left, right = minmax(near_zero, far)
    if f(left, nothing) * f(right, nothing) >= 0
        # No sign change on same side of 0 (e.g. A-stable methods).
        # Try bracket spanning 0: boundary is at z=0 for A-stable methods.
        left, right = minmax(far, -near_zero)
    end
    prob = IntervalNonlinearProblem{false}(f, (left, right))
    sol = solve(prob, Bisection(); kw...)
    return sol.u
end

"""
    imaginary_stability_interval(tab::ODERKTableau;
                                 initial_guess = length(tab) - 1)

Calculates the length of the imaginary stability interval, i.e.,
the size of the stability region on the imaginary axis.
See also [`stability_region`](@ref).
"""
function imaginary_stability_interval(
        tab::ODERKTableau;
        initial_guess = length(tab) - one(eltype(tab.A)),
        kw...
    )
    f = (x, p) -> abs(stability_region(im * x, tab)) - 1
    T = typeof(initial_guess)
    near_zero = one(T) / T(10)
    far = initial_guess * T(10)
    prob = IntervalNonlinearProblem{false}(f, (near_zero, far))
    sol = solve(prob, Bisection(); kw...)
    return sol.u
end

"""
    imaginary_stability_interval(alg::AbstractODEAlgorithm;
                                 initial_guess = 20.0)

Calculates the length of the imaginary stability interval, i.e.,
the size of the stability region on the imaginary axis.
See also [`stability_region`](@ref).
"""
function imaginary_stability_interval(
        alg::AbstractODEAlgorithm;
        initial_guess = 20.0,
        kw...
    )
    f = (x, p) -> abs(stability_region(im * x, alg)) - 1
    T = typeof(initial_guess)
    near_zero = one(T) / T(10)
    far = initial_guess * T(10)
    prob = IntervalNonlinearProblem{false}(f, (near_zero, far))
    sol = solve(prob, Bisection(); kw...)
    return sol.u
end

function RootedTrees.residual_order_condition(
        tab::ODERKTableau, order::Integer,
        reducer = nothing, mapper = x -> x^2;
        embedded::Bool = false
    )
    A, c = tab.A, tab.c
    b = embedded ? tab.αEEst : tab.α
    if reducer === nothing
        resid = map(RootedTreeIterator(order)) do t
            residual_order_condition(t, A, b, c)
        end
    else
        resid = mapreduce(reducer, RootedTreeIterator(order)) do t
            mapper(residual_order_condition(t, RungeKuttaMethod(A, b, c)))
        end
    end
    return resid
end

isfsal(tab::ExplicitRKTableau) = tab.fsal
isfsal(::ImplicitRKTableau) = nothing
function check_tableau(tab::ODERKTableau; tol::Real = 10eps(1.0))
    order = all(i -> residual_order_condition(tab, i, +, abs) < tol, 1:(tab.order))
    if !order
        error("Tableau's order is not correct.")
    end
    if tab.adaptiveorder != 0
        embedded_order = all(
            i -> residual_order_condition(
                tab, i, +, abs;
                embedded = true
            ) < tol,
            tab.adaptiveorder
        )
        if !embedded_order
            error("Tableau's embedded order is not correct.")
        end
    end
    return true
end
