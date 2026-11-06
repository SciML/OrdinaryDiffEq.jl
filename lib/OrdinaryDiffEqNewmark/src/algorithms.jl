"""
    NewmarkBeta

Classical Newmark-β method to solve second order ODEs, possibly in mass matrix form.
Local truncation errors are estimated with the estimate of Zienkiewicz and Xie.

## Adaptive time stepping

The Zienkiewicz-Xie estimate monitors both solution components. The position
channel uses the leading local error (β - 1/6)·dt²·(aₙ₊₁ - aₙ). The velocity
channel uses (γ - 1/2)·dt·(aₙ₊₁ - aₙ) with the coefficient magnitude clamped
at 1/100 from below: on and near the second-order family γ = 1/2 the true
coefficient vanishes, and the channel falls back to the surrogate
dt·(aₙ₊₁ - aₙ)/100, one order below the true velocity error there. This keeps
the velocity error tracking the tolerance at all frequencies — measured on
u'' = -ω²u for ω = 1..1000, tol = 1e-4..1e-10, the global error stays below
2.5·nsteps·tol everywhere, where without the clamp the velocity error reached
the full peak-to-peak range with a Success return at ω = 1000 — at the price
of dt ~ tol^(1/2) instead of tol^(1/3) where the surrogate binds: step counts
at ω = 1 grow by 1.05x (tol = 1e-4) up to 9.1x (tol = 1e-10). Combinations
with β within 1/24 of 1/6 and γ within 1/100 of 1/2 are rejected in adaptive
mode because neither leading error coefficient is measurable there; the
linear-acceleration method β = 1/6, γ = 1/2 remains available with
`adaptive = false`.

## References

Newmark, Nathan (1959), "A method of computation for structural dynamics",
Journal of the Engineering Mechanics Division, 85 (EM3) (3): 67–94, doi:
https://doi.org/10.1061/JMCEA3.0000098

Zienkiewicz, O. C., and Y. M. Xie. "A simple error estimator and adaptive
time stepping procedure for dynamic analysis." Earthquake engineering &
structural dynamics 20.9 (1991): 871-887, doi:
https://doi.org/10.1002/eqe.4290200907
"""
struct NewmarkBeta{PT, F, AD, Thread, CJ} <:
    OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm
    β::PT
    γ::PT
    nlsolve::F
    autodiff::AD
    thread::Thread
    concrete_jac::CJ
end

function NewmarkBeta(β, γ; kwargs...)
    return NewmarkBeta(; β, γ, kwargs...)
end

# Needed for remake
function NewmarkBeta(;
        β = 0.25, γ = 0.5,
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        nlsolve = NewtonRaphson(), thread = Serial()
    )
    autodiff = OrdinaryDiffEqCore._fixup_ad(autodiff)

    @assert concrete_jac === nothing "Using a user-defined Jacobian in Newmark-β is not yet possible."
    @assert 0.0 ≤ β ≤ 0.5 "Beta outside admissible range [0, 0.5]"
    @assert 0.0 ≤ γ ≤ 1.0 "Gamma outside admissible range [0, 1.0]"

    return NewmarkBeta(
        β, γ,
        nlsolve,
        autodiff,
        thread,
        _unwrap_val(concrete_jac)

    )
end

"""
    GeneralizedAlpha

Generalized-α method for second-order ODEs in mass-matrix form, due to Chung & Hulbert (1993).
Encompasses Newmark-β (αₘ = αf = 0) and HHT-α (αₘ = 0) as special cases.

The method evaluates the equations of motion at interpolated states:

    M · aₙ₊αₘ = f(uₙ₊αf, vₙ₊αf, tₙ₊αf)

where:

    aₙ₊αₘ = (1 - αₘ) · aₙ₊₁ + αₘ · aₙ
    uₙ₊αf = (1 - αf) · uₙ₊₁ + αf · uₙ
    vₙ₊αf = (1 - αf) · vₙ₊₁ + αf · vₙ

with the standard Newmark update formulas for uₙ₊₁ and vₙ₊₁.

## Constructors

    GeneralizedAlpha(; rho_inf)                  # spectral radius ρ∞ ∈ [0, 1]
    GeneralizedAlpha(αm, αf, β, γ)               # all four parameters directly
    GeneralizedAlpha(; alpha_hht)                 # HHT-α convenience (αₘ = 0)

### ρ∞ parameterization (recommended)

ρ∞ ∈ [0, 1] is the spectral radius at infinity (high-frequency damping).
ρ∞ = 1 → no algorithmic damping (identical to Newmark with γ = 1/2, β = 1/4).
ρ∞ = 0 → maximum algorithmic damping.

Parameters are set optimally (Chung & Hulbert 1993):

    αₘ = (2ρ∞ - 1) / (ρ∞ + 1)
    αf  = ρ∞        / (ρ∞ + 1)
    γ   = 1/2 - αₘ + αf
    β   = (1/2 + αf - αₘ)² / 4

### HHT-α convenience constructor

    GeneralizedAlpha(; alpha_hht = -0.1)   # α ∈ [-1/3, 0]

Sets αₘ = 0, αf = -α, γ = (1 - 2α)/2, β = (1 - α)²/4.

## Order and adaptive time stepping

This implementation re-derives the acceleration from `f` at the accepted state
on every step instead of carrying the algorithmic acceleration, so it is second
order iff γ(1 - αf) = (1 - αm)/2, which for αₘ ≠ αf differs from the textbook
condition γ = 1/2 - αₘ + αf; the ρ∞ and HHT parameterizations with damping
therefore converge at first order with a small leading constant (measured
fixed-step order 2.00 on the implementation locus).

The Zienkiewicz-Xie estimate is extended with the coefficients β - 1/(6κ)
(position) and γ - 1/(2κ) (velocity), κ = (1 - αf)/(1 - αm), the velocity
coefficient magnitude clamped at 1/100 from below so that the velocity channel
stays live on the family γ = 1/(2κ) where it otherwise vanishes. Large ρ∞
sits inside that clamp (ρ∞ = 0.9 has γ - 1/(2κ) ≈ 0.0026, ρ∞ = 1 is exactly
on the family); measured on u'' = -ω²u at ω = 100, ρ∞ = 0.9 holds the global
error below 0.7·nsteps·tol for tol = 1e-4..1e-6. Parameter sets with β
within 1/24 of 1/(6κ) and γ within 1/100 of 1/(2κ) are rejected in adaptive
mode because neither leading error coefficient is measurable there; they
remain available with `adaptive = false`.

## References

Chung, J., and Hulbert, G. M. (1993), "A time integration algorithm for structural
dynamics with improved numerical dissipation: The generalized-α method",
Journal of Applied Mechanics, 60(2): 371-375.
doi: https://doi.org/10.1115/1.2900803

Hilber, H. M., Hughes, T. J. R., and Taylor, R. L. (1977), "Improved numerical
dissipation for time integration algorithms in structural dynamics",
Earthquake Engineering & Structural Dynamics, 5(3): 283-292.
doi: https://doi.org/10.1002/eqe.4290050306
"""
struct GeneralizedAlpha{PT, F, AD, Thread, CJ} <:
    OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm
    αm::PT
    αf::PT
    β::PT
    γ::PT
    nlsolve::F
    autodiff::AD
    thread::Thread
    concrete_jac::CJ
end

# handles multiple ways the user can call the method: standard rho_inf, explicit 4 parameter construction, or the HHT-alpha convention
# The four parameters are also accepted as keywords because remake (reached
# through prepare_alg whenever a sublibrary with AD preparation is loaded
# alongside) reconstructs the algorithm from its field names.
function GeneralizedAlpha(
        αm_arg = nothing, αf_arg = nothing, β_arg = nothing, γ_arg = nothing;
        αm = αm_arg, αf = αf_arg, β = β_arg, γ = γ_arg,
        rho_inf = nothing,
        alpha_hht = nothing,
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        nlsolve = NewtonRaphson(),
        thread = Serial()
    )

    # rho inf case
    if !isnothing(rho_inf)
        @assert 0.0 ≤ rho_inf ≤ 1.0 "ρ∞ must be in [0, 1]; got $rho_inf"
        ρ = rho_inf
        αm = (2ρ - 1) / (ρ + 1)
        αf = ρ / (ρ + 1)
        γ = 1 / 2 - αm + αf
        β = (1 / 2 + αf - αm)^2 / 4

        # HHT-alpha case
    elseif !isnothing(alpha_hht)
        @assert -1 / 3 ≤ alpha_hht ≤ 0 "HHT α must be in [-1/3, 0]; got $alpha_hht"
        αm = zero(alpha_hht)
        αf = -alpha_hht
        γ = (1 - 2alpha_hht) / 2
        β = (1 - alpha_hht)^2 / 4

        # explicit 4 parameter case
    else
        @assert !isnothing(αm) && !isnothing(αf) && !isnothing(β) && !isnothing(γ) "Must provide either rho_inf, alpha_hht, or all four of αm, αf, β, γ"
    end

    autodiff = OrdinaryDiffEqCore._fixup_ad(autodiff)
    @assert concrete_jac === nothing "Using a user-defined Jacobian in GeneralizedAlpha is not yet possible."
    @assert αm ≤ αf "Unconditional stability requires αm ≤ αf; got αm=$αm, αf=$αf"
    @assert αf ≤ 0.5 "Unconditional stability requires αf ≤ 1/2; got αf=$αf"
    β_min = (1 / 2 + αf - αm)^2 / 4
    @assert β ≥ β_min "Unconditional stability requires β ≥ $(β_min); got β=$β"
    @assert 0.0 ≤ γ ≤ 1.0 "γ outside admissible range [0, 1]; got γ=$γ"

    return GeneralizedAlpha(αm, αf, β, γ, nlsolve, autodiff, thread, _unwrap_val(concrete_jac))
end
