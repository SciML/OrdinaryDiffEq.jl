"""
    NewmarkBeta

Classical Newmark-β method to solve second order ODEs, possibly in mass matrix form.
Local truncation errors are estimated with the estimate of Zienkiewicz and Xie.

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
        nlsolve = NewtonRaphson(), thread = Val{false}()
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
