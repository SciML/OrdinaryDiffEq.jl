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
struct NewmarkBeta{PT, F, CS, AD, FDT, ST, CJ, Thread} <:
    OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm{CS, AD, FDT, ST, CJ}
    β::PT
    γ::PT
    nlsolve::F
    autodiff::AD
    thread::Thread
end

function NewmarkBeta(β, γ; kwargs...)
    return NewmarkBeta(; β, γ, kwargs...)
end

# Needed for remake
function NewmarkBeta(;
        β = 0.25, γ = 0.5, chunk_size = Val{0}(),
        autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        nlsolve = NewtonRaphson(), thread = Val{false}()
    )
    AD_choice, chunk_size, diff_type = OrdinaryDiffEqCore._process_AD_choice(
        autodiff, chunk_size, diff_type
    )

    @assert concrete_jac === nothing "Using a aser-defined Jacobian in Newmark-β is not yet possible."
    @assert 0.0 ≤ β ≤ 0.5 "Beta outside admissible range [0, 0.5]"
    @assert 0.0 ≤ γ ≤ 1.0 "Gamma outside admissible range [0, 1.0]"

    return NewmarkBeta{
        typeof(β), typeof(nlsolve),
        _unwrap_val(chunk_size), typeof(AD_choice), autodiff, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(thread),
    }(
        β, γ,
        nlsolve,
        AD_choice,
        thread,
    )
end
