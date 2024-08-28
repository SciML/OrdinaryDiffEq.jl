
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
struct NewmarkBeta{PT, F, F2, P, CS, AD, FDT, ST, CJ} <:
    OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm{CS, AD, FDT, ST, CJ}
    β::PT
    γ::PT
    linsolve::F
    nlsolve::F2
    precs::P
end

function NewmarkBeta(β, γ; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
    concrete_jac = nothing, diff_type = Val{:forward},
    linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
    extrapolant = :linear)
    NewmarkBeta{
        typeof(β), typeof(linsolve), typeof(nlsolve), typeof(precs),
        _unwrap_val(chunk_size), _unwrap_val(autodiff), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(
        β, γ,
        linsolve,
        nlsolve,
        precs)
end

# Needed for remake
function NewmarkBeta(; β=0.25, γ=0.5, chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
    concrete_jac = nothing, diff_type = Val{:forward},
    linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
    extrapolant = :linear)
    NewmarkBeta{
        typeof(β), typeof(linsolve), typeof(nlsolve), typeof(precs),
        _unwrap_val(chunk_size), _unwrap_val(autodiff), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(
        β, γ,
        linsolve,
        nlsolve,
        precs)
end
