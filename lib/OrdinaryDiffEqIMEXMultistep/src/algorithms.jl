# IMEX Multistep methods

struct CNAB2{CS, AD, F, F2, P, FDT, ST, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
 linsolve::F
 nlsolve::F2
 precs::P
 extrapolant::Symbol
end

function CNAB2(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
     concrete_jac = nothing, diff_type = Val{:forward},
     linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
     extrapolant = :linear)
 CNAB2{
     _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
     typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(
     linsolve,
     nlsolve,
     precs,
     extrapolant)
end

struct CNLF2{CS, AD, F, F2, P, FDT, ST, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
 linsolve::F
 nlsolve::F2
 precs::P
 extrapolant::Symbol
end
function CNLF2(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
     concrete_jac = nothing, diff_type = Val{:forward},
     linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
     extrapolant = :linear)
 CNLF2{
     _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
     typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(
     linsolve,
     nlsolve,
     precs,
     extrapolant)
end