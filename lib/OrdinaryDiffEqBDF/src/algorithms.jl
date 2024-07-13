"""
Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton. Implicit-Explicit Methods for Time-
Dependent Partial Differential Equations. 1995 Society for Industrial and Applied Mathematics
Journal on Numerical Analysis, 32(3), pp 797-823, 1995. doi: https://doi.org/10.1137/0732037
"""
struct SBDF{CS, AD, F, F2, P, FDT, ST, CJ, K, T} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    order::Int
    ark::Bool
end

function SBDF(order; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, ark = false)
    SBDF{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol)}(linsolve,
        nlsolve,
        precs,
        κ,
        tol,
        extrapolant,
        order,
        ark)
end

# All keyword form needed for remake
function SBDF(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear,
        order, ark = false)
    SBDF{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol)}(linsolve,
        nlsolve,
        precs,
        κ,
        tol,
        extrapolant,
        order,
        ark)
end

"""
E. Alberdi Celayaa, J. J. Anza Aguirrezabalab, P. Chatzipantelidisc. Implementation of
an Adaptive BDF2 Formula and Comparison with The MATLAB Ode15s. Procedia Computer Science,
29, pp 1014-1026, 2014. doi: https://doi.org/10.1016/j.procs.2014.05.091

ABDF2: Multistep Method
An adaptive order 2 L-stable fixed leading coefficient multistep BDF method.
"""
struct ABDF2{CS, AD, F, F2, P, FDT, ST, CJ, K, T, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function ABDF2(; chunk_size = Val{0}(), autodiff = true, standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        κ = nothing, tol = nothing, linsolve = nothing, precs = DEFAULT_PRECS,
        nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :Standard, step_limiter! = trivial_limiter!)
    ABDF2{
        _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol), typeof(step_limiter!)}(linsolve, nlsolve, precs, κ, tol,
        smooth_est, extrapolant, controller, step_limiter!)
end

"""
QNDF1: Multistep Method
An adaptive order 1 quasi-constant timestep L-stable numerical differentiation function (NDF) method.
Optional parameter kappa defaults to Shampine's accuracy-optimal -0.1850.

See also `QNDF`.
"""
struct QNDF1{CS, AD, F, F2, P, FDT, ST, CJ, κType, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    kappa::κType
    controller::Symbol
    step_limiter!::StepLimiter
end

function QNDF1(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, kappa = -0.1850,
        controller = :Standard, step_limiter! = trivial_limiter!)
    QNDF1{
        _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(kappa), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        kappa,
        controller,
        step_limiter!)
end

"""
QNDF2: Multistep Method
An adaptive order 2 quasi-constant timestep L-stable numerical differentiation function (NDF) method.

See also `QNDF`.
"""
struct QNDF2{CS, AD, F, F2, P, FDT, ST, CJ, κType, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    kappa::κType
    controller::Symbol
    step_limiter!::StepLimiter
end

function QNDF2(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, kappa = -1 // 9,
        controller = :Standard, step_limiter! = trivial_limiter!)
    QNDF2{
        _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(kappa), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        kappa,
        controller,
        step_limiter!)
end

"""
QNDF: Multistep Method
An adaptive order quasi-constant timestep NDF method.
Utilizes Shampine's accuracy-optimal kappa values as defaults (has a keyword argument for a tuple of kappa coefficients).

@article{shampine1997matlab,
title={The matlab ode suite},
author={Shampine, Lawrence F and Reichelt, Mark W},
journal={SIAM journal on scientific computing},
volume={18},
number={1},
pages={1--22},
year={1997},
publisher={SIAM}
}
"""
struct QNDF{MO, CS, AD, F, F2, P, FDT, ST, CJ, K, T, κType, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    kappa::κType
    controller::Symbol
    step_limiter!::StepLimiter
end

function QNDF(; max_order::Val{MO} = Val{5}(), chunk_size = Val{0}(),
        autodiff = Val{true}(), standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, kappa = promote(-0.1850, -1 // 9, -0.0823, -0.0415, 0),
        controller = :Standard, step_limiter! = trivial_limiter!) where {MO}
    QNDF{MO, _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol), typeof(kappa), typeof(step_limiter!)}(
        max_order, linsolve, nlsolve, precs, κ, tol,
        extrapolant, kappa, controller, step_limiter!)
end

"""
MEBDF2: Multistep Method
The second order Modified Extended BDF method, which has improved stability properties over the standard BDF.
Fixed timestep only.
"""
struct MEBDF2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end
function MEBDF2(; chunk_size = Val{0}(), autodiff = true, standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant)
    MEBDF2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end