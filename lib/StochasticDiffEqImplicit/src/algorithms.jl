using ADTypes: AutoForwardDiff
using OrdinaryDiffEqCore: _fixup_ad, _unwrap_val

struct ImplicitEM{AD, F, F2, T2, T3, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T3
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
end
function ImplicitEM(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitEM(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant, new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac),
    )
end

STrapezoid(; kwargs...) = ImplicitEM(; theta = 1 / 2, kwargs...)
SImplicitMidpoint(; kwargs...) = ImplicitEM(; theta = 1 / 2, symplectic = true, kwargs...)

struct ImplicitEulerHeun{AD, F, N, T2, T3, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T3
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
end
function ImplicitEulerHeun(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitEulerHeun(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac),
    )
end

struct ImplicitRKMil{AD, F, N, T2, T3, interpretation, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T3
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
end
function ImplicitRKMil(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3, interpretation = SciMLBase.AlgorithmInterpretation.Ito
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitRKMil{
        typeof(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(symplectic ? 1 / 2 : theta), typeof(new_jac_conv_bound),
        interpretation, typeof(_unwrap_val(concrete_jac)),
    }(
        linsolve, nlsolve, symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac),
    )
end

struct ISSEM{AD, F, N, T2, T3, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T3
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
end
function ISSEM(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
    )
    autodiff = _fixup_ad(autodiff)
    return ISSEM(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac),
    )
end

struct ISSEulerHeun{AD, F, N, T2, T3, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T3
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
end
function ISSEulerHeun(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
    )
    autodiff = _fixup_ad(autodiff)
    return ISSEulerHeun(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac),
    )
end

struct SKenCarp{AD, F, N, T2, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    smooth_est::Bool
    extrapolant::Symbol
    new_jac_conv_bound::T2
    ode_error_est::Bool
    autodiff::AD
    concrete_jac::CJ
end

function SKenCarp(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :min_correct,
        new_jac_conv_bound = 1.0e-3,
        ode_error_est = true
    )
    autodiff = _fixup_ad(autodiff)
    return SKenCarp(
        linsolve, nlsolve, smooth_est, extrapolant, new_jac_conv_bound,
        ode_error_est,
        autodiff, _unwrap_val(concrete_jac),
    )
end
