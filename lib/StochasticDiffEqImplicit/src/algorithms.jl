using ADTypes: AutoForwardDiff
using OrdinaryDiffEqCore: _fixup_ad, _unwrap_val

struct ImplicitEM{AD, F, F2, P, T2, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
    controller::Symbol
end
function ImplicitEM(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitEM(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant, new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac), controller
    )
end

STrapezoid(; kwargs...) = ImplicitEM(; theta = 1 / 2, kwargs...)
SImplicitMidpoint(; kwargs...) = ImplicitEM(; theta = 1 / 2, symplectic = true, kwargs...)

struct ImplicitEulerHeun{AD, F, N, T2, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
    controller::Symbol
end
function ImplicitEulerHeun(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitEulerHeun(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac), controller
    )
end

struct ImplicitRKMil{AD, F, N, T2, interpretation, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
    controller::Symbol
end
function ImplicitRKMil(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive, interpretation = SciMLBase.AlgorithmInterpretation.Ito
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitRKMil{
        typeof(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(symplectic ? 1 / 2 : theta), typeof(interpretation),
        _unwrap_val(concrete_jac),
    }(
        linsolve, nlsolve, symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac), controller
    )
end

struct ISSEM{AD, F, N, T2, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
    controller::Symbol
end
function ISSEM(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    autodiff = _fixup_ad(autodiff)
    return ISSEM(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac), controller
    )
end

struct ISSEulerHeun{AD, F, N, T2, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
    controller::Symbol
end
function ISSEulerHeun(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    autodiff = _fixup_ad(autodiff)
    return ISSEulerHeun(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac), controller
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
    controller::Symbol
end

function SKenCarp(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :min_correct,
        new_jac_conv_bound = 1.0e-3, controller = :Predictive,
        ode_error_est = true
    )
    autodiff = _fixup_ad(autodiff)
    return SKenCarp(
        linsolve, nlsolve, smooth_est, extrapolant, new_jac_conv_bound,
        ode_error_est,
        autodiff, _unwrap_val(concrete_jac), controller
    )
end
