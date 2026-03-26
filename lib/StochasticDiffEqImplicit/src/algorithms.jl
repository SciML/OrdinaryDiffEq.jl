struct ImplicitEM{AD, F, F2, P, T2} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
    autodiff::AD
    concrete_jac::Union{Nothing, Bool}
    controller::Symbol
end
function ImplicitEM(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitEM(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant, new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac), controller
    )
end

STrapezoid(; kwargs...) = ImplicitEM(; theta = 1 / 2, kwargs...)
SImplicitMidpoint(; kwargs...) = ImplicitEM(; theta = 1 / 2, symplectic = true, kwargs...)

struct ImplicitEulerHeun{AD, F, P, N, T2} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
    autodiff::AD
    concrete_jac::Union{Nothing, Bool}
    controller::Symbol
end
function ImplicitEulerHeun(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitEulerHeun(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac), controller
    )
end

struct ImplicitRKMil{AD, F, P, N, T2, interpretation} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
    autodiff::AD
    concrete_jac::Union{Nothing, Bool}
    controller::Symbol
end
function ImplicitRKMil(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive, interpretation = SciMLBase.AlgorithmInterpretation.Ito
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitRKMil{
        typeof(autodiff), typeof(linsolve), typeof(precs), typeof(nlsolve),
        typeof(symplectic ? 1 / 2 : theta), typeof(interpretation),
    }(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac), controller
    )
end

struct ISSEM{AD, F, P, N, T2} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
    autodiff::AD
    concrete_jac::Union{Nothing, Bool}
    controller::Symbol
end
function ISSEM(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    autodiff = _fixup_ad(autodiff)
    return ISSEM(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac), controller
    )
end

struct ISSEulerHeun{AD, F, P, N, T2} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
    autodiff::AD
    concrete_jac::Union{Nothing, Bool}
    controller::Symbol
end
function ISSEulerHeun(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    autodiff = _fixup_ad(autodiff)
    return ISSEulerHeun(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac), controller
    )
end

struct SKenCarp{AD, F, P, N, T2} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    new_jac_conv_bound::T2
    ode_error_est::Bool
    autodiff::AD
    concrete_jac::Union{Nothing, Bool}
    controller::Symbol
end

function SKenCarp(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :min_correct,
        new_jac_conv_bound = 1.0e-3, controller = :Predictive,
        ode_error_est = true
    )
    autodiff = _fixup_ad(autodiff)
    return SKenCarp(
        linsolve, nlsolve, precs, smooth_est, extrapolant, new_jac_conv_bound,
        ode_error_est,
        autodiff, _unwrap_val(concrete_jac), controller
    )
end
