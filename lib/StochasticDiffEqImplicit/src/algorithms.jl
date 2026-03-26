struct ImplicitEM{CS, AD, F, F2, P, FDT, ST, CJ, T2, Controller} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::F2
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
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
    return ImplicitEM{
        OrdinaryDiffEqCore._ad_chunksize_int(autodiff), typeof(autodiff),
        typeof(linsolve), typeof(nlsolve), typeof(precs),
        OrdinaryDiffEqCore._ad_fdtype(autodiff),
        true,
        SciMLBase._unwrap_val(concrete_jac),
        typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant, new_jac_conv_bound, symplectic
    )
end

STrapezoid(; kwargs...) = ImplicitEM(; theta = 1 / 2, kwargs...)
SImplicitMidpoint(; kwargs...) = ImplicitEM(; theta = 1 / 2, symplectic = true, kwargs...)

struct ImplicitEulerHeun{CS, AD, F, P, FDT, ST, CJ, N, T2, Controller} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
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
    return ImplicitEulerHeun{
        OrdinaryDiffEqCore._ad_chunksize_int(autodiff), typeof(autodiff),
        typeof(linsolve), typeof(precs),
        OrdinaryDiffEqCore._ad_fdtype(autodiff),
        true,
        SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve),
        typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic
    )
end

struct ImplicitRKMil{CS, AD, F, P, FDT, ST, CJ, N, T2, Controller, interpretation} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
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
    return ImplicitRKMil{
        OrdinaryDiffEqCore._ad_chunksize_int(autodiff), typeof(autodiff),
        typeof(linsolve), typeof(precs),
        OrdinaryDiffEqCore._ad_fdtype(autodiff),
        true,
        SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve), typeof(new_jac_conv_bound),
        controller, interpretation,
    }(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic
    )
end

struct ISSEM{CS, AD, F, P, FDT, ST, CJ, N, T2, Controller} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
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
    return ISSEM{
        OrdinaryDiffEqCore._ad_chunksize_int(autodiff), typeof(autodiff),
        typeof(linsolve), typeof(precs),
        OrdinaryDiffEqCore._ad_fdtype(autodiff),
        true,
        SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve),
        typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic
    )
end

struct ISSEulerHeun{CS, AD, F, P, FDT, ST, CJ, N, T2, Controller} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::N
    precs::P
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T2
    symplectic::Bool
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
    return ISSEulerHeun{
        OrdinaryDiffEqCore._ad_chunksize_int(autodiff), typeof(autodiff),
        typeof(linsolve), typeof(precs),
        OrdinaryDiffEqCore._ad_fdtype(autodiff),
        true,
        SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve), typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic
    )
end

struct SKenCarp{CS, AD, F, P, FDT, ST, CJ, N, T2, Controller} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ, Controller}
    linsolve::F
    nlsolve::N
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    new_jac_conv_bound::T2
    ode_error_est::Bool
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
    return SKenCarp{
        OrdinaryDiffEqCore._ad_chunksize_int(autodiff), typeof(autodiff),
        typeof(linsolve), typeof(precs),
        OrdinaryDiffEqCore._ad_fdtype(autodiff),
        true, SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve), typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs, smooth_est, extrapolant, new_jac_conv_bound,
        ode_error_est
    )
end
