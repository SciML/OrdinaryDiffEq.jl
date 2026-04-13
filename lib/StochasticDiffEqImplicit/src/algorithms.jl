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
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    return ImplicitEM{
        chunk_size, autodiff,
        typeof(linsolve), typeof(nlsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag),
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
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    return ImplicitEulerHeun{
        chunk_size, autodiff,
        typeof(linsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag),
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
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive, interpretation = SciMLBase.AlgorithmInterpretation.Ito
    )
    return ImplicitRKMil{
        chunk_size, autodiff,
        typeof(linsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag),
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
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    return ISSEM{
        chunk_size, autodiff,
        typeof(linsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag),
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
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
        controller = :Predictive
    )
    return ISSEulerHeun{
        chunk_size, autodiff,
        typeof(linsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag),
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
        chunk_size = 0, autodiff = true, diff_type = Val{:central},
        standardtag = Val{true}(), concrete_jac = nothing,
        precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :min_correct,
        new_jac_conv_bound = 1.0e-3, controller = :Predictive,
        ode_error_est = true
    )
    return SKenCarp{
        chunk_size, autodiff, typeof(linsolve), typeof(precs), diff_type,
        SciMLBase._unwrap_val(standardtag), SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve), typeof(new_jac_conv_bound), controller,
    }(
        linsolve, nlsolve, precs, smooth_est, extrapolant, new_jac_conv_bound,
        ode_error_est
    )
end
