hairer1999stiff = """@article{hairer1999stiff,
title={Stiff differential equations solved by Radau methods},
author={Hairer, Ernst and Wanner, Gerhard},
journal={Journal of Computational and Applied Mathematics},
volume={111},
number={1-2},
pages={93--111},
year={1999},
publisher={Elsevier}}"""

AdaptiveRadauPaper = """@article{AdaptiveRadauPaper,
author={Ekanathan, Shreyas and Smith, Oscar and Rackauckas, Christopher},
booktitle={2025 IEEE High Performance Extreme Computing Conference (HPEC)}, 
title={A Fully Adaptive Radau Method for the Efficient Solution of Stiff Ordinary Differential Equations at Low Tolerances}, 
year={2025},
pages={1-9},
doi={10.1109/HPEC67600.2025.11196706}}"""

extra_keyword_description = """
- `extrapolant`: TBD
- `smooth_est`: TBD
- `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`"""
extra_keyword_default = """
extrapolant = :dense,
smooth_est = true,
step_limiter! = trivial_limiter!"""

@doc differentiation_rk_docstring(
    "An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency.
Similar to Hairer's SEULEX.",
    "RadauIIA3",
    "Fully-Implicit Runge-Kutta Method.";
    references = hairer1999stiff,
    extra_keyword_description = extra_keyword_description,
    extra_keyword_default = extra_keyword_default
)
struct RadauIIA3{CS, AD, F, P, FDT, ST, CJ, Tol, C1, C2, StepLimiter} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    extrapolant::Symbol
    κ::Tol
    maxiters::Int
    fast_convergence_cutoff::C1
    new_W_γdt_cutoff::C2
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
end

function RadauIIA3(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return RadauIIA3{
        _ad_chunksize_int(autodiff), typeof(autodiff), typeof(linsolve),
        typeof(precs), _ad_fdtype(autodiff), true, _unwrap_val(concrete_jac),
        typeof(κ), typeof(fast_convergence_cutoff),
        typeof(new_W_γdt_cutoff), typeof(step_limiter!),
    }(
        linsolve,
        precs,
        extrapolant,
        κ,
        maxiters,
        fast_convergence_cutoff,
        new_W_γdt_cutoff,
        controller,
        step_limiter!,
        autodiff
    )
end

@doc differentiation_rk_docstring(
    "An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency. 5th order method with excellent numerical stability. Good for highly stiff systems, problems requiring high-order implicit integration, systems with complex eigenvalue structures. Best for low tolerance stiff problems (<1e-9).",
    "RadauIIA5",
    "Fully-Implicit Runge-Kutta Method.";
    references = hairer1999stiff,
    extra_keyword_description = extra_keyword_description,
    extra_keyword_default = extra_keyword_default
)
struct RadauIIA5{CS, AD, F, P, FDT, ST, CJ, Tol, C1, C2, StepLimiter} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    κ::Tol
    maxiters::Int
    fast_convergence_cutoff::C1
    new_W_γdt_cutoff::C2
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
end

function RadauIIA5(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10, smooth_est = true,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return RadauIIA5{
        _ad_chunksize_int(autodiff), typeof(autodiff), typeof(linsolve),
        typeof(precs), _ad_fdtype(autodiff), true, _unwrap_val(concrete_jac),
        typeof(κ), typeof(fast_convergence_cutoff),
        typeof(new_W_γdt_cutoff), typeof(step_limiter!),
    }(
        linsolve,
        precs,
        smooth_est,
        extrapolant,
        κ,
        maxiters,
        fast_convergence_cutoff,
        new_W_γdt_cutoff,
        controller,
        step_limiter!,
        autodiff
    )
end

@doc differentiation_rk_docstring(
    "An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency.
Similar to Hairer's SEULEX.",
    "RadauIIA9",
    "Fully-Implicit Runge-Kutta Method.";
    references = hairer1999stiff,
    extra_keyword_description = extra_keyword_description,
    extra_keyword_default = extra_keyword_default
)
struct RadauIIA9{CS, AD, F, P, FDT, ST, CJ, Tol, C1, C2, StepLimiter} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    κ::Tol
    maxiters::Int
    fast_convergence_cutoff::C1
    new_W_γdt_cutoff::C2
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
end

function RadauIIA9(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10, smooth_est = true,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return RadauIIA9{
        _ad_chunksize_int(autodiff), typeof(autodiff), typeof(linsolve),
        typeof(precs), _ad_fdtype(autodiff), true, _unwrap_val(concrete_jac),
        typeof(κ), typeof(fast_convergence_cutoff),
        typeof(new_W_γdt_cutoff), typeof(step_limiter!),
    }(
        linsolve,
        precs,
        smooth_est,
        extrapolant,
        κ,
        maxiters,
        fast_convergence_cutoff,
        new_W_γdt_cutoff,
        controller,
        step_limiter!,
        autodiff
    )
end
@doc differentiation_rk_docstring(
    "An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency.
     Fully autonomous derivation of tableau for arbitrary order and order adaptivity.",
    "AdaptiveRadau",
    "Fully-Implicit Runge-Kutta Method.";
    references = AdaptiveRadauPaper,
    extra_keyword_description = extra_keyword_description,
    extra_keyword_default = extra_keyword_default
)
struct AdaptiveRadau{CS, AD, F, P, FDT, ST, CJ, Tol, C1, C2, StepLimiter, TO} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    κ::Tol
    maxiters::Int
    fast_convergence_cutoff::C1
    new_W_γdt_cutoff::C2
    controller::Symbol
    step_limiter!::StepLimiter
    min_order::Int
    max_order::Int
    threading::TO
    autodiff::AD
end

function AdaptiveRadau(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        min_order = 5, max_order = 13, threading = false,
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10, smooth_est = true,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return AdaptiveRadau{
        _ad_chunksize_int(autodiff), typeof(autodiff), typeof(linsolve),
        typeof(precs), _ad_fdtype(autodiff), true, _unwrap_val(concrete_jac),
        typeof(κ), typeof(fast_convergence_cutoff),
        typeof(new_W_γdt_cutoff), typeof(step_limiter!), typeof(threading),
    }(
        linsolve,
        precs,
        smooth_est,
        extrapolant,
        κ,
        maxiters,
        fast_convergence_cutoff,
        new_W_γdt_cutoff,
        controller,
        step_limiter!, min_order, max_order, threading,
        autodiff
    )
end
