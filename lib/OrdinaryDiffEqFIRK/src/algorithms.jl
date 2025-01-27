
hairer1999stiff = """@article{hairer1999stiff,
    title={Stiff differential equations solved by Radau methods},
    author={Hairer, Ernst and Wanner, Gerhard},
    journal={Journal of Computational and Applied Mathematics},
    volume={111},
    number={1-2},
    pages={93--111},
    year={1999},
    publisher={Elsevier}}"""

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
    extra_keyword_default = extra_keyword_default)
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

function RadauIIA3(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10,
        step_limiter! = trivial_limiter!)

    AD_choice = _process_AD_choice(autodiff, chunk_size, diff_type)

    RadauIIA3{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(fast_convergence_cutoff),
        typeof(new_W_γdt_cutoff), typeof(step_limiter!)}(linsolve,
        precs,
        extrapolant,
        κ,
        maxiters,
        fast_convergence_cutoff,
        new_W_γdt_cutoff,
        controller,
        step_limiter!,
        AD_choice)
end

@doc differentiation_rk_docstring(
    "An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency.
Similar to Hairer's SEULEX.",
    "RadauIIA5",
    "Fully-Implicit Runge-Kutta Method.";
    references = hairer1999stiff,
    extra_keyword_description = extra_keyword_description,
    extra_keyword_default = extra_keyword_default)
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

function RadauIIA5(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10, smooth_est = true,
        step_limiter! = trivial_limiter!)

    AD_choice = _process_AD_choice(autodiff, chunk_size, diff_type)

    RadauIIA5{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(fast_convergence_cutoff),
        typeof(new_W_γdt_cutoff), typeof(step_limiter!)}(linsolve,
        precs,
        smooth_est,
        extrapolant,
        κ,
        maxiters,
        fast_convergence_cutoff,
        new_W_γdt_cutoff,
        controller,
        step_limiter!,
        AD_choice)
end

@doc differentiation_rk_docstring(
    "An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency.
Similar to Hairer's SEULEX.",
    "RadauIIA9",
    "Fully-Implicit Runge-Kutta Method.";
    references = hairer1999stiff,
    extra_keyword_description = extra_keyword_description,
    extra_keyword_default = extra_keyword_default)
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

function RadauIIA9(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10, smooth_est = true,
        step_limiter! = trivial_limiter!)

    AD_choice = _process_AD_choice(autodiff, chunk_size, diff_type)

    RadauIIA9{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(fast_convergence_cutoff),
        typeof(new_W_γdt_cutoff), typeof(step_limiter!)}(linsolve,
        precs,
        smooth_est,
        extrapolant,
        κ,
        maxiters,
        fast_convergence_cutoff,
        new_W_γdt_cutoff,
        controller,
        step_limiter!,
        AD_choice)
end

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

function AdaptiveRadau(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
     standardtag = Val{true}(), concrete_jac = nothing,
     diff_type = Val{:forward}, min_order = 5, max_order = 13, threading = false,
     linsolve = nothing, precs = DEFAULT_PRECS,
     extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
     new_W_γdt_cutoff = 1 // 5,
     controller = :Predictive, κ = nothing, maxiters = 10, smooth_est = true,
     step_limiter! = trivial_limiter!)

    AD_choice = _process_AD_choice(autodiff, chunk_size, diff_type)

 AdaptiveRadau{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
     typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
     typeof(κ), typeof(fast_convergence_cutoff),
     typeof(new_W_γdt_cutoff), typeof(step_limiter!), typeof(threading)}(linsolve,
     precs,
     smooth_est,
     extrapolant,
     κ,
     maxiters,
     fast_convergence_cutoff,
     new_W_γdt_cutoff,
     controller,
     step_limiter!, min_order, max_order, threading, 
     AD_choice)
end

