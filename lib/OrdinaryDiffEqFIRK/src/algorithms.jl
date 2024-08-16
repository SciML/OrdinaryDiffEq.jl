function implicit_rk_docstring(description::String,
    name::String;
    references::String = "",
    extra_keyword_description::String = "",
    extra_keyword_default::String = "")
    keyword_default = """
    chunk_size = Val{0}(),
    autodiff = true,
    standardtag = Val{true}(),
    concrete_jac = nothing,
    diff_type = Val{:forward},
    linsolve = nothing,
    precs = DEFAULT_PRECS,
    """ * extra_keyword_default

    keyword_default_description = """
    - `chunk_size`: TBD
    - `autodiff`: TBD
    - `standardtag`: TBD
    - `concrete_jac`: TBD
    - `diff_type`: TBD
    - `linsolve`: TBD
    - `precs`: TBD
    """ * extra_keyword_description

    generic_solver_docstring(
        description, name, "Fully-Implicit Runge-Kutta Method.", references,
        keyword_default_description, keyword_default
    )
end
# FIRK Methods

"""
@article{hairer1999stiff,
title={Stiff differential equations solved by Radau methods},
author={Hairer, Ernst and Wanner, Gerhard},
journal={Journal of Computational and Applied Mathematics},
volume={111},
number={1-2},
pages={93--111},
year={1999},
publisher={Elsevier}
}

RadauIIA3: Fully-Implicit Runge-Kutta Method
An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency.
"""
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
end

function RadauIIA3(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10,
        step_limiter! = trivial_limiter!)
    RadauIIA3{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
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
        step_limiter!)
end

"""
@article{hairer1999stiff,
title={Stiff differential equations solved by Radau methods},
author={Hairer, Ernst and Wanner, Gerhard},
journal={Journal of Computational and Applied Mathematics},
volume={111},
number={1-2},
pages={93--111},
year={1999},
publisher={Elsevier}
}

RadauIIA5: Fully-Implicit Runge-Kutta Method
An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency.
"""
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
end

function RadauIIA5(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10, smooth_est = true,
        step_limiter! = trivial_limiter!)
    RadauIIA5{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
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
        step_limiter!)
end

@doc implicit_rk_docstring("An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency.
    Similar to Hairer's SEULEX.",
    "RadauIIA9";
    references = """@article{hairer1999stiff,
    title={Stiff differential equations solved by Radau methods},
    author={Hairer, Ernst and Wanner, Gerhard},
    journal={Journal of Computational and Applied Mathematics},
    volume={111},
    number={1-2},
    pages={93--111},
    year={1999},
    publisher={Elsevier}}""",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `fast_convergence_cutoff`: TBD
    - `new_W_γdt_cutoff`: TBD
    - `controller`: TBD
    - `κ`: TBD
    - `maxiters`: TBD
    - `smooth_est`: TBD
    - `step_limiter!`: TBD""",
    extra_keyword_default = """
    extrapolant = :dense,
    fast_convergence_cutoff = 1 // 5,
    new_W_γdt_cutoff = 1 // 5,
    controller = :Predictive,
    κ = nothing,
    maxiters = 10,
    smooth_est = true,
    step_limiter! = trivial_limiter!
    """)
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
end

function RadauIIA9(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10, smooth_est = true,
        step_limiter! = trivial_limiter!)
    RadauIIA9{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
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
        step_limiter!)
end
