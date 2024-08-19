function SDIRK_docstring(description::String,
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
        nlsolve = NLNewton(),
        """ * extra_keyword_default

    keyword_default_description = """
        - `chunk_size`: TBD
        - `autodiff`: TBD
        - `standardtag`: TBD
        - `concrete_jac`: TBD
        - `diff_type`: TBD
        - `linsolve`: TBD
        - `precs`: TBD
        - `nlsolve`: TBD
        """ * extra_keyword_description

    generic_solver_docstring(
        description, name, "SDIRK Method.", references,
        keyword_default_description, keyword_default
    )
end

@doc SDIRK_docstring("A 1st order implicit solver. A-B-L-stable.
    Adaptive timestepping through a divided differences estimate via memory.
    Strong-stability preserving (SSP).",
    "ImplicitEuler";
    references = "TBD",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :constant,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """)
struct ImplicitEuler{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end

function ImplicitEuler(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant,
        controller = :PI, step_limiter! = trivial_limiter!)
    ImplicitEuler{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve, precs, extrapolant, controller, step_limiter!)
end

@doc SDIRK_docstring("A second order A-stable symplectic and symmetric implicit solver.
    Good for highly stiff equations which need symplectic integration.",
    "ImplicitMidpoint";
    references = "TBD",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `step_limiter!`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """)
struct ImplicitMidpoint{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    step_limiter!::StepLimiter
end

function ImplicitMidpoint(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, step_limiter! = trivial_limiter!)
    ImplicitMidpoint{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        step_limiter!)
end


@doc SDIRK_docstring("""Second order A-stable symmetric ESDIRK method.
    "Almost symplectic" without numerical dampening.
    Also known as Crank-Nicolson when applied to PDEs. Adaptive timestepping via divided
    differences approximation to the second derivative terms in the local truncation error
    estimate (the SPICE approximation strategy).""",
    "Trapezoid";
    references = "Andre Vladimirescu. 1994. The Spice Book. John Wiley & Sons, Inc., New York, NY, USA.",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """)
struct Trapezoid{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end

function Trapezoid(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    Trapezoid{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        controller,
        step_limiter!)
end

@doc SDIRK_docstring("A second order A-B-L-S-stable one-step ESDIRK method.
    Includes stiffness-robust error estimates for accurate adaptive timestepping,
    smoothed derivatives for highly stiff and oscillatory problems.",
    "TRBDF2";
    references = "@article{hosea1996analysis,
    title={Analysis and implementation of TR-BDF2},
    author={Hosea, ME and Shampine, LF},
    journal={Applied Numerical Mathematics},
    volume={20},
    number={1-2},
    pages={21--37},
    year={1996},
    publisher={Elsevier}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """)
struct TRBDF2{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end

function TRBDF2(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    TRBDF2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

TruncatedStacktraces.@truncate_stacktrace TRBDF2

@doc SDIRK_docstring("SDIRK2: SDIRK Method An A-B-L stable 2nd order SDIRK method",
    "SDIRK2";
    references = "@article{hindmarsh2005sundials,
    title={{SUNDIALS}: Suite of nonlinear and differential/algebraic equation solvers},
    author={Hindmarsh, Alan C and Brown, Peter N and Grant, Keith E and Lee, Steven L and Serban, Radu and Shumaker, Dan E and Woodward, Carol S},
    journal={ACM Transactions on Mathematical Software (TOMS)},
    volume={31},
    number={3},
    pages={363--396},
    year={2005},
    publisher={ACM}}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """)
struct SDIRK2{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end

function SDIRK2(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    SDIRK2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(
        linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller,
        step_limiter!)
end

@doc SDIRK_docstring("Description TBD",
    "SDIRK22";
    references = "TBD",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """)
struct SDIRK22{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end

function SDIRK22(;
        chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    Trapezoid{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        controller,
        step_limiter!)
end

@doc SDIRK_docstring("Description TBD",
    "SSPSDIRK2";
    references = "TBD",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :constant,
    controller = :PI,
    """)
struct SSPSDIRK2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ} # Not adaptive
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
end

function SSPSDIRK2(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :constant,
        controller = :PI)
    SSPSDIRK2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller)
end


@doc SDIRK_docstring("An A-L stable stiffly-accurate 3rd order ESDIRK method.",
    "Kvaerno3";
    references = "@article{kvaerno2004singly,
    title={Singly diagonally implicit Runge--Kutta methods with an explicit first stage},
    author={Kv{\\ae}rn{\\o}, Anne},
    journal={BIT Numerical Mathematics},
    volume={44},
    number={3},
    pages={489--502},
    year={2004},
    publisher={Springer}}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """)
struct Kvaerno3{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function Kvaerno3(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    Kvaerno3{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

@doc SDIRK_docstring("An A-L stable stiffly-accurate 3rd order ESDIRK method with splitting.",
    "KenCarp3";
    references = "@book{kennedy2001additive,
    title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
    author={Kennedy, Christopher Alan},
    year={2001},
    publisher={National Aeronautics and Space Administration, Langley Research Center}}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """)
struct KenCarp3{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function KenCarp3(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    KenCarp3{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

@doc SDIRK_docstring("Description TBD.",
    "CFNLIRK3";
    references = "TBD",
    extra_keyword_description = """
    - `extrapolant`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """)
struct CFNLIRK3{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end
function CFNLIRK3(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    CFNLIRK3{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

@doc SDIRK_docstring("An A-L stable 4th order SDIRK method.",
    "Cash4";
    references = "@article{hindmarsh2005sundials,
    title={{SUNDIALS}: Suite of nonlinear and differential/algebraic equation solvers},
    author={Hindmarsh, Alan C and Brown, Peter N and Grant, Keith E and Lee, Steven L and Serban, Radu and Shumaker, Dan E and Woodward, Carol S},
    journal={ACM Transactions on Mathematical Software (TOMS)},
    volume={31},
    number={3},
    pages={363--396},
    year={2005},
    publisher={ACM}}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `embedding`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    embedding = 3,
    """)
struct Cash4{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    embedding::Int
    controller::Symbol
end
function Cash4(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, embedding = 3)
    Cash4{
        _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(
        linsolve,
        nlsolve,
        precs,
        smooth_est,
        extrapolant,
        embedding,
        controller)
end

@doc SDIRK_docstring("Description TBD.",
    "SFSDIRK4";
    references = "TBD",
    extra_keyword_description = """
    - `extrapolant`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """)
struct SFSDIRK4{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end
function SFSDIRK4(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    SFSDIRK4{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

@doc SDIRK_docstring("Description TBD.",
    "SFSDIRK5";
    references = "TBD",
    extra_keyword_description = """
    - `extrapolant`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """)
struct SFSDIRK5{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end

function SFSDIRK5(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    SFSDIRK5{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

@doc SDIRK_docstring("Description TBD.",
    "SFSDIRK6";
    references = "TBD",
    extra_keyword_description = """
    - `extrapolant`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """)
struct SFSDIRK6{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end

function SFSDIRK6(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    SFSDIRK6{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end


@doc SDIRK_docstring("Description TBD.",
    "SFSDIRK7";
    references = "TBD",
    extra_keyword_description = """
    - `extrapolant`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """)
struct SFSDIRK7{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end

function SFSDIRK7(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    SFSDIRK7{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

@doc SDIRK_docstring("Description TBD.",
    "SFSDIRK8";
    references = "TBD",
    extra_keyword_description = """
    - `extrapolant`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """)
struct SFSDIRK8{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end

function SFSDIRK8(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    SFSDIRK8{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

@doc SDIRK_docstring("An A-L stable 4th order SDIRK method.",
    "Hairer4";
    references = "E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
    differential-algebraic problems. Computational mathematics (2nd revised ed.),
    Springer (1996)",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    """)
struct Hairer4{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
end
function Hairer4(;
        chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    Hairer4{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller)
end

@doc SDIRK_docstring("An A-L stable 4th order SDIRK method.",
    "Hairer42";
    references = "E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
    differential-algebraic problems. Computational mathematics (2nd revised ed.),
    Springer (1996)",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    """)
struct Hairer42{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
end
function Hairer42(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    Hairer42{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller)
end

@doc SDIRK_docstring("An A-L stable stiffly-accurate 4th order ESDIRK method.",
    "Kvaerno4";
    references = "@article{kvaerno2004singly,
    title={Singly diagonally implicit Runge--Kutta methods with an explicit first stage},
    author={Kv{\\ae}rn{\\o}, Anne},
    journal={BIT Numerical Mathematics},
    volume={44},
    number={3},
    pages={489--502},
    year={2004},
    publisher={Springer}}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """)
struct Kvaerno4{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function Kvaerno4(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    Kvaerno4{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

@doc SDIRK_docstring("An A-L stable stiffly-accurate 5th order ESDIRK method.",
    "Kvaerno5";
    references = "@article{kvaerno2004singly,
    title={Singly diagonally implicit Runge--Kutta methods with an explicit first stage},
    author={Kv{\\ae}rn{\\o}, Anne},
    journal={BIT Numerical Mathematics},
    volume={44},
    number={3},
    pages={489--502},
    year={2004},
    publisher={Springer}}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """)
struct Kvaerno5{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function Kvaerno5(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    Kvaerno5{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

@doc SDIRK_docstring("An A-L stable stiffly-accurate 4th order ESDIRK method with splitting.",
    "KenCarp4";
    references = "@book{kennedy2001additive,
    title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
    author={Kennedy, Christopher Alan},
    year={2001},
    publisher={National Aeronautics and Space Administration, Langley Research Center}}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """)
struct KenCarp4{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function KenCarp4(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    KenCarp4{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

TruncatedStacktraces.@truncate_stacktrace KenCarp4

@doc SDIRK_docstring("An A-L stable stiffly-accurate 4th order seven-stage ESDIRK method with splitting.",
    "KenCarp47";
    references = "@article{kennedy2019higher,
    title={Higher-order additive Runge--Kutta schemes for ordinary differential equations},
    author={Kennedy, Christopher A and Carpenter, Mark H},
    journal={Applied Numerical Mathematics},
    volume={136},
    pages={183--205},
    year={2019},
    publisher={Elsevier}}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    """)
struct KenCarp47{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
end
function KenCarp47(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    KenCarp47{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller)
end

@doc SDIRK_docstring("An A-L stable stiffly-accurate 5th order ESDIRK method with splitting.",
    "KenCarp5";
    references = "@book{kennedy2001additive,
    title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
    author={Kennedy, Christopher Alan},
    year={2001},
    publisher={National Aeronautics and Space Administration, Langley Research Center}}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """)
struct KenCarp5{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function KenCarp5(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    KenCarp5{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

@doc SDIRK_docstring("An A-L stable stiffly-accurate 5th order eight-stage ESDIRK method with splitting.",
    "KenCarp58";
    references = "@article{kennedy2019higher,
    title={Higher-order additive Runge--Kutta schemes for ordinary differential equations},
    author={Kennedy, Christopher A and Carpenter, Mark H},
    journal={Applied Numerical Mathematics},
    volume={136},
    pages={183--205},
    year={2019},
    publisher={Elsevier}}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    """)
struct KenCarp58{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
end
function KenCarp58(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    KenCarp58{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller)
end

# `smooth_est` is not necessary, as the embedded method is also L-stable
@doc SDIRK_docstring("Description TBD.",
    "ESDIRK54I8L2SA";
    references = "TBD",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    controller = :PI,
    """)
struct ESDIRK54I8L2SA{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function ESDIRK54I8L2SA(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    ESDIRK54I8L2SA{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller)
end

@doc SDIRK_docstring("Description TBD.",
    "ESDIRK436L2SA2";
    references = """@article{Kennedy2019DiagonallyIR,
    title={Diagonally implicit Runge–Kutta methods for stiff ODEs},
    author={Christopher A. Kennedy and Mark H. Carpenter},
    journal={Applied Numerical Mathematics},
    year={2019},
    volume={146},
    pages={221-244}
    }""",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    controller = :PI,
    """)
struct ESDIRK436L2SA2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function ESDIRK436L2SA2(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    ESDIRK436L2SA2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller)
end

@doc SDIRK_docstring("Description TBD.",
    "ESDIRK437L2SA";
    references = """@article{Kennedy2019DiagonallyIR,
    title={Diagonally implicit Runge–Kutta methods for stiff ODEs},
    author={Christopher A. Kennedy and Mark H. Carpenter},
    journal={Applied Numerical Mathematics},
    year={2019},
    volume={146},
    pages={221-244}
    }""",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    controller = :PI,
    """)
struct ESDIRK437L2SA{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function ESDIRK437L2SA(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    ESDIRK437L2SA{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller)
end

@doc SDIRK_docstring("Description TBD.",
    "ESDIRK547L2SA2";
    references = """@article{Kennedy2019DiagonallyIR,
    title={Diagonally implicit Runge–Kutta methods for stiff ODEs},
    author={Christopher A. Kennedy and Mark H. Carpenter},
    journal={Applied Numerical Mathematics},
    year={2019},
    volume={146},
    pages={221-244}
    }""",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    controller = :PI,
    """)
struct ESDIRK547L2SA2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function ESDIRK547L2SA2(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    ESDIRK547L2SA2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller)
end

@doc SDIRK_docstring("Description TBD.
    Currently has STABILITY ISSUES, causing it to fail the adaptive tests.
    Check issue https://github.com/SciML/OrdinaryDiffEq.jl/issues/1933 for more details.",
    "ESDIRK659L2SA";
    references = """@article{Kennedy2019DiagonallyIR,
    title={Diagonally implicit Runge–Kutta methods for stiff ODEs},
    author={Christopher A. Kennedy and Mark H. Carpenter},
    journal={Applied Numerical Mathematics},
    year={2019},
    volume={146},
    pages={221-244}
    }""",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    controller = :PI,
    """)
struct ESDIRK659L2SA{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function ESDIRK659L2SA(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    ESDIRK659L2SA{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller)
end
