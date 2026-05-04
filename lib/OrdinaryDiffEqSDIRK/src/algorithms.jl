function SDIRK_docstring(
        description::String,
        name::String;
        references::String = "",
        extra_keyword_description::String = "",
        extra_keyword_default::String = ""
    )
    keyword_default = """
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing,
        nlsolve = NLNewton(),
        """ * extra_keyword_default

    keyword_default_description = """
        - `autodiff`: Uses [ADTypes.jl](https://sciml.github.io/ADTypes.jl/stable/) 
            to specify whether to use automatic differentiation via
            [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) or finite
            differencing via [FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl). 
            Defaults to `AutoForwardDiff()` for automatic differentiation, which by default uses
            `chunksize = 0`, and thus uses the internal ForwardDiff.jl algorithm for the choice.
            To use `FiniteDiff.jl`, the `AutoFiniteDiff()` ADType can be used, which has a keyword argument
            `fdtype` with default value `Val{:forward}()`, and alternatives `Val{:central}()` and `Val{:complex}()`.
        - `concrete_jac`: Specifies whether a Jacobian should be constructed. Defaults to
            `nothing`, which means it will be chosen true/false depending on circumstances
            of the solver, such as whether a Krylov subspace method is used for `linsolve`.
        - `linsolve`: Any [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) compatible linear solver.
          For example, to use [KLU.jl](https://github.com/JuliaSparse/KLU.jl), specify
          `$name(linsolve = KLUFactorization()`).
           When `nothing` is passed, uses `DefaultLinearSolver`.
        - `nlsolve`: nonlinear solver algorithm used for solving the implicit system.
            """ * extra_keyword_description

    return generic_solver_docstring(
        description, name, "SDIRK Method.", references,
        keyword_default_description, keyword_default
    )
end

abstract type OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm <: OrdinaryDiffEqNewtonAdaptiveAlgorithm end

@doc SDIRK_docstring(
    "A 1st order implicit solver. A-B-L-stable. Adaptive timestepping through a divided differences estimate. Strong-stability preserving (SSP). Good for highly stiff equations.",
    "ImplicitEuler";
    references = "@book{wanner1996solving,
    title={Solving ordinary differential equations II},
    author={Wanner, Gerhard and Hairer, Ernst},
    volume={375},
    year={1996},
    publisher={Springer Berlin Heidelberg New York}}",
    extra_keyword_description = """
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    extrapolant = :constant,
    step_limiter! = trivial_limiter!,
    """
)
struct ImplicitEuler{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end

function ImplicitEuler(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return ImplicitEuler(
        linsolve,
        nlsolve, extrapolant, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "A second order A-stable symplectic and symmetric implicit solver. Excellent for Hamiltonian systems and highly stiff equations.",
    "ImplicitMidpoint";
    references = "@book{wanner1996solving,
    title={Solving ordinary differential equations II},
    author={Wanner, Gerhard and Hairer, Ernst},
    volume={375},
    year={1996},
    publisher={Springer Berlin Heidelberg New York}}",
    extra_keyword_description = """
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct ImplicitMidpoint{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end

function ImplicitMidpoint(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear, step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return ImplicitMidpoint(
        linsolve,
        nlsolve,
        extrapolant,
        step_limiter!, autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc SDIRK_docstring(
    "A second order A-stable symmetric ESDIRK method. 'Almost symplectic' without numerical dampening.",
    "Trapezoid";
    references = "Andre Vladimirescu. 1994. The Spice Book. John Wiley & Sons, Inc., New York, NY, USA.",
    extra_keyword_description = """
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct Trapezoid{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end

function Trapezoid(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return Trapezoid(
        linsolve,
        nlsolve,
        extrapolant,
        step_limiter!,
        autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "A second order A-B-L-S-stable one-step ESDIRK method. Includes stiffness-robust error estimates for accurate adaptive timestepping, smoothed derivatives for highly stiff and oscillatory problems. Good for high tolerances (>1e-2) on stiff problems.",
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
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct TRBDF2{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end

function TRBDF2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return TRBDF2(
        linsolve, nlsolve,
        smooth_est, extrapolant, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@truncate_stacktrace TRBDF2

@doc SDIRK_docstring(
    "SDIRK2: SDIRK Method An A-B-L stable 2nd order SDIRK method",
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
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct SDIRK2{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end

function SDIRK2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return SDIRK2(
        linsolve, nlsolve, smooth_est, extrapolant,
        step_limiter!,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc SDIRK_docstring(
    "An SDIRK method.",
    "SDIRK22";
    references = "@techreport{kennedy2016diagonally,
    title={Diagonally implicit Runge-Kutta methods for ordinary differential equations. A review},
    author={Kennedy, Christopher A and Carpenter, Mark H},
    year={2016}}",
    extra_keyword_description = """
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct SDIRK22{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end

function SDIRK22(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return Trapezoid(
        linsolve,
        nlsolve,
        extrapolant,
        step_limiter!,
        autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    """SSPSDIRK is an SSP-optimized SDIRK method,
    so it's an implicit SDIRK method for handling stiffness but if the `dt` is below the SSP `coefficient * dt`,
    then the SSP property of the SSP integrators (the other page) is satisfied.
    As such this is a method which is expected to be good on advection-dominated cases where an explicit SSP integrator would be used,
    but where reaction equations are sufficient stiff to justify implicit integration.""",
    "SSPSDIRK2";
    references = "@article{ketcheson2009optimal,
    title={Optimal implicit strong stability preserving Runge--Kutta methods},
    author={Ketcheson, David I and Macdonald, Colin B and Gottlieb, Sigal},
    journal={Applied Numerical Mathematics},
    volume={59},
    number={2},
    pages={373--392},
    year={2009},
    publisher={Elsevier}}",
    extra_keyword_description = """
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
        """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :constant,
    """
)
struct SSPSDIRK2{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm # Not adaptive
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end

function SSPSDIRK2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :constant,
    )
    autodiff = _fixup_ad(autodiff)

    return SSPSDIRK2(
        linsolve, nlsolve, smooth_est, extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "An A-L stable stiffly-accurate 3rd order ESDIRK method.",
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
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct Kvaerno3{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function Kvaerno3(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return Kvaerno3(
        linsolve, nlsolve,
        smooth_est, extrapolant, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "An A-L stable stiffly-accurate 3rd order ESDIRK method with splitting.",
    "KenCarp3";
    references = "@book{kennedy2001additive,
    title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
    author={Kennedy, Christopher Alan},
    year={2001},
    publisher={National Aeronautics and Space Administration, Langley Research Center}}",
    extra_keyword_description = """
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct KenCarp3{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function KenCarp3(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return KenCarp3(
        linsolve, nlsolve,
        smooth_est, extrapolant, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "Third order method.",
    "CFNLIRK3";
    references = "@article{calvo2001linearly,
    title={Linearly implicit Runge--Kutta methods for advection--reaction--diffusion equations},
    author={Calvo, MP and De Frutos, J and Novo, J},
    journal={Applied Numerical Mathematics},
    volume={37},
    number={4},
    pages={535--549},
    year={2001},
    publisher={Elsevier}}",
    extra_keyword_description = """
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """
)
struct CFNLIRK3{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end
function CFNLIRK3(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear
    )
    autodiff = _fixup_ad(autodiff)

    return CFNLIRK3(
        linsolve,
        nlsolve,
        extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc SDIRK_docstring(
    "An A-L stable 4th order SDIRK method.",
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
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `embedding`: which embedded error estimate to use for step size control.
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    embedding = 3,
    """
)
struct Cash4{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    embedding::Int
    autodiff::AD
    concrete_jac::CJ
end
function Cash4(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        embedding = 3
    )
    autodiff = _fixup_ad(autodiff)

    return Cash4(
        linsolve,
        nlsolve,
        smooth_est,
        extrapolant,
        embedding,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc SDIRK_docstring(
    "Method of order 4.",
    "SFSDIRK4";
    references = "@article{ferracina2008strong,
    title={Strong stability of singly-diagonally-implicit Runge--Kutta methods},
    author={Ferracina, Luca and Spijker, MN},
    journal={Applied Numerical Mathematics},
    volume={58},
    number={11},
    pages={1675--1686},
    year={2008},
    publisher={Elsevier}}",
    extra_keyword_description = """
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """
)
struct SFSDIRK4{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end
function SFSDIRK4(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear
    )
    autodiff = _fixup_ad(autodiff)

    return SFSDIRK4(
        linsolve,
        nlsolve,
        extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc SDIRK_docstring(
    "Method of order 5.",
    "SFSDIRK5";
    references = "@article{ferracina2008strong,
    title={Strong stability of singly-diagonally-implicit Runge--Kutta methods},
    author={Ferracina, Luca and Spijker, MN},
    journal={Applied Numerical Mathematics},
    volume={58},
    number={11},
    pages={1675--1686},
    year={2008},
    publisher={Elsevier}}",
    extra_keyword_description = """
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """
)
struct SFSDIRK5{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end

function SFSDIRK5(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear
    )
    autodiff = _fixup_ad(autodiff)

    return SFSDIRK5(
        linsolve,
        nlsolve,
        extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc SDIRK_docstring(
    "Method of order 6.",
    "SFSDIRK6";
    references = "@article{ferracina2008strong,
    title={Strong stability of singly-diagonally-implicit Runge--Kutta methods},
    author={Ferracina, Luca and Spijker, MN},
    journal={Applied Numerical Mathematics},
    volume={58},
    number={11},
    pages={1675--1686},
    year={2008},
    publisher={Elsevier}}",
    extra_keyword_description = """
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """
)
struct SFSDIRK6{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end

function SFSDIRK6(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear
    )
    autodiff = _fixup_ad(autodiff)

    return SFSDIRK6(
        linsolve,
        nlsolve,
        extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc SDIRK_docstring(
    "Method of order 7.",
    "SFSDIRK7";
    references = "@article{ferracina2008strong,
    title={Strong stability of singly-diagonally-implicit Runge--Kutta methods},
    author={Ferracina, Luca and Spijker, MN},
    journal={Applied Numerical Mathematics},
    volume={58},
    number={11},
    pages={1675--1686},
    year={2008},
    publisher={Elsevier}}",
    extra_keyword_description = """
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """
)
struct SFSDIRK7{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end

function SFSDIRK7(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear
    )
    autodiff = _fixup_ad(autodiff)

    return SFSDIRK7(
        linsolve,
        nlsolve,
        extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc SDIRK_docstring(
    "Method of order 8.",
    "SFSDIRK8";
    references = "@article{ferracina2008strong,
    title={Strong stability of singly-diagonally-implicit Runge--Kutta methods},
    author={Ferracina, Luca and Spijker, MN},
    journal={Applied Numerical Mathematics},
    volume={58},
    number={11},
    pages={1675--1686},
    year={2008},
    publisher={Elsevier}}",
    extra_keyword_description = """
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    """,
    extra_keyword_default = """
    extrapolant = :linear,
    """
)
struct SFSDIRK8{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end

function SFSDIRK8(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear
    )
    autodiff = _fixup_ad(autodiff)

    return SFSDIRK8(
        linsolve,
        nlsolve,
        extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc SDIRK_docstring(
    "An A-L stable 4th order SDIRK method.",
    "Hairer4";
    references = "E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
    differential-algebraic problems. Computational mathematics (2nd revised ed.),
    Springer (1996)",
    extra_keyword_description = """
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
        """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    """
)
struct Hairer4{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end
function Hairer4(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
    )
    autodiff = _fixup_ad(autodiff)

    return Hairer4(
        linsolve, nlsolve, smooth_est, extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "An A-L stable 4th order SDIRK method.",
    "Hairer42";
    references = "E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
    differential-algebraic problems. Computational mathematics (2nd revised ed.),
    Springer (1996)",
    extra_keyword_description = """
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
        """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    """
)
struct Hairer42{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end
function Hairer42(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
    )
    autodiff = _fixup_ad(autodiff)

    return Hairer42(
        linsolve, nlsolve, smooth_est, extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "An A-L stable stiffly-accurate 4th order ESDIRK method.",
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
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct Kvaerno4{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function Kvaerno4(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return Kvaerno4(
        linsolve, nlsolve,
        smooth_est, extrapolant, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "An A-L stable stiffly-accurate 5th order ESDIRK method.",
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
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct Kvaerno5{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function Kvaerno5(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return Kvaerno5(
        linsolve, nlsolve,
        smooth_est, extrapolant, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "An A-L stable stiffly-accurate 4th order ESDIRK method with splitting. Includes splitting capabilities. Recommended for medium tolerance stiff problems (>1e-8).",
    "KenCarp4";
    references = "@book{kennedy2001additive,
    title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
    author={Kennedy, Christopher Alan},
    year={2001},
    publisher={National Aeronautics and Space Administration, Langley Research Center}}",
    extra_keyword_description = """
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct KenCarp4{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function KenCarp4(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return KenCarp4(
        linsolve, nlsolve,
        smooth_est, extrapolant, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@truncate_stacktrace KenCarp4

@doc SDIRK_docstring(
    "An A-L stable stiffly-accurate 4th order seven-stage ESDIRK method with splitting.",
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
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
        """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    """
)
struct KenCarp47{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function KenCarp47(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return KenCarp47(
        linsolve, nlsolve, smooth_est, extrapolant,
        controller, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "An A-L stable stiffly-accurate 5th order ESDIRK method with splitting.",
    "KenCarp5";
    references = "@book{kennedy2001additive,
    title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
    author={Kennedy, Christopher Alan},
    year={2001},
    publisher={National Aeronautics and Space Administration, Langley Research Center}}",
    extra_keyword_description = """
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct KenCarp5{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function KenCarp5(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return KenCarp5(
        linsolve, nlsolve,
        smooth_est, extrapolant, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "An A-L stable stiffly-accurate 5th order eight-stage ESDIRK method with splitting.",
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
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
        """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    """
)
struct KenCarp58{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function KenCarp58(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return KenCarp58(
        linsolve, nlsolve, smooth_est, extrapolant,
        controller, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

# `smooth_est` is not necessary, as the embedded method is also L-stable
@doc SDIRK_docstring(
    "Optimized ESDIRK tableaus.
Updates of the original KenCarp tableau expected to achieve lower error for the same steps in theory,
but are still being fully evaluated in context.",
    "ESDIRK54I8L2SA";
    references = """@article{Kennedy2019DiagonallyIR,
    title={Diagonally implicit Runge–Kutta methods for stiff ODEs},
    author={Christopher A. Kennedy and Mark H. Carpenter},
    journal={Applied Numerical Mathematics},
    year={2019},
    volume={146},
    pages={221-244}
    }""",
    extra_keyword_description = """
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
        """,
    extra_keyword_default = """
    extrapolant = :linear,
    """
)
struct ESDIRK54I8L2SA{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function ESDIRK54I8L2SA(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return ESDIRK54I8L2SA(
        linsolve, nlsolve, smooth_est, extrapolant,
        controller, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "Optimized ESDIRK tableaus.
Updates of the original KenCarp tableau expected to achieve lower error for the same steps in theory,
but are still being fully evaluated in context.",
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
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
        """,
    extra_keyword_default = """
    extrapolant = :linear,
    """
)
struct ESDIRK436L2SA2{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function ESDIRK436L2SA2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return ESDIRK436L2SA2(
        linsolve, nlsolve, smooth_est, extrapolant,
        controller, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "Optimized ESDIRK tableaus.
Updates of the original KenCarp tableau expected to achieve lower error for the same steps in theory,
but are still being fully evaluated in context.",
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
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
        """,
    extra_keyword_default = """
    extrapolant = :linear,
    """
)
struct ESDIRK437L2SA{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function ESDIRK437L2SA(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return ESDIRK437L2SA(
        linsolve, nlsolve, smooth_est, extrapolant,
        controller, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "Optimized ESDIRK tableaus.
Updates of the original KenCarp tableau expected to achieve lower error for the same steps in theory,
but are still being fully evaluated in context.",
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
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
        """,
    extra_keyword_default = """
    extrapolant = :linear,
    """
)
struct ESDIRK547L2SA2{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function ESDIRK547L2SA2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return ESDIRK547L2SA2(
        linsolve, nlsolve, smooth_est, extrapolant,
        controller, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "Optimized ESDIRK tableaus.
Updates of the original KenCarp tableau expected to achieve lower error for the same steps in theory,
but are still being fully evaluated in context.
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
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
        """,
    extra_keyword_default = """
    extrapolant = :linear,
    """
)
struct ESDIRK659L2SA{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function ESDIRK659L2SA(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return ESDIRK659L2SA(
        linsolve, nlsolve, smooth_est, extrapolant,
        controller, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc SDIRK_docstring(
    "3rd order L-stable IMEX ARK method. Uses a generic tableau-driven implementation that supports both split and non-split forms.",
    "ARS343";
    references = "@article{ascher1997implicit,
    title={Implicit-explicit Runge-Kutta methods for time-dependent partial differential equations},
    author={Ascher, Uri M and Ruuth, Steven J and Spiteri, Raymond J},
    journal={Applied Numerical Mathematics},
    volume={25},
    number={2-3},
    pages={151--167},
    year={1997},
    publisher={Elsevier}}",
    extra_keyword_description = """
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    smooth_est = true,
    extrapolant = :linear,
    controller = :PI,
    step_limiter! = trivial_limiter!,
    """
)
struct ARS343{AD, F, F2, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveESDIRKAlgorithm
    linsolve::F
    nlsolve::F2
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function ARS343(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return ARS343(
        linsolve, nlsolve, smooth_est, extrapolant,
        controller, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end
