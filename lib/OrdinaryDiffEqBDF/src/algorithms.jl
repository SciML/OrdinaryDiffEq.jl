"""Helper to generate docstrings for BDF-family algorithms."""
function BDF_docstring(
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
        """ * "\n" * extra_keyword_default

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
            """ * "/n" * extra_keyword_description
    return generic_solver_docstring(
        description, name, "Multistep Method.", references,
        keyword_default_description, keyword_default
    )
end

@doc BDF_docstring(
    "An adaptive order 2 L-stable fixed leading coefficient multistep BDF method.",
    "ABDF2",
    references = """
    E. Alberdi Celayaa, J. J. Anza Aguirrezabalab, P. Chatzipantelidisc. Implementation of
    an Adaptive BDF2 Formula and Comparison with The MATLAB Ode15s. Procedia Computer Science,
    29, pp 1014-1026, 2014. doi: https://doi.org/10.1016/j.procs.2014.05.091
    """,
    extra_keyword_description = """
    - `κ`: coefficient for the order and stability control of the BDF method. When `nothing`, the default value is used.
    - `tol`: tolerance for the nonlinear solver. When `nothing`, uses the default tolerance.
    - `nlsolve`: nonlinear solver algorithm used for solving the implicit system.
    - `smooth_est`: whether to use a smoothed estimate for error control.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    κ = nothing,
    tol = nothing,
    nlsolve = NLNewton(),
    smooth_est = true,
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    """
)
struct ABDF2{AD, F, F2, K, T, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    κ::K
    tol::T
    smooth_est::Bool
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end
function ABDF2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        κ = nothing, tol = nothing, linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return ABDF2(
        linsolve, nlsolve, κ, tol,
        smooth_est, extrapolant, step_limiter!, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc BDF_docstring(
    "Implicit-explicit (IMEX) method designed for SplitODEFunction equations,
which reduce the size of the implicit handling to a subset of the equations.
It's similar to the additive Runge-Kutta methods in splitting mode,
like `KenCarp4`, but instead using a multistep BDF approach",
    "SBDF",
    references = """@article{ascher1995implicit,
    title={Implicit-explicit methods for time-dependent partial differential equations},
    author={Ascher, Uri M and Ruuth, Steven J and Wetton, Brian TR},
    journal={SIAM Journal on Numerical Analysis},
    volume={32},
    number={3},
    pages={797--823},
    year={1995},
    publisher={SIAM}}
    """,
    extra_keyword_description = """
    - `κ`: coefficient for the order and stability control of the BDF method. When `nothing`, the default value is used.
    - `tol`: tolerance for the nonlinear solver. When `nothing`, uses the default tolerance.
    - `nlsolve`: nonlinear solver algorithm used for solving the implicit system.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `ark`: whether to use an additive Runge-Kutta formulation.
    """,
    extra_keyword_default = """
    κ = nothing,
    tol = nothing,
    nlsolve = NLNewton(),
    extrapolant = :linear,
    ark = false,
    order,
    """
)
struct SBDF{AD, F, F2, K, T, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    κ::K
    tol::T
    extrapolant::Symbol
    order::Int
    ark::Bool
    autodiff::AD
    concrete_jac::CJ
end

function SBDF(
        order; autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, ark = false
    )
    autodiff = _fixup_ad(autodiff)

    return SBDF(
        linsolve,
        nlsolve,
        κ,
        tol,
        extrapolant,
        order,
        ark,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

# All keyword form needed for remake
function SBDF(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear,
        order, ark = false
    )
    autodiff = _fixup_ad(autodiff)

    return SBDF(
        linsolve,
        nlsolve,
        κ,
        tol,
        extrapolant,
        order,
        ark,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

"""
    SBDF2(;kwargs...)

The two-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

See also `SBDF`.
"""
SBDF2(; kwargs...) = SBDF(2; kwargs...)

"""
    SBDF3(;kwargs...)

The three-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

See also `SBDF`.
"""
SBDF3(; kwargs...) = SBDF(3; kwargs...)

"""
    SBDF4(;kwargs...)

The four-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

See also `SBDF`.
"""
SBDF4(; kwargs...) = SBDF(4; kwargs...)

@doc BDF_docstring(
    "An adaptive order 1 quasi-constant timestep L-stable numerical differentiation function method.",
    "QNDF1",
    references = """@article{shampine1997matlab,
    title={The matlab ode suite},
    author={Shampine, Lawrence F and Reichelt, Mark W},
    journal={SIAM journal on scientific computing},
    volume={18},
    number={1},
    pages={1--22},
    year={1997},
    publisher={SIAM}
    }""",
    extra_keyword_description = """
    - `nlsolve`: nonlinear solver algorithm used for solving the implicit system.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `kappa`: coefficient for the BDF error estimator.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    nlsolve = NLNewton(),
    extrapolant = :linear,
    kappa = -0.1850,
    step_limiter! = trivial_limiter!,
    """
)
struct QNDF1{AD, F, F2, κType, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    kappa::κType
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end

function QNDF1(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear, kappa = -37 // 200,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return QNDF1(
        linsolve,
        nlsolve,
        extrapolant,
        kappa,
        step_limiter!,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc BDF_docstring(
    "An adaptive order 2 quasi-constant timestep L-stable numerical differentiation function (NDF) method.",
    "QNDF2",
    references = """@article{shampine1997matlab,
    title={The matlab ode suite},
    author={Shampine, Lawrence F and Reichelt, Mark W},
    journal={SIAM journal on scientific computing},
    volume={18},
    number={1},
    pages={1--22},
    year={1997},
    publisher={SIAM}
    }""",
    extra_keyword_description = """
    - `nlsolve`: nonlinear solver algorithm used for solving the implicit system.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `kappa`: coefficient for the BDF error estimator.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    nlsolve = NLNewton(),
    extrapolant = :linear,
    kappa =  -1 // 9,
    step_limiter! = trivial_limiter!,
    """
)
struct QNDF2{AD, F, F2, κType, StepLimiter, CJ} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    kappa::κType
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
end

function QNDF2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear, kappa = -1 // 9,
        step_limiter! = trivial_limiter!
    )
    autodiff = _fixup_ad(autodiff)

    return QNDF2(
        linsolve,
        nlsolve,
        extrapolant,
        kappa,
        step_limiter!,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc BDF_docstring(
    "An adaptive order quasi-constant timestep NDF method. Similar to MATLAB's ode15s. Uses Shampine's accuracy-optimal coefficients. Performance improves with larger, more complex ODEs. Good for medium to highly stiff problems. Recommended for large systems (>1000 ODEs).",
    "QNDF",
    references = """@article{shampine1997matlab,
    title={The matlab ode suite},
    author={Shampine, Lawrence F and Reichelt, Mark W},
    journal={SIAM journal on scientific computing},
    volume={18},
    number={1},
    pages={1--22},
    year={1997},
    publisher={SIAM}
    }""",
    extra_keyword_description = """
    - `κ`: coefficient for the order and stability control of the BDF method. When `nothing`, the default value is used.
    - `tol`: tolerance for the nonlinear solver. When `nothing`, uses the default tolerance.
    - `nlsolve`: nonlinear solver algorithm used for solving the implicit system.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `kappa`: coefficient for the BDF error estimator.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    κ = nothing,
    tol = nothing,
    nlsolve = NLNewton(),
    extrapolant = :linear,
    kappa =  promote(-0.1850, -1 // 9, -0.0823, -0.0415, 0),
    step_limiter! = trivial_limiter!,
    """
)
struct QNDF{MO, AD, F, F2, K, T, κType, StepLimiter, CJ, QT} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    κ::K
    tol::T
    extrapolant::Symbol
    kappa::κType
    step_limiter!::StepLimiter
    autodiff::AD
    concrete_jac::CJ
    qmax::QT
    qsteady_min::QT
    qsteady_max::QT
end

function QNDF(;
        max_order::Val{MO} = Val{5}(),
        autodiff = AutoForwardDiff(), concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, kappa = (
            -37 // 200, -1 // 9, -823 // 10000, -83 // 2000, 0 // 1,
        ),
        step_limiter! = trivial_limiter!,
        qsteady_min = 1 // 1, qsteady_max = 2 // 1, qmax = 10 // 1,
    ) where {MO}
    autodiff = _fixup_ad(autodiff)

    return QNDF(
        max_order, linsolve, nlsolve, κ, tol,
        extrapolant, kappa, step_limiter!, autodiff,
        _unwrap_val(concrete_jac), qmax, qsteady_min, qsteady_max,
    )
end

@truncate_stacktrace QNDF

@doc BDF_docstring(
    "The second order Modified Extended BDF method,
    which has improved stability properties over the standard BDF.
    Fixed timestep only.",
    "MEBDF2",
    references = """@article{cash2000modified,
    title={Modified extended backward differentiation formulae for the numerical solution of stiff initial value problems in ODEs and DAEs},
    author={Cash, JR},
    journal={Journal of Computational and Applied Mathematics},
    volume={125},
    number={1-2},
    pages={117--130},
    year={2000},
    publisher={Elsevier}}""",
    extra_keyword_description = """
    - `nlsolve`: nonlinear solver algorithm used for solving the implicit system.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    """,
    extra_keyword_default = """
    nlsolve = NLNewton(),
    extrapolant = :constant,
    """
)
struct MEBDF2{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end
function MEBDF2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant
    )
    autodiff = _fixup_ad(autodiff)

    return MEBDF2(
        linsolve,
        nlsolve,
        extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

@doc BDF_docstring(
    "An adaptive order quasi-constant timestep NDF method.
Fixed leading coefficient BDF.
Utilizes Shampine's accuracy-optimal kappa values as defaults (has a keyword argument for a tuple of kappa coefficients).",
    "FBDF",
    references = """@article{shampine2002solving,
    title={Solving 0= F (t, y (t), y′(t)) in Matlab},
    author={Shampine, Lawrence F},
    year={2002},
    publisher={Walter de Gruyter GmbH \\& Co. KG}}""",
    extra_keyword_description = """
    - `κ`: coefficient for the order and stability control of the BDF method. When `nothing`, the default value is used.
    - `tol`: tolerance for the nonlinear solver. When `nothing`, uses the default tolerance.
    - `nlsolve`: nonlinear solver algorithm used for solving the implicit system.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    - `max_order`: maximum order of the adaptive-order BDF method.
    - `stald`: Enable Stability Limit Detection (STALD) for BDF orders 3-5. Default: `true`.
    - `stald_rrcut`: STALD cutoff for characteristic root magnitude. Default: `0.98`.
    - `stald_vrrtol`: STALD tolerance for variance of ratios. Default: `1e-4`.
    - `stald_vrrt2`: STALD secondary variance tolerance. Default: `5e-4`.
    - `stald_sqtol`: STALD tolerance for quartic residual. Default: `1e-3`.
    - `stald_rrtol`: STALD tolerance for rr cross-verification. Default: `1e-2`.
    - `stald_tiny`: STALD tiny value to avoid division by zero. Default: `1e-90`.
    """,
    extra_keyword_default = """
    κ = nothing,
    tol = nothing,
    nlsolve = NLNewton(),
    extrapolant = :linear,
    step_limiter! = trivial_limiter!,
    max_order::Val{MO} = Val{5}(),
    stald = true,
    stald_rrcut = 0.98,
    stald_vrrtol = 1e-4,
    stald_vrrt2 = 5e-4,
    stald_sqtol = 1e-3,
    stald_rrtol = 1e-2,
    stald_tiny = 1e-90,
    """
)
struct FBDF{MO, AD, F, F2, K, T, StepLimiter, CJ, QT} <:
    OrdinaryDiffEqNewtonAdaptiveAlgorithm
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    κ::K
    tol::T
    extrapolant::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
    stald::Bool
    stald_rrcut::Float64
    stald_vrrtol::Float64
    stald_vrrt2::Float64
    stald_sqtol::Float64
    stald_rrtol::Float64
    stald_tiny::Float64
    concrete_jac::CJ
    qmax::QT
    qsteady_min::QT
    qsteady_max::QT
end

function FBDF(;
        max_order::Val{MO} = Val{5}(),
        autodiff = AutoForwardDiff(), concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, step_limiter! = trivial_limiter!,
        stald = true,
        stald_rrcut = 0.98,
        stald_vrrtol = 1.0e-4,
        stald_vrrt2 = 5.0e-4,
        stald_sqtol = 1.0e-3,
        stald_rrtol = 1.0e-2,
        stald_tiny = 1.0e-90,
        qsteady_min = 9 // 10, qsteady_max = 2 // 1, qmax = 10 // 1,
    ) where {MO}
    autodiff = _fixup_ad(autodiff)

    return FBDF(
        max_order, linsolve, nlsolve, κ, tol, extrapolant,
        step_limiter!, autodiff,
        stald, Float64(stald_rrcut), Float64(stald_vrrtol), Float64(stald_vrrt2),
        Float64(stald_sqtol), Float64(stald_rrtol), Float64(stald_tiny),
        _unwrap_val(concrete_jac),
        qmax, qsteady_min, qsteady_max,
    )
end

@truncate_stacktrace FBDF

"""
QBDF1: Multistep Method

An alias of `QNDF1` with κ=0.
"""
QBDF1(; kwargs...) = QNDF1(; kappa = 0, kwargs...)

"""
QBDF2: Multistep Method

An alias of `QNDF2` with κ=0.
"""
QBDF2(; kwargs...) = QNDF2(; kappa = 0, kwargs...)

"""
QBDF: Multistep Method

An alias of `QNDF` with κ=0.
"""
QBDF(; kwargs...) = QNDF(; kappa = tuple(0 // 1, 0 // 1, 0 // 1, 0 // 1, 0 // 1), kwargs...)

"""
    IMEXEuler(;kwargs...)

The one-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

When applied to a `SplitODEProblem` of the form

```
u'(t) = f1(u) + f2(u)
```

The default `IMEXEuler()` method uses an update of the form

```
unew = uold + dt * (f1(unew) + f2(uold))
```

See also `SBDF`, `IMEXEulerARK`.
"""
IMEXEuler(; kwargs...) = SBDF(1; kwargs...)

"""
    IMEXEulerARK(;kwargs...)

The one-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

When applied to a `SplitODEProblem` of the form

```
u'(t) = f1(u) + f2(u)
```

A classical additive Runge-Kutta method in the sense of
[Araújo, Murua, Sanz-Serna (1997)](https://doi.org/10.1137/S0036142995292128)
consisting of the implicit and the explicit Euler method given by

```
y1   = uold + dt * f1(y1)
unew = uold + dt * (f1(unew) + f2(y1))
```

See also `SBDF`, `IMEXEuler`.
"""
IMEXEulerARK(; kwargs...) = SBDF(1; ark = true, kwargs...)

@doc BDF_docstring(
    "1st order A-L and stiffly stable adaptive implicit Euler. Implicit Euler for implicit DAE form.
It uses an apriori error estimator for adaptivity based on a finite differencing approximation from SPICE.",
    "DImplicitEuler",
    extra_keyword_description = """
    - `nlsolve`: nonlinear solver algorithm used for solving the implicit system.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    """,
    extra_keyword_default = """
    nlsolve = NLNewton(),
    extrapolant = :constant,
    """
)
struct DImplicitEuler{AD, F, F2, CJ} <: DAEAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end
function DImplicitEuler(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
    )
    autodiff = _fixup_ad(autodiff)

    return DImplicitEuler(
        linsolve,
        nlsolve, extrapolant, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc BDF_docstring(
    "2nd order A-L stable adaptive BDF method. Fully implicit implementation of BDF2.",
    "DABDF2",
    references = """@article{celaya2014implementation,
    title={Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s},
    author={Celaya, E Alberdi and Aguirrezabala, JJ Anza and Chatzipantelidis, Panagiotis},
    journal={Procedia Computer Science},
    volume={29},
    pages={1014--1026},
    year={2014},
    publisher={Elsevier}}""",
    extra_keyword_description = """
    - `nlsolve`: nonlinear solver algorithm used for solving the implicit system.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    """,
    extra_keyword_default = """
    nlsolve = NLNewton(),
    extrapolant = :constant,
    """
)
struct DABDF2{AD, F, F2, CJ} <: DAEAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end
function DABDF2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
    )
    autodiff = _fixup_ad(autodiff)

    return DABDF2(
        linsolve,
        nlsolve, extrapolant, autodiff,
        _unwrap_val(concrete_jac)
    )
end

@doc BDF_docstring(
    "Fixed-leading coefficient adaptive-order adaptive-time BDF method. Fully implicit implementation of FBDF based on Shampine's",
    "DFBDF",
    references = """@article{shampine2002solving,
    title={Solving 0= F (t, y (t), y′(t)) in Matlab},
    author={Shampine, Lawrence F},
    year={2002},
    publisher={Walter de Gruyter GmbH and Co. KG}
    }""",
    extra_keyword_description = """
    - `κ`: coefficient for the order and stability control of the BDF method. When `nothing`, the default value is used.
    - `tol`: tolerance for the nonlinear solver. When `nothing`, uses the default tolerance.
    - `nlsolve`: nonlinear solver algorithm used for solving the implicit system.
    - `extrapolant`: extrapolation method used for the initial guess in the nonlinear solve.
    - `max_order`: maximum order of the adaptive-order BDF method.
    """,
    extra_keyword_default = """
    κ = nothing,
    tol = nothing,
    nlsolve = NLNewton(),
    extrapolant = :linear,
    max_order::Val{MO} = Val{5}(),
    """
)
struct DFBDF{MO, AD, F, F2, K, T, CJ, QT} <: DAEAlgorithm
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    κ::K
    tol::T
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
    qmax::QT
    qsteady_min::QT
    qsteady_max::QT
end
function DFBDF(;
        max_order::Val{MO} = Val{5}(),
        autodiff = AutoForwardDiff(), concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear,
        qsteady_min = 1 // 1, qsteady_max = 2 // 1, qmax = 10 // 1,
    ) where {MO}
    autodiff = _fixup_ad(autodiff)

    return DFBDF(
        max_order, linsolve, nlsolve, κ, tol, extrapolant,
        autodiff,
        _unwrap_val(concrete_jac),
        qmax, qsteady_min, qsteady_max,
    )
end

@truncate_stacktrace DFBDF
