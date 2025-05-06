function BDF_docstring(description::String,
        name::String;
        references::String = "",
        extra_keyword_description::String = "",
        extra_keyword_default::String = "")
    keyword_default = """
        chunk_size = Val{0}(),
        autodiff = AutoForwardDiff(),
        standardtag = Val{true}(),
        concrete_jac = nothing,
        linsolve = nothing,
        precs = DEFAULT_PRECS,
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
    - `standardtag`: Specifies whether to use package-specific tags instead of the
        ForwardDiff default function-specific tags. For more information, see
        [this blog post](https://www.stochasticlifestyle.com/improved-forwarddiff-jl-stacktraces-with-package-tags/).
        Defaults to `Val{true}()`.
    - `concrete_jac`: Specifies whether a Jacobian should be constructed. Defaults to
        `nothing`, which means it will be chosen true/false depending on circumstances
        of the solver, such as whether a Krylov subspace method is used for `linsolve`.
    - `linsolve`: Any [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) compatible linear solver.
      For example, to use [KLU.jl](https://github.com/JuliaSparse/KLU.jl), specify
      `$name(linsolve = KLUFactorization()`).
       When `nothing` is passed, uses `DefaultLinearSolver`.
    - `precs`: Any [LinearSolve.jl-compatible preconditioner](https://docs.sciml.ai/LinearSolve/stable/basics/Preconditioners/)
      can be used as a left or right preconditioner.
      Preconditioners are specified by the `Pl,Pr = precs(W,du,u,p,t,newW,Plprev,Prprev,solverdata)`
      function where the arguments are defined as:
        - `W`: the current Jacobian of the nonlinear system. Specified as either
            ``I - \\gamma J`` or ``I/\\gamma - J`` depending on the algorithm. This will
            commonly be a `WOperator` type defined by OrdinaryDiffEq.jl. It is a lazy
            representation of the operator. Users can construct the W-matrix on demand
            by calling `convert(AbstractMatrix,W)` to receive an `AbstractMatrix` matching
            the `jac_prototype`.
        - `du`: the current ODE derivative
        - `u`: the current ODE state
        - `p`: the ODE parameters
        - `t`: the current ODE time
        - `newW`: a `Bool` which specifies whether the `W` matrix has been updated since
            the last call to `precs`. It is recommended that this is checked to only
            update the preconditioner when `newW == true`.
        - `Plprev`: the previous `Pl`.
        - `Prprev`: the previous `Pr`.
        - `solverdata`: Optional extra data the solvers can give to the `precs` function.
            Solver-dependent and subject to change.
      The return is a tuple `(Pl,Pr)` of the LinearSolve.jl-compatible preconditioners.
      To specify one-sided preconditioning, simply return `nothing` for the preconditioner
      which is not used. Additionally, `precs` must supply the dispatch:
      ```julia
      Pl, Pr = precs(W, du, u, p, t, ::Nothing, ::Nothing, ::Nothing, solverdata)
      ```
      which is used in the solver setup phase to construct the integrator
      type with the preconditioners `(Pl,Pr)`.
      The default is `precs=DEFAULT_PRECS` where the default preconditioner function
      is defined as:
      ```julia
      DEFAULT_PRECS(W, du, u, p, t, newW, Plprev, Prprev, solverdata) = nothing, nothing
      ```
        """ * "/n" * extra_keyword_description
    generic_solver_docstring(
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
    - `κ`: TBD
    - `tol`: TBD
    - `nlsolve`: TBD
    - `smooth_est`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    κ = nothing,
    tol = nothing,
    nlsolve = NLNewton(),
    smooth_est = true,
    extrapolant = :linear,
    controller = :Standard,
    step_limiter! = trivial_limiter!,
    """)
struct ABDF2{CS, AD, F, F2, P, FDT, ST, CJ, K, T, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
end
function ABDF2(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        κ = nothing, tol = nothing, linsolve = nothing, precs = DEFAULT_PRECS,
        nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :Standard, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    ABDF2{
        _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol), typeof(step_limiter!)}(linsolve, nlsolve, precs, κ, tol,
        smooth_est, extrapolant, controller, step_limiter!, AD_choice)
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
    - `κ`: TBD
    - `tol`: TBD
    - `nlsolve`: TBD
    - `extrapolant`: TBD
    - `ark`: TBD
    """,
    extra_keyword_default = """
    κ = nothing,
    tol = nothing,
    nlsolve = NLNewton(),
    extrapolant = :linear,
    ark = false,
    order,
    """)
struct SBDF{CS, AD, F, F2, P, FDT, ST, CJ, K, T} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    order::Int
    ark::Bool
    autodiff::AD
end

function SBDF(order; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, ark = false)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    SBDF{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol)}(linsolve,
        nlsolve,
        precs,
        κ,
        tol,
        extrapolant,
        order,
        ark,
        AD_choice)
end

# All keyword form needed for remake
function SBDF(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear,
        order, ark = false)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    SBDF{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol)}(linsolve,
        nlsolve,
        precs,
        κ,
        tol,
        extrapolant,
        order,
        ark,
        AD_choice)
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
    "An adaptive order 1 quasi-constant timestep L-stable numerical differentiation function (NDF) method.
Optional parameter kappa defaults to Shampine's accuracy-optimal -0.1850.",
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
    - `nlsolve`: TBD
    - `extrapolant`: TBD
    - `kappa`: TBD
    - `controller`: TBD
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    nlsolve = NLNewton(),
    extrapolant = :linear,
    kappa = -0.1850,
    controller = :Standard,
    step_limiter! = trivial_limiter!,
    """)
struct QNDF1{CS, AD, F, F2, P, FDT, ST, CJ, κType, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    kappa::κType
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
end

function QNDF1(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, kappa = -37 // 200,
        controller = :Standard, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    QNDF1{
        _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(kappa), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        kappa,
        controller,
        step_limiter!,
        AD_choice)
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
    - `nlsolve`: TBD
    - `extrapolant`: TBD
    - `kappa`: TBD
    - `controller`: TBD
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    nlsolve = NLNewton(),
    extrapolant = :linear,
    kappa =  -1 // 9,
    controller = :Standard,
    step_limiter! = trivial_limiter!,
    """)
struct QNDF2{CS, AD, F, F2, P, FDT, ST, CJ, κType, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    kappa::κType
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
end

function QNDF2(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, kappa = -1 // 9,
        controller = :Standard, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    QNDF2{
        _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(kappa), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        kappa,
        controller,
        step_limiter!,
        AD_choice)
end

@doc BDF_docstring(
    "An adaptive order quasi-constant timestep NDF method.
Utilizes Shampine's accuracy-optimal kappa values as defaults (has a keyword argument for a tuple of kappa coefficients).",
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
    - `κ`: TBD
    - `tol`: TBD
    - `nlsolve`: TBD
    - `extrapolant`: TBD
    - `kappa`: TBD
    - `controller`: TBD
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    """,
    extra_keyword_default = """
    κ = nothing,
    tol = nothing,
    nlsolve = NLNewton(),
    extrapolant = :linear,
    kappa =  promote(-0.1850, -1 // 9, -0.0823, -0.0415, 0),
    controller = :Standard,
    step_limiter! = trivial_limiter!,
    """)
struct QNDF{MO, CS, AD, F, F2, P, FDT, ST, CJ, K, T, κType, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    kappa::κType
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
end

function QNDF(; max_order::Val{MO} = Val{5}(), chunk_size = Val{0}(),
        autodiff = AutoForwardDiff(), standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, kappa = (
            -37 // 200, -1 // 9, -823 // 10000, -83 // 2000, 0 // 1),
        controller = :Standard, step_limiter! = trivial_limiter!) where {MO}
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    QNDF{MO, _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol), typeof(kappa), typeof(step_limiter!)}(
        max_order, linsolve, nlsolve, precs, κ, tol,
        extrapolant, kappa, controller, step_limiter!, AD_choice)
end

TruncatedStacktraces.@truncate_stacktrace QNDF

@doc BDF_docstring("The second order Modified Extended BDF method,
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
    - `nlsolve`: TBD
    - `extrapolant`: TBD
    """,
    extra_keyword_default = """
    nlsolve = NLNewton(),
    extrapolant = :constant,
    """)
struct MEBDF2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    autodiff::AD
end
function MEBDF2(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    MEBDF2{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        AD_choice)
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
    - `κ`: TBD
    - `tol`: TBD
    - `nlsolve`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
    - `max_order`: TBD
    """,
    extra_keyword_default = """
    κ = nothing,
    tol = nothing,
    nlsolve = NLNewton(),
    extrapolant = :linear,
    controller = :Standard,
    step_limiter! = trivial_limiter!,
    max_order::Val{MO} = Val{5}(),
    """)
struct FBDF{MO, CS, AD, F, F2, P, FDT, ST, CJ, K, T, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
end

function FBDF(; max_order::Val{MO} = Val{5}(), chunk_size = Val{0}(),
        autodiff = AutoForwardDiff(), standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, controller = :Standard, step_limiter! = trivial_limiter!) where {MO}
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    FBDF{MO, _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol), typeof(step_limiter!)}(
        max_order, linsolve, nlsolve, precs, κ, tol, extrapolant,
        controller, step_limiter!, AD_choice)
end

TruncatedStacktraces.@truncate_stacktrace FBDF

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
    "Implicit Euler for implicit DAE form.
It uses an apriori error estimator for adaptivity based on a finite differencing approximation from SPICE.",
    "DImplicitEuler",
    extra_keyword_description = """
    - `nlsolve`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    nlsolve = NLNewton(),
    extrapolant = :constant,
    controller = :Standard,
    """)
struct DImplicitEuler{CS, AD, F, F2, P, FDT, ST, CJ} <: DAEAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
    autodiff::AD
end
function DImplicitEuler(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant,
        controller = :Standard)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    DImplicitEuler{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve, precs, extrapolant, controller, AD_choice)
end

@doc BDF_docstring("Fully implicit implementation of BDF2.",
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
    - `nlsolve`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    """,
    extra_keyword_default = """
    nlsolve = NLNewton(),
    extrapolant = :constant,
    controller = :Standard,
    """)
struct DABDF2{CS, AD, F, F2, P, FDT, ST, CJ} <: DAEAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
    autodiff::AD
end
function DABDF2(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant,
        controller = :Standard)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    DABDF2{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve, precs, extrapolant, controller, AD_choice)
end

#=
struct DBDF{CS,AD,F,F2,P,FDT,ST,CJ} <: DAEAlgorithm{CS,AD,FDT,ST,CJ}
  linsolve::F
  nlsolve::F2
  precs::P
  extrapolant::Symbol
end

DBDF(;chunk_size=Val{0}(),autodiff=Val{true}(), standardtag = Val{true}(), concrete_jac = nothing,diff_type=Val{:forward},
     linsolve=nothing,precs = DEFAULT_PRECS,nlsolve=NLNewton(),extrapolant=:linear) =
     DBDF{_unwrap_val(chunk_size),_unwrap_val(autodiff),typeof(linsolve),typeof(nlsolve),typeof(precs),diff_type,_unwrap_val(standardtag),_unwrap_val(concrete_jac)}(
     linsolve,nlsolve,precs,extrapolant)
=#

@doc BDF_docstring("Fully implicit implementation of FBDF based on Shampine's",
    "DFBDF",
    references = """@article{shampine2002solving,
    title={Solving 0= F (t, y (t), y′(t)) in Matlab},
    author={Shampine, Lawrence F},
    year={2002},
    publisher={Walter de Gruyter GmbH and Co. KG}
    }""",
    extra_keyword_description = """
    - `κ`: TBD
    - `tol`: TBD
    - `nlsolve`: TBD
    - `extrapolant`: TBD
    - `controller`: TBD
    - `max_order`: TBD
    """,
    extra_keyword_default = """
    κ = nothing,
    tol = nothing,
    nlsolve = NLNewton(),
    extrapolant = :linear,
    controller = :Standard,
    max_order::Val{MO} = Val{5}(),
    """)
struct DFBDF{MO, CS, AD, F, F2, P, FDT, ST, CJ, K, T} <: DAEAlgorithm{CS, AD, FDT, ST, CJ}
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    controller::Symbol
    autodiff::AD
end
function DFBDF(; max_order::Val{MO} = Val{5}(), chunk_size = Val{0}(),
        autodiff = AutoForwardDiff(), standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, controller = :Standard) where {MO}
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    DFBDF{MO, _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol)}(max_order, linsolve, nlsolve, precs, κ, tol, extrapolant,
        controller, AD_choice)
end

TruncatedStacktraces.@truncate_stacktrace DFBDF
