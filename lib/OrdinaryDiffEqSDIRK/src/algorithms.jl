function SDIRK_docstring(description::String,
        name::String;
        references::String = "",
        extra_keyword_description::String = "",
        extra_keyword_default::String = "")
    keyword_default = """
        chunk_size = Val{0}(),
        autodiff = AutoForwardDiff(),
        standardtag = Val{true}(),
        concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing,
        precs = DEFAULT_PRECS,
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
    references = "@book{wanner1996solving,
    title={Solving ordinary differential equations II},
    author={Wanner, Gerhard and Hairer, Ernst},
    volume={375},
    year={1996},
    publisher={Springer Berlin Heidelberg New York}}",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
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
    autodiff::AD
end

function ImplicitEuler(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    ImplicitEuler{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve, precs, extrapolant, controller, step_limiter!, AD_choice)
end

@doc SDIRK_docstring("A second order A-stable symplectic and symmetric implicit solver.
    Good for highly stiff equations which need symplectic integration.",
    "ImplicitMidpoint";
    references = "@book{wanner1996solving,
    title={Solving ordinary differential equations II},
    author={Wanner, Gerhard and Hairer, Ernst},
    volume={375},
    year={1996},
    publisher={Springer Berlin Heidelberg New York}}",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
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
    autodiff::AD
end

function ImplicitMidpoint(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    ImplicitMidpoint{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        step_limiter!, AD_choice)
end

@doc SDIRK_docstring(
    """Second order A-stable symmetric ESDIRK method.
"Almost symplectic" without numerical dampening.
Also known as Crank-Nicolson when applied to PDEs. Adaptive timestepping via divided
differences approximation to the second derivative terms in the local truncation error
estimate (the SPICE approximation strategy).""",
    "Trapezoid";
    references = "Andre Vladimirescu. 1994. The Spice Book. John Wiley & Sons, Inc., New York, NY, USA.",
    extra_keyword_description = """
    - `extrapolant`: TBD
    - `controller`: TBD
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
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
    autodiff::AD
end

function Trapezoid(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    Trapezoid{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        controller,
        step_limiter!,
        AD_choice)
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
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
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
    autodiff::AD
end

function TRBDF2(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    TRBDF2{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!, AD_choice)
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
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
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
    autodiff::AD
end

function SDIRK2(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    SDIRK2{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(
        linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller,
        step_limiter!,
        AD_choice)
end

@doc SDIRK_docstring("Description TBD",
    "SDIRK22";
    references = "@techreport{kennedy2016diagonally,
    title={Diagonally implicit Runge-Kutta methods for ordinary differential equations. A review},
    author={Kennedy, Christopher A and Carpenter, Mark H},
    year={2016}}",
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
    """)
struct SDIRK22{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
end

function SDIRK22(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    Trapezoid{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        controller,
        step_limiter!,
        AD_choice)
end

@doc SDIRK_docstring(
    """SSPSDIRK is an SSP-optimized SDIRK method,
so it's an implicit SDIRK method for handling stiffness but if the `dt` is below the SSP `coefficient * dt`,
then the SSP property of the SSP integrators (the other page) is satisified.
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
    autodiff::AD
end

function SSPSDIRK2(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :constant,
        controller = :PI)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    SSPSDIRK2{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller, AD_choice)
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
    - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
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
    autodiff::AD
end
function Kvaerno3(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    Kvaerno3{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!, AD_choice)
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
    autodiff::AD
end
function KenCarp3(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    KenCarp3{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!, AD_choice)
end

@doc SDIRK_docstring("Third order method.",
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
    autodiff::AD
end
function CFNLIRK3(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    CFNLIRK3{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        AD_choice)
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
    autodiff::AD
end
function Cash4(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, embedding = 3)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    Cash4{
        _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(
        linsolve,
        nlsolve,
        precs,
        smooth_est,
        extrapolant,
        embedding,
        controller,
        AD_choice)
end

@doc SDIRK_docstring("Method of order 4.",
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
    autodiff::AD
end
function SFSDIRK4(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    SFSDIRK4{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        AD_choice)
end

@doc SDIRK_docstring("Method of order 5.",
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
    autodiff::AD
end

function SFSDIRK5(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    SFSDIRK5{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        AD_choice)
end

@doc SDIRK_docstring("Method of order 6.",
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
    autodiff::AD
end

function SFSDIRK6(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    SFSDIRK6{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        AD_choice)
end

@doc SDIRK_docstring("Method of order 7.",
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
    autodiff::AD
end

function SFSDIRK7(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    SFSDIRK7{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        AD_choice)
end

@doc SDIRK_docstring("Method of order 8.",
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
    autodiff::AD
end

function SFSDIRK8(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    SFSDIRK8{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        AD_choice)
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
    autodiff::AD
end
function Hairer4(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    Hairer4{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller, AD_choice)
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
    autodiff::AD
end
function Hairer42(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    Hairer42{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller, AD_choice)
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
    autodiff::AD
end
function Kvaerno4(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    Kvaerno4{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!, AD_choice)
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
    autodiff::AD
end
function Kvaerno5(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    Kvaerno5{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!, AD_choice)
end

@doc SDIRK_docstring(
    "An A-L stable stiffly-accurate 4th order ESDIRK method with splitting.",
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
    autodiff::AD
end
function KenCarp4(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    KenCarp4{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!, AD_choice)
end

TruncatedStacktraces.@truncate_stacktrace KenCarp4

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
    autodiff::AD
end
function KenCarp47(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    KenCarp47{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller, AD_choice)
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
    autodiff::AD
end
function KenCarp5(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    KenCarp5{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!, AD_choice)
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
    autodiff::AD
end
function KenCarp58(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    KenCarp58{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller, AD_choice)
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
    autodiff::AD
end
function ESDIRK54I8L2SA(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    ESDIRK54I8L2SA{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller, AD_choice)
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
    autodiff::AD
end
function ESDIRK436L2SA2(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    ESDIRK436L2SA2{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller, AD_choice)
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
    autodiff::AD
end
function ESDIRK437L2SA(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    ESDIRK437L2SA{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller, AD_choice)
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
    autodiff::AD
end
function ESDIRK547L2SA2(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    ESDIRK547L2SA2{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller, AD_choice)
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
    autodiff::AD
end
function ESDIRK659L2SA(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    ESDIRK659L2SA{_unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller, AD_choice)
end
