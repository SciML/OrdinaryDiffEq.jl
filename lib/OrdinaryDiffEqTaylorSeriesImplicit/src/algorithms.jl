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
    "ImplicitTaylor1";
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
struct ImplicitTaylor1{T, CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    μ::T
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
    autodiff::AD
end

function ImplicitTaylor1(; μ = 1.0, chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(),
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant,
        controller = :PI, step_limiter! = trivial_limiter!)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

    ImplicitTaylor1{typeof(μ), _unwrap_val(chunk_size), typeof(AD_choice), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(μ, linsolve,
        nlsolve, precs, extrapolant, controller, step_limiter!, AD_choice)
end
