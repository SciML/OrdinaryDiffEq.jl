module OrdinaryDiffEqRosenbrock

import OrdinaryDiffEqCore: alg_order, alg_adaptive_order, isWmethod, isfsal, _unwrap_val,
                           DEFAULT_PRECS, OrdinaryDiffEqRosenbrockAlgorithm, @cache,
                           alg_cache, initialize!, @unpack,
                           calculate_residuals!, OrdinaryDiffEqMutableCache,
                           OrdinaryDiffEqConstantCache, _ode_interpolant, _ode_interpolant!,
                           _vec, _reshape, perform_step!, trivial_limiter!,
                           OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
                           OrdinaryDiffEqRosenbrockAlgorithm, generic_solver_docstring,
                           namify, initialize!, perform_step!, get_fsalfirstlast,
                           constvalue, only_diagonal_mass_matrix,
                           calculate_residuals, has_stiff_interpolation, ODEIntegrator,
                           resize_non_user_cache!, _ode_addsteps!, full_cache,
                           DerivativeOrderNotPossibleError, _bool_to_ADType,
                           _process_AD_choice, LinearAliasSpecifier
using MuladdMacro, FastBroadcast, RecursiveArrayTools
import MacroTools
using MacroTools: @capture
using DiffEqBase: @def
import LinearSolve
import LinearSolve: UniformScaling
import ForwardDiff
using FiniteDiff
using LinearAlgebra: mul!, diag, diagm, I, Diagonal, norm
import ADTypes: AutoForwardDiff, AbstractADType
import OrdinaryDiffEqCore

using OrdinaryDiffEqDifferentiation: TimeDerivativeWrapper, TimeGradientWrapper,
                                     UDerivativeWrapper, UJacobianWrapper,
                                     wrapprecs, calc_tderivative, build_grad_config,
                                     build_jac_config, issuccess_W, jacobian2W!,
                                     resize_jac_config!, resize_grad_config!,
                                     calc_W, calc_rosenbrock_differentiation!, build_J_W,
                                     UJacobianWrapper, dolinsolve

using Reexport
@reexport using DiffEqBase

import OrdinaryDiffEqCore: alg_autodiff
import OrdinaryDiffEqCore

function rosenbrock_wolfbrandt_docstring(description::String,
        name::String;
        references::String = "",
        extra_keyword_description = "",
        extra_keyword_default = "",
        with_step_limiter = false)
    keyword_default = """
    chunk_size = Val{0}(),
    standardtag = Val{true}(),
    autodiff = AutoForwardDiff(),
    concrete_jac = nothing,
    diff_type = Val{:forward}(),
    linsolve = nothing,
    precs = DEFAULT_PRECS,
    """ * extra_keyword_default

    keyword_default_description = """
    - `standardtag`: Specifies whether to use package-specific tags instead of the
        ForwardDiff default function-specific tags. For more information, see
        [this blog post](https://www.stochasticlifestyle.com/improved-forwarddiff-jl-stacktraces-with-package-tags/).
        Defaults to `Val{true}()`.
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
    """ * extra_keyword_description

    if with_step_limiter
        keyword_default *= "step_limiter! = OrdinaryDiffEq.trivial_limiter!,\n"
        keyword_default_description *= "- `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`\n"
    end

    generic_solver_docstring(
        description, name, "Rosenbrock-Wanner-W(olfbrandt) Method. ", references,
        keyword_default_description, keyword_default
    )
end

function rosenbrock_docstring(description::String,
        name::String;
        references::String = "",
        extra_keyword_description = "",
        extra_keyword_default = "",
        with_step_limiter = false)
    keyword_default = """
    - `standardtag`: Specifies whether to use package-specific tags instead of the
        ForwardDiff default function-specific tags. For more information, see
        [this blog post](https://www.stochasticlifestyle.com/improved-forwarddiff-jl-stacktraces-with-package-tags/).
        Defaults to `Val{true}()`.
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
    """ * extra_keyword_default

    keyword_default_description = """
    - `chunk_size`: TBD
    - `standardtag`: TBD
    - `autodiff`: boolean to control if the Jacobian should be computed via AD or not
    - `concrete_jac`: function of the form `jac!(J, u, p, t)`
    - `diff_type`: TBD
    - `linsolve`: custom solver for the inner linear systems
    - `precs`: custom preconditioner for the inner linear solver
    """ * extra_keyword_description

    if with_step_limiter
        keyword_default *= "step_limiter! = OrdinaryDiffEq.trivial_limiter!,\n"
        keyword_default_description *= "- `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`\n"
    end

    generic_solver_docstring(
        description, name, "Rosenbrock-Wanner Method. ", references,
        keyword_default_description, keyword_default
    )
end

include("algorithms.jl")
include("alg_utils.jl")
include("generic_rosenbrock.jl")
include("rosenbrock_caches.jl")
include("rosenbrock_tableaus.jl")
include("interp_func.jl")
include("rosenbrock_interpolants.jl")
include("stiff_addsteps.jl")
include("rosenbrock_perform_step.jl")
include("integrator_interface.jl")

import PrecompileTools
import Preferences
PrecompileTools.@compile_workload begin
    lorenz = OrdinaryDiffEqCore.lorenz
    lorenz_oop = OrdinaryDiffEqCore.lorenz_oop
    solver_list = [Rosenbrock23(), Rodas5P()]
    prob_list = []

    if Preferences.@load_preference("PrecompileDefaultSpecialize", true)
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0)))
        push!(prob_list, ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileAutoSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.AutoSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileFunctionWrapperSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(lorenz, [1.0; 0.0; 0.0],
                (0.0, 1.0), Float64[]))
    end

    if Preferences.@load_preference("PrecompileNoSpecialize", false)
        push!(prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0)))
        push!(prob_list,
            ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0),
                Float64[]))
    end

    for prob in prob_list, solver in solver_list
        solve(prob, solver)(5.0)
    end

    prob_list = nothing
    solver_list = nothing
end

export Rosenbrock23, Rosenbrock32, RosShamp4, Veldd4, Velds4, GRK4T, GRK4A,
       Ros4LStab, ROS3P, Rodas3, Rodas23W, Rodas3P, Rodas4, Rodas42, Rodas4P, Rodas4P2,
       Rodas5, Rodas5P, Rodas5Pe, Rodas5Pr,
       RosenbrockW6S4OS, ROS34PW1a, ROS34PW1b, ROS34PW2, ROS34PW3, ROS34PRw, ROS3PRL,
       ROS3PRL2, ROK4a,
       ROS2, ROS2PR, ROS2S, ROS3, ROS3PR, Scholz4_7

end
