module OrdinaryDiffEqRosenbrock

import OrdinaryDiffEq: alg_order, alg_adaptive_order, isWmethod, isfsal, _unwrap_val,
                       DEFAULT_PRECS, OrdinaryDiffEqRosenbrockAlgorithm, @cache, alg_cache, initialize!, @unpack,
                       calc_W, calculate_residuals!, calc_rosenbrock_differentiation!, OrdinaryDiffEqMutableCache,
                       build_J_W, UJacobianWrapper, OrdinaryDiffEqConstantCache, _ode_interpolant, _ode_interpolant!,
                       _vec, _reshape, perform_step!, trivial_limiter!, dolinsolve, OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
                       OrdinaryDiffEqRosenbrockAlgorithm, generic_solver_docstring, namify, initialize!, perform_step!,
                       constvalue, TimeDerivativeWrapper, TimeGradientWrapper, UDerivativeWrapper, UJacobianWrapper,
                       wrapprecs, alg_autodiff, calc_tderivative, build_grad_config, build_jac_config,
                       issuccess_W, calculate_residuals, has_stiff_interpolation,
                       resize_non_user_cache!, _ode_addsteps!
using TruncatedStacktraces, MuladdMacro, FastBroadcast, DiffEqBase, RecursiveArrayTools
import MacroTools
using MacroTools: @capture
using DiffEqBase: @def
import LinearSolve
import ForwardDiff
using LinearAlgebra: mul!, diag, diagm, I, Diagonal

function rosenbrock_wanner_docstring(description::String,
        name::String;
        references::String = "",
        extra_keyword_description = "",
        extra_keyword_default = "",
        with_step_limiter = false)
    keyword_default = """
    autodiff = Val{true}(),
    concrete_jac = nothing,
    linsolve = nothing,
    precs = DEFAULT_PRECS,
    """ * extra_keyword_default

    keyword_default_description = """
    - `autodiff`: boolean to control if the Jacobian should be computed via AD or not
    - `concrete_jac`: function of the form `jac!(J, u, p, t)`
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

function rosenbrock_docstring(description::String,
        name::String;
        references::String = "",
        extra_keyword_description = "",
        extra_keyword_default = "",
        with_step_limiter = false)
    keyword_default = """
    autodiff = Val{true}(),
    concrete_jac = nothing,
    linsolve = nothing,
    precs = DEFAULT_PRECS,
    """ * extra_keyword_default

    keyword_default_description = """
    - `autodiff`: boolean to control if the Jacobian should be computed via AD or not
    - `concrete_jac`: function of the form `jac!(J, u, p, t)`
    - `linsolve`: custom solver for the inner linear systems
    - `precs`: custom preconditioner for the inner linear solver
    """ * extra_keyword_description

    if with_step_limiter
        keyword_default *= "step_limiter! = OrdinaryDiffEq.trivial_limiter!,\n"
        keyword_default_description *= "- `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`\n"
    end

    generic_solver_docstring(
        description, name, "Rosenbrock Method. ", references,
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

export Rosenbrock23, Rosenbrock32, RosShamp4, Veldd4, Velds4, GRK4T, GRK4A,
       Ros4LStab, ROS3P, Rodas3, Rodas23W, Rodas3P, Rodas4, Rodas42, Rodas4P, Rodas4P2,
       Rodas5, Rodas5P, Rodas5Pe, Rodas5Pr,
       RosenbrockW6S4OS, ROS34PW1a, ROS34PW1b, ROS34PW2, ROS34PW3, ROS34PRw, ROS3PRL,
       ROS3PRL2, ROK4a,
       ROS2, ROS2PR, ROS2S, ROS3, ROS3PR, Scholz4_7

end
