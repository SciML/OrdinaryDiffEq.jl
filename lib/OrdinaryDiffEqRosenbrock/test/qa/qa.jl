using SciMLTesting, OrdinaryDiffEqRosenbrock, Test

# The names ignored below are internal (non-public) cross-package interface
# names. OrdinaryDiffEqRosenbrock is a solver sublibrary that, by design,
# consumes the internal solver interfaces of its sibling packages
# OrdinaryDiffEqCore and OrdinaryDiffEqDifferentiation (and a few DiffEqBase /
# SciMLBase / Base internals). These names have no public home to migrate to,
# so they are listed explicitly rather than blanket-skipped via `ei_broken`, so
# that any name later made public will surface as an Unexpected Pass and can be
# dropped from these lists. See SciML/OrdinaryDiffEq.jl#3776.

# Internal names imported via a non-owner sibling package (re-exported through
# OrdinaryDiffEqCore / OrdinaryDiffEqDifferentiation) or otherwise not public in
# the module they are imported from.
const ROSENBROCK_INTERNAL_EXPLICIT_IMPORTS = (
    Symbol("@cache"), Symbol("@def"),
    :DerivativeOrderNotPossibleError, :DifferentialVarsUndefined, :ODEIntegrator,
    :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqMutableCache,
    :OrdinaryDiffEqRosenbrockAdaptiveAlgorithm, :OrdinaryDiffEqRosenbrockAlgorithm,
    :TimeDerivativeWrapper, :TimeGradientWrapper, :UDerivativeWrapper,
    :UJacobianWrapper, :WOperator,
    :_fixup_ad, :_ode_addsteps!, :_ode_interpolant, :_ode_interpolant!,
    :_reshape, :_unwrap_val, :_vec,
    :alg_adaptive_order, :alg_autodiff, :alg_cache,
    :build_J_W, :build_grad_config, :build_jac_config,
    :calc_rosenbrock_differentiation, :calc_rosenbrock_differentiation!,
    :calc_tderivative, :calculate_residuals, :calculate_residuals!,
    :constvalue, :copyat_or_push!, :dolinsolve, :generic_solver_docstring,
    :get_fsalfirstlast, :has_stiff_interpolation, :initialize!, :isWmethod,
    :isfsal, :issuccess_W, :jacobian2W!, :only_diagonal_mass_matrix,
    :perform_step!, :resize_J_W!, :resize_grad_config!, :resize_jac_config!,
    :trivial_limiter!, :wrapprecs,
)

# Internal names accessed via qualified access on a non-owner / non-public name.
const ROSENBROCK_INTERNAL_QUALIFIED_ACCESSES = (
    Symbol("@set"),               # SciMLBase.@set (owned by Accessors, non-public re-export)
    :FunctionWrapperSpecialize,   # SciMLBase internal (precompile workload only)
    :NoSpecialize,                # SciMLBase internal (precompile workload only)
    :get_EEst, :set_EEst!, :increment_nf!,  # OrdinaryDiffEqCore internal
    :lorenz, :lorenz_oop,         # OrdinaryDiffEqCore precompile test helpers
    :setindex,                    # Base.setindex (non-public)
)

run_qa(
    OrdinaryDiffEqRosenbrock;
    explicit_imports = true,
    ei_kwargs = (
        all_explicit_imports_via_owners = (; ignore = ROSENBROCK_INTERNAL_EXPLICIT_IMPORTS),
        all_explicit_imports_are_public = (; ignore = ROSENBROCK_INTERNAL_EXPLICIT_IMPORTS),
        all_qualified_accesses_via_owners = (; ignore = ROSENBROCK_INTERNAL_QUALIFIED_ACCESSES),
        all_qualified_accesses_are_public = (; ignore = ROSENBROCK_INTERNAL_QUALIFIED_ACCESSES),
    ),
)
