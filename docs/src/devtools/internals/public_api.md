# [Developer Extension API](@id Developer-Extension-API)

This page documents the version-controlled API intended for solver authors and
OrdinaryDiffEq monorepo subpackages. These names are not application-facing
OrdinaryDiffEq user API; user-facing DiffEqBase and OrdinaryDiffEqCore APIs are
documented under the API section.

!!! warning "Developer API, not user API"

    Do not build application code against these hooks. They exist so solver
    packages can extend common traits, controllers, interpolation hooks, and
    initialization protocols without depending on undocumented implementation
    details. Concrete per-algorithm caches and low-level nonlinear-solve helpers
    remain internal unless listed below; the listed automatic-switch caches and
    nonlinear-solver hooks are shared with sibling integrator packages.

## DiffEqBase solver hooks

```@docs
DiffEqBase.CallbackCache
DiffEqBase.EvalFunc
DiffEqBase.OrdinaryDiffEqTag
DiffEqBase.apply_callback!
DiffEqBase.apply_discrete_callback!
DiffEqBase.calculate_residuals
DiffEqBase.calculate_residuals!
DiffEqBase.check_prob_alg_pairing
DiffEqBase.default_factorize
DiffEqBase.finalize!
DiffEqBase.find_callback_time
DiffEqBase.find_first_continuous_callback
DiffEqBase.get_condition
DiffEqBase.get_tstops
DiffEqBase.get_tstops_array
DiffEqBase.get_tstops_max
DiffEqBase.initialize!
DiffEqBase.max_vector_callback_length
DiffEqBase.max_vector_callback_length_int
DiffEqBase.merge_problem_kwargs
DiffEqBase.prepare_alg
DiffEqBase.prob2dtmin
DiffEqBase.stripunits
DiffEqBase.timedepentdtmin
```

## Solver code-generation utilities

These macros are versioned for solver packages, not application code. Generated
tableau bindings are local implementation details. Foldability and loop-policy
annotations may only be applied when their documented effect assumptions hold.

```@docs
DiffEqBase.@tight_loop_macros
OrdinaryDiffEqCore.@OnDemandTableauExtract
OrdinaryDiffEqCore.@fold
```

## Algorithm type hierarchy

Every solver subtypes one of these abstract algorithm types.

### Minimal algorithm contract

Solver packages should subtype the narrowest applicable algorithm type, define
the order, and then implement the cache and stepping hooks. The generic trait
defaults are intentional: an explicit fixed-step method needs only the methods
below before adding its cache and stepping implementation.

```julia
using OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm

struct MyEuler <: OrdinaryDiffEqAlgorithm end
OrdinaryDiffEqCore.alg_order(::MyEuler) = 1
OrdinaryDiffEqCore.isfsal(::MyEuler) = false

@assert !OrdinaryDiffEqCore.isadaptive(MyEuler())
@assert !OrdinaryDiffEqCore.isimplicit(MyEuler())
```

The complete solver implementation then supplies:

1. `alg_cache(alg, ...)`, returning an `OrdinaryDiffEqConstantCache` or
   `OrdinaryDiffEqMutableCache`.
2. `initialize!(integrator, cache)`, initializing FSAL and cache state.
3. `perform_step!(integrator, cache, repeat_step = false)`, writing the candidate
   state to `integrator.u` and the error estimate to `integrator.EEst`.

Adaptive methods subtype `OrdinaryDiffEqAdaptiveAlgorithm` and must provide an
embedded error estimate suitable for the controller. Implicit methods subtype
`OrdinaryDiffEqImplicitAlgorithm` and must provide the nonlinear-solver cache
and stage operations expected by the selected `nlsolve` configuration.

```@docs
OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqAdaptiveAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqCompositeAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqImplicitAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqAdaptiveImplicitAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqNewtonAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqNewtonAdaptiveAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqRosenbrockAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqRosenbrockAdaptiveAlgorithm
OrdinaryDiffEqCore.NewtonAlgorithm
OrdinaryDiffEqCore.RosenbrockAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqExponentialAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqAdaptiveExponentialAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqLinearExponentialAlgorithm
OrdinaryDiffEqCore.ExponentialAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm
OrdinaryDiffEqCore.DAEAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqPartitionedAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqAdaptivePartitionedAlgorithm
OrdinaryDiffEqCore.PartitionedAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqImplicitSecondOrderAlgorithm
OrdinaryDiffEqCore.OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm
OrdinaryDiffEqCore.ImplicitSecondOrderAlgorithm
```

### SDE / RODE algorithm hierarchy

Subtyped by StochasticDiffEq.jl and defined in the core so shared machinery can
dispatch on these algorithms.

```@docs
OrdinaryDiffEqCore.StochasticDiffEqAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqAdaptiveAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqCompositeAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqNewtonAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqNewtonAdaptiveAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqRODEAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqRODEAdaptiveAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqRODECompositeAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqJumpAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqJumpAdaptiveAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqJumpNewtonAdaptiveAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqJumpDiffusionAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqJumpDiffusionAdaptiveAlgorithm
OrdinaryDiffEqCore.StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm
```

## Composite algorithms and automatic switching

```@docs
OrdinaryDiffEqCore.CompositeAlgorithm
OrdinaryDiffEqCore.CompositeCache
OrdinaryDiffEqCore.AutoSwitchCache
OrdinaryDiffEqCore.isautoswitch
OrdinaryDiffEqCore.default_autoswitch
OrdinaryDiffEqCore.unwrap_alg
OrdinaryDiffEqCore.isdefaultalg
OrdinaryDiffEqCore.is_composite_algorithm
OrdinaryDiffEqCore.is_composite_cache
OrdinaryDiffEqCore.is_constant_cache
```

## Algorithm trait functions

Solver sublibraries specialize these to describe their algorithms.

```@docs
OrdinaryDiffEqCore.alg_order
OrdinaryDiffEqCore.alg_maximum_order
OrdinaryDiffEqCore.alg_adaptive_order
OrdinaryDiffEqCore.alg_stability_size
OrdinaryDiffEqCore.alg_extrapolates
OrdinaryDiffEqCore.alg_can_repeat_jac
OrdinaryDiffEqCore.alg_autodiff
OrdinaryDiffEqCore.alg_difftype
OrdinaryDiffEqCore.isfsal
OrdinaryDiffEqCore.fsal_typeof
OrdinaryDiffEqCore.isimplicit
OrdinaryDiffEqCore.isadaptive
OrdinaryDiffEqCore.isdtchangeable
OrdinaryDiffEqCore.ismultistep
OrdinaryDiffEqCore.dt_required
OrdinaryDiffEqCore.uses_uprev
OrdinaryDiffEqCore.has_autodiff
OrdinaryDiffEqCore.has_special_newton_error
OrdinaryDiffEqCore.has_dtnew_modification
OrdinaryDiffEqCore.has_stiff_interpolation
OrdinaryDiffEqCore.allows_null_u0
OrdinaryDiffEqCore.isaposteriori
OrdinaryDiffEqCore.isdiscretealg
OrdinaryDiffEqCore.isdiscretecache
OrdinaryDiffEqCore.isdp8
OrdinaryDiffEqCore.isesdirk
OrdinaryDiffEqCore.isfirk
OrdinaryDiffEqCore.isnewton
OrdinaryDiffEqCore.isWmethod
OrdinaryDiffEqCore.is_mass_matrix_alg
OrdinaryDiffEqCore.issplit
OrdinaryDiffEqCore.only_diagonal_mass_matrix
OrdinaryDiffEqCore.standardtag
OrdinaryDiffEqCore.concrete_jac
OrdinaryDiffEqCore.fac_default_gamma
OrdinaryDiffEqCore.default_linear_interpolation
```

## Order / step-size / autodiff-config accessors

```@docs
OrdinaryDiffEqCore.get_current_alg_order
OrdinaryDiffEqCore.get_current_adaptive_order
OrdinaryDiffEqCore.get_current_alg_autodiff
OrdinaryDiffEqCore.get_chunksize
OrdinaryDiffEqCore._get_fdtype
OrdinaryDiffEqCore._get_fwd_chunksize
OrdinaryDiffEqCore._get_fwd_chunksize_int
OrdinaryDiffEqCore._fixup_ad
OrdinaryDiffEqCore.diffdir
OrdinaryDiffEqCore.error_constant
OrdinaryDiffEqCore.constvalue
OrdinaryDiffEqCore.unitfulvalue
```

## Enums and status types

```@docs
OrdinaryDiffEqCore.COEFFICIENT_MULTISTEP
OrdinaryDiffEqCore.CompiledFloats
OrdinaryDiffEqCore.Convergence
OrdinaryDiffEqCore.DifferentialVarsUndefined
OrdinaryDiffEqCore.DIRK
OrdinaryDiffEqCore.Divergence
OrdinaryDiffEqCore.FastConvergence
OrdinaryDiffEqCore.GLM
OrdinaryDiffEqCore.MethodType
OrdinaryDiffEqCore.NLStatus
OrdinaryDiffEqCore.NORDSIECK_MULTISTEP
OrdinaryDiffEqCore.SlowConvergence
OrdinaryDiffEqCore.TryAgain
OrdinaryDiffEqCore.VerySlowConvergence
```

## Nonlinear solver interface

The public nonlinear-solver algorithms are documented on the OrdinaryDiffEqCore
API page. These core abstractions and W-matrix hooks are the solver-author
extension points.

```@docs
OrdinaryDiffEqCore.AbstractNLSolver
OrdinaryDiffEqCore.AbstractNLSolverAlgorithm
OrdinaryDiffEqCore.AbstractNLSolverCache
OrdinaryDiffEqCore.nlsolve_f
OrdinaryDiffEqCore.get_W
OrdinaryDiffEqCore.set_new_W!
OrdinaryDiffEqCore.set_W_γdt!
OrdinaryDiffEqCore.get_new_W_γdt_cutoff
OrdinaryDiffEqCore.isfirstcall
OrdinaryDiffEqCore.isfirststage
OrdinaryDiffEqCore.isJcurrent
OrdinaryDiffEqCore.resize_J_W!
OrdinaryDiffEqCore.resize_nlsolver!
OrdinaryDiffEqCore.default_nlsolve
```

### OrdinaryDiffEqNonlinearSolve driver hooks

These hooks are a version-controlled, developer-only contract for implicit
solver packages. Cache constructors create an opaque nonlinear solver with
`build_nlsolver`; step implementations mark stages, drive the solve, classify
failure, and query optional workspace capabilities through the remaining
functions. `compute_step!` and `initial_η` are extension points for sibling
packages that implement an `AbstractNLSolver` subtype. The Anderson helpers are
shared implementations, but their concrete workspace types and fields remain
internal.

```@docs
OrdinaryDiffEqNonlinearSolve.build_nlsolver
OrdinaryDiffEqNonlinearSolve.nlsolve!
OrdinaryDiffEqNonlinearSolve.nlsolvefail
OrdinaryDiffEqNonlinearSolve.markfirststage!
OrdinaryDiffEqNonlinearSolve.du_alias_or_new
OrdinaryDiffEqNonlinearSolve.can_smooth_est
OrdinaryDiffEqNonlinearSolve.compute_step!
OrdinaryDiffEqNonlinearSolve.initial_η
OrdinaryDiffEqNonlinearSolve.anderson
OrdinaryDiffEqNonlinearSolve.anderson!
```

## Jacobian / W-matrix / differentiation configuration

Provided by `OrdinaryDiffEqDifferentiation`.

```@docs
OrdinaryDiffEqDifferentiation.build_J_W
OrdinaryDiffEqDifferentiation.build_uf
OrdinaryDiffEqDifferentiation.build_jac_config
OrdinaryDiffEqDifferentiation.build_grad_config
OrdinaryDiffEqDifferentiation.calc_J
OrdinaryDiffEqDifferentiation.calc_J!
OrdinaryDiffEqDifferentiation.calc_tderivative
OrdinaryDiffEqDifferentiation.calc_tderivative!
OrdinaryDiffEqDifferentiation.calc_rosenbrock_differentiation
OrdinaryDiffEqDifferentiation.calc_rosenbrock_differentiation!
OrdinaryDiffEqDifferentiation.jacobian!
OrdinaryDiffEqDifferentiation.jacobian2W!
OrdinaryDiffEqDifferentiation.update_W!
OrdinaryDiffEqDifferentiation.resize_jac_config!
OrdinaryDiffEqDifferentiation.resize_grad_config!
OrdinaryDiffEqDifferentiation.dolinsolve
OrdinaryDiffEqDifferentiation.wrapprecs
OrdinaryDiffEqDifferentiation.default_krylov_warm_start
OrdinaryDiffEqDifferentiation.is_always_new
OrdinaryDiffEqDifferentiation.islinearfunction
OrdinaryDiffEqDifferentiation.issuccess_W
OrdinaryDiffEqDifferentiation.drain_jvp_count!
OrdinaryDiffEqDifferentiation.jvp_counter
OrdinaryDiffEqDifferentiation.set_linear_reltol!
```

## Integrator step, cache construction, and initialization hooks

The cache hierarchy and `@cache` macro are developer API for solver packages
that implement new algorithms. They define internal solver storage, not
application-facing cache objects.

The `StochasticDiffEq*Cache` trio is the SDE/RODE analogue, subtyped by the
StochasticDiffEq.jl solver sublibraries.

```@docs
OrdinaryDiffEqCore.OrdinaryDiffEqCache
OrdinaryDiffEqCore.OrdinaryDiffEqConstantCache
OrdinaryDiffEqCore.OrdinaryDiffEqMutableCache
OrdinaryDiffEqCore.StochasticDiffEqCache
OrdinaryDiffEqCore.StochasticDiffEqConstantCache
OrdinaryDiffEqCore.StochasticDiffEqMutableCache
OrdinaryDiffEqCore.DefaultCache
OrdinaryDiffEqCore.strip_cache
OrdinaryDiffEqCore.@cache
OrdinaryDiffEqCore.alg_cache
OrdinaryDiffEqCore.get_fsalfirstlast
OrdinaryDiffEqCore.get_fresh_jacobian
OrdinaryDiffEqCore.perform_step!
OrdinaryDiffEqCore.apply_step!
```

```@autodocs
Modules = [SciMLBase]
Public = true
Private = false
Filter = x -> x in [
    SciMLBase.postamble!,
    SciMLBase.last_step_failed,
    SciMLBase.check_error,
    SciMLBase.check_error!
]
```

```@docs
OrdinaryDiffEqCore.set_discontinuity
OrdinaryDiffEqCore.increment_accept!
OrdinaryDiffEqCore.increment_reject!
OrdinaryDiffEqCore.increment_nf!
OrdinaryDiffEqCore.ode_determine_initdt
OrdinaryDiffEqCore._determine_initdt
OrdinaryDiffEqCore._ode_init
OrdinaryDiffEqCore._initialize_dae!
OrdinaryDiffEqCore.find_algebraic_vars_eqs
OrdinaryDiffEqCore.get_differential_vars
OrdinaryDiffEqCore.handle_callback_modifiers!
OrdinaryDiffEqCore.resolve_stage_step_limiters
OrdinaryDiffEqCore.trivial_limiter!
OrdinaryDiffEqCore.DEOptions
OrdinaryDiffEqCore.DummyController
```

### Time-stop and saving queues

Custom integrator initialization and stepping loops use these hooks to preserve
the standard `tstops`, `saveat`, derivative-discontinuity, and time-step-bound
semantics. They are versioned developer API, not user-facing solver controls.

```@docs
OrdinaryDiffEqCore.initialize_tstops
OrdinaryDiffEqCore.initialize_saveat
OrdinaryDiffEqCore.initialize_d_discontinuities
OrdinaryDiffEqCore.fix_dt_at_bounds!
OrdinaryDiffEqCore.handle_tstop!
```

## Dense output / interpolation

```@docs
OrdinaryDiffEqCore.OrdinaryDiffEqInterpolation
OrdinaryDiffEqCore.InterpolationData
OrdinaryDiffEqCore.DerivativeOrderNotPossibleError
OrdinaryDiffEqCore.ode_interpolant
OrdinaryDiffEqCore.ode_interpolant!
OrdinaryDiffEqCore.hermite_interpolant
OrdinaryDiffEqCore.hermite_interpolant!
OrdinaryDiffEqCore.interpolation_differential_vars
OrdinaryDiffEqCore.current_interpolant
OrdinaryDiffEqCore.current_extrapolant
OrdinaryDiffEqCore.current_extrapolant!
OrdinaryDiffEqCore.ode_addsteps!
OrdinaryDiffEqCore._ode_interpolant
OrdinaryDiffEqCore._ode_interpolant!
OrdinaryDiffEqCore._ode_addsteps!
```

## Noise-process hooks

Used by SDE/RODE solver sublibraries; no-ops for pure ODEs.

```@docs
OrdinaryDiffEqCore.accept_noise!
OrdinaryDiffEqCore.reject_noise!
OrdinaryDiffEqCore.save_noise!
OrdinaryDiffEqCore.reinit_noise!
OrdinaryDiffEqCore.noise_curt
OrdinaryDiffEqCore.is_noise_saveable
```

## Docstring builders

Helpers that build consistent algorithm docstrings.

```@docs
OrdinaryDiffEqCore.generic_solver_docstring
OrdinaryDiffEqCore.explicit_rk_docstring
OrdinaryDiffEqCore.differentiation_rk_docstring
```

## Stochastic solver extension API

These names are a version-controlled contract for sibling stochastic solver
packages. They are not application-facing solver API; use the stochastic
algorithm pages to select methods for `solve`.

```@docs
StochasticDiffEqCore.AbstractJ
StochasticDiffEqCore.AbstractJCommute
StochasticDiffEqCore.AbstractJDiagonal
StochasticDiffEqCore.DiffEqNLSolveTag
StochasticDiffEqCore.IICommutative
StochasticDiffEqCore.IIFNLSolveFunc
StochasticDiffEqCore.IILevyArea
StochasticDiffEqCore.Ihat2
StochasticDiffEqCore.IteratedIntegralAlgorithm_iip
StochasticDiffEqCore.IteratedIntegralApprox
StochasticDiffEqCore.JCommute_iip
StochasticDiffEqCore.JCommute_oop
StochasticDiffEqCore.JDiagonal_iip
StochasticDiffEqCore.JDiagonal_oop
StochasticDiffEqCore.NLSOLVEJL_SETUP
StochasticDiffEqCore.SDEAlgTypes
StochasticDiffEqCore.SDEIntegrator
StochasticDiffEqCore.SDEOptions
StochasticDiffEqCore.DiffCache
StochasticDiffEqCore.@cache
StochasticDiffEqCore.StochasticCompositeAlgorithm
StochasticDiffEqCore.StochasticCompositeCache
StochasticDiffEqCore.TauLeapingDrift
StochasticDiffEqCore._resolve_rng
StochasticDiffEqCore._sde_init
StochasticDiffEqCore._z_prototype
StochasticDiffEqCore.addat_noise!
StochasticDiffEqCore.alg_cache
StochasticDiffEqCore.alg_compatible
StochasticDiffEqCore.alg_control_rate
StochasticDiffEqCore.alg_mass_matrix_compatible
StochasticDiffEqCore.alg_stability_size
StochasticDiffEqCore.alg_can_repeat_jac
StochasticDiffEqCore.alg_needs_extra_process
StochasticDiffEqCore.calc_threepoint_random
StochasticDiffEqCore.calc_threepoint_random!
StochasticDiffEqCore.calc_twopoint_random
StochasticDiffEqCore.calc_twopoint_random!
StochasticDiffEqCore.concrete_prob
StochasticDiffEqCore.deleteat_noise!
StochasticDiffEqCore.delta_default
StochasticDiffEqCore.determine_chunksize
StochasticDiffEqCore.fill_new_noise_caches!
StochasticDiffEqCore.get_chunksize
StochasticDiffEqCore.get_current_alg_order
StochasticDiffEqCore.get_Jalg
StochasticDiffEqCore.get_iterated_I
StochasticDiffEqCore.get_iterated_I!
StochasticDiffEqCore.is_split_step
StochasticDiffEqCore.resize_noise!
StochasticDiffEqCore.supports_regular_jumps
StochasticDiffEqCore.AutoAlgSwitch
StochasticDiffEqCore.AutoSwitch
StochasticDiffEqCore.unwrap_alg
```

### High-order stochastic tableau construction

```@docs
StochasticDiffEqHighOrder.checkSRAOrder
StochasticDiffEqHighOrder.checkSRIOrder
StochasticDiffEqHighOrder.constructExplicitSKenCarp
StochasticDiffEqHighOrder.constructSKenCarp
StochasticDiffEqHighOrder.constructSOSRA
StochasticDiffEqHighOrder.constructSOSRA2
StochasticDiffEqHighOrder.constructSRA1
StochasticDiffEqHighOrder.constructSRA2
StochasticDiffEqHighOrder.constructSRA3
StochasticDiffEqHighOrder.constructSRIOpt1
StochasticDiffEqHighOrder.constructSRIOpt2
StochasticDiffEqHighOrder.constructSRIW1
StochasticDiffEqHighOrder.constructSRIW2
StochasticDiffEqHighOrder.du_cache
StochasticDiffEqHighOrder.u_cache
StochasticDiffEqHighOrder.user_cache
```

### ESDIRK-IMEX implementation hooks

The following types are developer API for sibling implicit solver packages.
They are not stable application-facing cache layouts.

```@docs
OrdinaryDiffEqSDIRK.ESDIRKIMEXCache
OrdinaryDiffEqSDIRK.ESDIRKIMEXConstantCache
OrdinaryDiffEqSDIRK.ImplicitEulerESDIRKIMEXTableau
```

### Cross-sublibrary cache hooks

These opaque cache types are developer-only extension contracts for sibling
solver packages. Application code must construct algorithms and call `solve`,
rather than depend on cache fields or constructors.

```@docs
OrdinaryDiffEqLowOrderRK.BS3Cache
OrdinaryDiffEqLowOrderRK.BS3ConstantCache
OrdinaryDiffEqLowOrderRK.RK4Cache
OrdinaryDiffEqLowOrderRK.RK4ConstantCache
OrdinaryDiffEqRosenbrock.RosenbrockMutableCache
OrdinaryDiffEqTsit5.Tsit5Cache
OrdinaryDiffEqTsit5.Tsit5ConstantCache
```
