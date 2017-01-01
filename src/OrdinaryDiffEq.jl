__precompile__()

module OrdinaryDiffEq

  using DiffEqBase
  import DiffEqBase: solve, solve!, init
  using Parameters, GenericSVD, ForwardDiff, InplaceOps, RecursiveArrayTools,
        Ranges, NLsolve, RecipesBase, Juno, Calculus, Roots, DataStructures

  import Base: linspace
  import ForwardDiff.Dual

  #Constants

  const TEST_FLOPS_CUTOFF = 1e10

  include("misc_utils.jl")
  include("algorithms.jl")
  include("alg_utils.jl")
  include("integrators/integrator_utils.jl")
  include("integrators/fixed_timestep_integrators.jl")
  include("integrators/explicit_rk_integrator.jl")
  include("integrators/low_order_rk_integrators.jl")
  include("integrators/high_order_rk_integrators.jl")
  include("integrators/verner_rk_integrators.jl")
  include("integrators/feagin_rk_integrators.jl")
  include("integrators/implicit_integrators.jl")
  include("integrators/rosenbrock_integrators.jl")
  include("integrators/threaded_rk_integrators.jl")
  include("constants.jl")
  include("callbacks.jl")
  include("caches.jl")
  include("solve.jl")
  include("initdt.jl")
  include("dense.jl")
  include("integrators/unrolled_tableaus.jl")


  #General Functions
  export solve, solve!, init

  #Callback Necessary
  export ode_addsteps!, ode_interpolant, DIFFERENTIALEQUATIONSJL_SPECIALDENSEALGS,
        @ode_callback, @ode_event, @ode_change_cachesize, @ode_change_deleteat,
        @ode_terminate, ode_savevalues!, copyat_or_push!, isspecialdense,
        ode_postable!,isfsal

  export constructDP5, constructVern6, constructDP8, constructDormandPrince, constructFeagin10,
        constructFeagin12, constructFeagin14

  # Reexport the Alg Types

  export OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveAlgorithm,
        Euler, Midpoint, RK4, ExplicitRK, BS3, BS5, DP5, DP5Threaded, Tsit5,
        DP8, Vern6, Vern7, Vern8, TanYam7, TsitPap8, Vern9, ImplicitEuler,
        Trapezoid, Rosenbrock23, Rosenbrock32, Feagin10, Feagin12, Feagin14
end # module
