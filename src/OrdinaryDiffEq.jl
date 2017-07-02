__precompile__()

module OrdinaryDiffEq

  using Reexport
  @reexport using DiffEqBase

  using Compat

  # Interfaces
  import DiffEqBase: solve, solve!, init, step!, build_solution, initialize!

  # Internal utils
  import DiffEqBase: realtype, ODE_DEFAULT_NORM, ODE_DEFAULT_ISOUTOFDOMAIN, ODE_DEFAULT_PROG_MESSAGE, ODE_DEFAULT_UNSTABLE_CHECK

  import RecursiveArrayTools: chain

  using Parameters, GenericSVD, ForwardDiff, RecursiveArrayTools,
        NLsolve, Juno, Calculus, Roots, DataStructures

  import Base: linspace

  import Base: start, next, done, eltype

  import ForwardDiff.Dual

  # Integrator Interface
  import DiffEqBase: resize!,deleteat!,addat!,full_cache,user_cache,u_cache,du_cache,
                     resize_non_user_cache!,deleteat_non_user_cache!,addat_non_user_cache!,
                     terminate!,get_du, get_dt,get_proposed_dt,set_proposed_dt!,
                     u_modified!,savevalues!,add_tstop!,add_saveat!,set_reltol!,
                     set_abstol!

  macro tight_loop_macros(ex)
   :($(esc(ex)))
  end

  const CompiledFloats = Union{Float32,Float64}

  include("misc_utils.jl")
  include("algorithms.jl")
  include("alg_utils.jl")
  include("caches.jl")
  include("tableaus/low_order_rk_tableaus.jl")
  include("tableaus/high_order_rk_tableaus.jl")
  include("tableaus/symplectic_tableaus.jl")
  include("tableaus/verner_tableaus.jl")
  include("tableaus/feagin_tableaus.jl")
  include("integrators/type.jl")
  include("integrators/integrator_utils.jl")
  include("integrators/fixed_timestep_integrators.jl")
  include("integrators/symplectic_integrators.jl")
  include("integrators/split_integrators.jl")
  include("integrators/iif_integrators.jl")
  include("integrators/exponential_rk_integrators.jl")
  include("integrators/explicit_rk_integrator.jl")
  include("integrators/low_order_rk_integrators.jl")
  include("integrators/high_order_rk_integrators.jl")
  include("integrators/verner_rk_integrators.jl")
  include("integrators/feagin_rk_integrators.jl")
  include("integrators/ssprk_integrators.jl")
  include("integrators/implicit_integrators.jl")
  include("integrators/rosenbrock_integrators.jl")
  include("integrators/threaded_rk_integrators.jl")
  include("integrators/integrator_interface.jl")
  include("integrators/composite_integrator.jl")
  include("iterator_interface.jl")
  include("constants.jl")
  include("callbacks.jl")
  include("composite_solution.jl")
  include("solve.jl")
  include("initdt.jl")
  include("interp_func.jl")
  include("dense/generic_dense.jl")
  include("dense/interpolants.jl")
  include("dense/stiff_addsteps.jl")
  include("dense/low_order_rk_addsteps.jl")
  include("dense/verner_addsteps.jl")
  include("dense/high_order_rk_addsteps.jl")
  include("derivative_wrappers.jl")

  #General Functions
  export solve, solve!, init, step!

  #Callback Necessary
  export ode_addsteps!, ode_interpolant,
        terminate!, savevalues!, copyat_or_push!, isfsal

  export constructDP5, constructVern6, constructDP8, constructDormandPrince, constructFeagin10,
        constructFeagin12, constructFeagin14

  # Reexport the Alg Types

  export OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqCompositeAlgorithm,
        Discrete, FunctionMap, Euler, Midpoint, SSPRK22, SSPRK33, SSPRK432, SSPRK104, RK4, ExplicitRK, BS3, BS5,
        DP5, DP5Threaded, Tsit5, DP8, Vern6, Vern7, Vern8, TanYam7, TsitPap8, Vern9, ImplicitEuler,
        Trapezoid, Rosenbrock23, Rosenbrock32, Feagin10, Feagin12, Feagin14,
        CompositeAlgorithm

  export IIF1, IIF2

  export LawsonEuler, NorsettEuler

  export SymplecticEuler, VelocityVerlet, Ruth3, McAte3

  export SplitEuler
end # module
