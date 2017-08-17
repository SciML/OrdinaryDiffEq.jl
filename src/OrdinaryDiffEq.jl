__precompile__()

module OrdinaryDiffEq

  using Reexport
  @reexport using DiffEqBase

  using Compat

  using MuladdMacro

  # Interfaces
  import DiffEqBase: solve, solve!, init, step!, build_solution, initialize!

  # Internal utils
  import DiffEqBase: ODE_DEFAULT_NORM, ODE_DEFAULT_ISOUTOFDOMAIN, ODE_DEFAULT_PROG_MESSAGE, ODE_DEFAULT_UNSTABLE_CHECK

  import RecursiveArrayTools: chain

  using Parameters, GenericSVD, ForwardDiff, RecursiveArrayTools,
        NLsolve, Juno, Calculus, Roots, DataStructures

  import Base: linspace

  import Base: start, next, done, eltype

  import ForwardDiff.Dual

  # Required by temporary fix in not in-place methods with 12+ broadcasts
  import StaticArrays: SArray

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

  include("caches/basic_caches.jl")
  include("caches/low_order_rk_caches.jl")
  include("caches/high_order_rk_caches.jl")
  include("caches/ssprk_caches.jl")
  include("caches/feagin_caches.jl")
  include("caches/verner_caches.jl")
  include("caches/sdirk_caches.jl")
  include("caches/generic_implicit_caches.jl")
  include("caches/linear_caches.jl")
  include("caches/linear_nonlinear_caches.jl")
  include("caches/symplectic_caches.jl")
  include("caches/rosenbrock_caches.jl")
  include("caches/rkn_caches.jl")

  include("tableaus/low_order_rk_tableaus.jl")
  include("tableaus/high_order_rk_tableaus.jl")
  include("tableaus/symplectic_tableaus.jl")
  include("tableaus/verner_tableaus.jl")
  include("tableaus/feagin_tableaus.jl")
  include("tableaus/rosenbrock_tableaus.jl")
  include("tableaus/sdirk_tableaus.jl")
  include("tableaus/rkn_tableaus.jl")

  include("integrators/type.jl")
  include("integrators/integrator_utils.jl")
  include("integrators/fixed_timestep_integrators.jl")
  include("integrators/symplectic_integrators.jl")
  include("integrators/rkn_integrators.jl")
  include("integrators/split_integrators.jl")
  include("integrators/linear_integrators.jl")
  include("integrators/iif_integrators.jl")
  include("integrators/exponential_rk_integrators.jl")
  include("integrators/explicit_rk_integrator.jl")
  include("integrators/low_order_rk_integrators.jl")
  include("integrators/high_order_rk_integrators.jl")
  include("integrators/verner_rk_integrators.jl")
  include("integrators/feagin_rk_integrators.jl")
  include("integrators/ssprk_integrators.jl")
  include("integrators/sdirk_integrators.jl")
  include("integrators/generic_implicit_integrators.jl")
  include("integrators/rosenbrock_integrators.jl")
  include("integrators/threaded_rk_integrators.jl")
  include("integrators/integrator_interface.jl")
  include("integrators/composite_integrator.jl")

  include("dense/generic_dense.jl")
  include("dense/interpolants.jl")
  include("dense/rosenbrock_interpolants.jl")
  include("dense/stiff_addsteps.jl")
  include("dense/low_order_rk_addsteps.jl")
  include("dense/verner_addsteps.jl")
  include("dense/high_order_rk_addsteps.jl")
  include("derivative_wrappers.jl")

  include("iterator_interface.jl")
  include("constants.jl")
  include("callbacks.jl")
  include("composite_solution.jl")
  include("solve.jl")
  include("initdt.jl")
  include("interp_func.jl")

  #General Functions
  export solve, solve!, init, step!

  #Callback Necessary
  export ode_addsteps!, ode_interpolant,
        terminate!, savevalues!, copyat_or_push!, isfsal

  export constructDP5, constructVern6, constructDP8,
         constructDormandPrince, constructFeagin10,
         constructFeagin12, constructFeagin14

  # Reexport the Alg Types

  export Discrete, FunctionMap, Euler, Heun, Ralston, Midpoint, SSPRK22,
         SSPRK33, SSPRK432, SSPRK104, RK4, ExplicitRK,
         OwrenZen3, OwrenZen4, OwrenZen5, BS3, BS5,
         DP5, DP5Threaded, Tsit5, DP8, Vern6, Vern7, Vern8, TanYam7, TsitPap8,
         Vern9,Feagin10, Feagin12, Feagin14, CompositeAlgorithm

  export ImplicitEuler, Trapezoid, TRBDF2, SDIRK2, Kvaerno3, KenCarp3, Cash4,
         Hairer4, Hairer42, SSPSDIRK2, Kvaerno4, Kvaerno5, KenCarp4, KenCarp5

  export GenericImplicitEuler, GenericTrapezoid

  export LinearImplicitEuler, StrangSplitting

  export Rosenbrock23, Rosenbrock32, RosShamp4, Veldd4, Velds4, GRK4T, GRK4A,
         Ros4LStab, ROS3P, Rodas3, Rodas4, Rodas42, Rodas4P, Rodas5

  export IIF1, IIF2

  export LawsonEuler, NorsettEuler

  export SymplecticEuler, VelocityVerlet, VerletLeapfrog, PseudoVerletLeapfrog,
         McAte2, Ruth3, McAte3, CandyRoz4, McAte4, McAte42, McAte5,
         CalvoSanz4, Yoshida6, KahanLi6, McAte8, KahanLi8, SofSpa10

  export SplitEuler

  export Nystrom4, Nystrom4VelocityIndependent, Nystrom5VelocityIndependent,
         IRKN3, IRKN4, DPRKN6, DPRKN8
end # module
