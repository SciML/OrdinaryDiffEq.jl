__precompile__()

module OrdinaryDiffEq

  using Reexport
  @reexport using DiffEqBase

  using Compat

  using MuladdMacro

  using LinearAlgebra

  # Interfaces
  import DiffEqBase: solve!, step!, initialize!

  # Internal utils
  import DiffEqBase: ODE_DEFAULT_NORM, ODE_DEFAULT_ISOUTOFDOMAIN, ODE_DEFAULT_PROG_MESSAGE, ODE_DEFAULT_UNSTABLE_CHECK

  using DiffEqOperators: DiffEqArrayOperator

  import RecursiveArrayTools: chain, recursivecopy!

  using Parameters, GenericSVD, ForwardDiff, RecursiveArrayTools,
        NLsolve, Juno, Roots, DataStructures, DiffEqDiffTools

  import ForwardDiff.Dual

  # Required by temporary fix in not in-place methods with 12+ broadcasts
  # `MVector` is used by Nordsieck forms
  import StaticArrays: SArray, MVector, SVector

  # Integrator Interface
  import DiffEqBase: resize!,deleteat!,addat!,full_cache,user_cache,u_cache,du_cache,
                     resize_non_user_cache!,deleteat_non_user_cache!,addat_non_user_cache!,
                     terminate!,get_du, get_dt,get_proposed_dt,set_proposed_dt!,
                     u_modified!,savevalues!,add_tstop!,add_saveat!,set_reltol!,
                     set_abstol!, postamble!, last_step_failed,
                     isautodifferentiable

  using DiffEqBase: check_error!

  macro tight_loop_macros(ex)
   :($(esc(ex)))
  end

  const CompiledFloats = Union{Float32,Float64}

  include("misc_utils.jl")
  include("algorithms.jl")

  include("caches/basic_caches.jl")
  include("caches/low_order_rk_caches.jl")
  include("caches/high_order_rk_caches.jl")
  include("caches/ssprk_caches.jl")
  include("caches/feagin_caches.jl")
  include("caches/verner_caches.jl")
  include("caches/sdirk_caches.jl")
  include("caches/kencarp_kvaerno_caches.jl")
  include("caches/generic_implicit_caches.jl")
  include("caches/linear_caches.jl")
  include("caches/linear_nonlinear_caches.jl")
  include("caches/symplectic_caches.jl")
  include("caches/rosenbrock_caches.jl")
  include("caches/rkn_caches.jl")
  include("caches/adams_bashforth_moulton_caches.jl")
  include("caches/nordsieck_caches.jl")
  include("caches/bdf_caches.jl")
  include("caches/rkc_caches.jl")
  include("caches/euler_imex_caches.jl")

  include("alg_utils.jl")

  include("tableaus/low_order_rk_tableaus.jl")
  include("tableaus/high_order_rk_tableaus.jl")
  include("tableaus/symplectic_tableaus.jl")
  include("tableaus/verner_tableaus.jl")
  include("tableaus/feagin_tableaus.jl")
  include("tableaus/rosenbrock_tableaus.jl")
  include("tableaus/sdirk_tableaus.jl")
  include("tableaus/rkn_tableaus.jl")
  include("tableaus/rkc_tableaus.jl")

  include("integrators/type.jl")
  include("integrators/integrator_utils.jl")
  include("integrators/integrator_interface.jl")

  include("perform_step/fixed_timestep_perform_step.jl")
  include("perform_step/symplectic_perform_step.jl")
  include("perform_step/rkn_perform_step.jl")
  include("perform_step/split_perform_step.jl")
  include("perform_step/linear_perform_step.jl")
  include("perform_step/iif_perform_step.jl")
  include("perform_step/exponential_rk_perform_step.jl")
  include("perform_step/explicit_rk_perform_step.jl")
  include("perform_step/low_order_rk_perform_step.jl")
  include("perform_step/high_order_rk_perform_step.jl")
  include("perform_step/verner_rk_perform_step.jl")
  include("perform_step/feagin_rk_perform_step.jl")
  include("perform_step/ssprk_perform_step.jl")
  include("perform_step/sdirk_perform_step.jl")
  include("perform_step/kencarp_kvaerno_perform_step.jl")
  include("perform_step/generic_implicit_perform_step.jl")
  include("perform_step/rosenbrock_perform_step.jl")
  include("perform_step/threaded_rk_perform_step.jl")
  include("perform_step/composite_perform_step.jl")
  include("perform_step/adams_bashforth_moulton_perform_step.jl")
  include("perform_step/nordsieck_perform_step.jl")
  include("perform_step/bdf_perform_step.jl")
  include("perform_step/rkc_perform_step.jl")
  include("perform_step/euler_imex_perform_step.jl")

  include("dense/generic_dense.jl")
  include("dense/interpolants.jl")
  include("dense/rosenbrock_interpolants.jl")
  include("dense/stiff_addsteps.jl")
  include("dense/low_order_rk_addsteps.jl")
  include("dense/verner_addsteps.jl")
  include("dense/high_order_rk_addsteps.jl")

  include("derivative_utils.jl")
  include("nlsolve_utils.jl")
  include("nordsieck_utils.jl")
  include("adams_utils.jl")
  include("bdf_utils.jl")
  include("exponential_utils.jl")
  include("rkc_utils.jl")
  include("derivative_wrappers.jl")
  include("iterator_interface.jl")
  include("constants.jl")
  include("callbacks.jl")
  include("composite_solution.jl")
  include("solve.jl")
  include("initdt.jl")
  include("interp_func.jl")
  include("composite_algs.jl")

  #General Functions
  export solve, solve!, init, step!

  export OrdinaryDiffEqAlgorithm

  #Callback Necessary
  export ode_addsteps!, ode_interpolant,
        terminate!, savevalues!, copyat_or_push!, isfsal

  export constructDormandPrince

  # Reexport the Alg Types

  export Discrete, FunctionMap, Euler, Heun, Ralston, Midpoint, SSPRK22,
         SSPRK33, SSPRK53, SSPRK63, SSPRK73, SSPRK83, SSPRK432, SSPRK932,
         SSPRK54, SSPRK104, RK4, ExplicitRK, OwrenZen3, OwrenZen4, OwrenZen5,
         BS3, BS5, CarpenterKennedy2N54,
         DP5, DP5Threaded, Tsit5, DP8, Vern6, Vern7, Vern8, TanYam7, TsitPap8,
         Vern9,Feagin10, Feagin12, Feagin14, CompositeAlgorithm, Anas5

  export ImplicitEuler, ImplicitMidpoint, Trapezoid, TRBDF2, SDIRK2, Kvaerno3,
         KenCarp3, Cash4, Hairer4, Hairer42, SSPSDIRK2, Kvaerno4, Kvaerno5,
         KenCarp4, KenCarp5

  export GenericImplicitEuler, GenericTrapezoid

  export LinearImplicitEuler, MidpointSplitting

  export Rosenbrock23, Rosenbrock32, RosShamp4, Veldd4, Velds4, GRK4T, GRK4A,
         Ros4LStab, ROS3P, Rodas3, Rodas4, Rodas42, Rodas4P, Rodas5

  export GenericIIF1, GenericIIF2

  export LawsonEuler, NorsettEuler, ETD1, ETDRK2, ETDRK3, ETDRK4, HochOst4, Exp4, EPIRK4s3A, EPIRK4s3B,
         EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2, ETD2

  export SymplecticEuler, VelocityVerlet, VerletLeapfrog, PseudoVerletLeapfrog,
         McAte2, Ruth3, McAte3, CandyRoz4, McAte4, McAte42, McAte5,
         CalvoSanz4, Yoshida6, KahanLi6, McAte8, KahanLi8, SofSpa10

  export SplitEuler

  export Nystrom4, Nystrom4VelocityIndependent, Nystrom5VelocityIndependent,
         IRKN3, IRKN4, DPRKN6, DPRKN8, DPRKN12, ERKN4, ERKN5

  export ROCK2

  export AB3, AB4, AB5, ABM32, ABM43, ABM54

  export VCAB3, VCAB4, VCAB5, VCABM3, VCABM4, VCABM5

  export VCABM

  export IMEXEuler, CNAB2, CNLF2

  export AN5, JVODE, JVODE_Adams

  export ABDF2, QNDF1, QBDF1, QNDF2, QBDF2, QNDF

  export AutoSwitch, AutoTsit5, AutoDP5,
         AutoVern6, AutoVern7, AutoVern8, AutoVern9
end # module
