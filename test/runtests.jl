using OrdinaryDiffEq
using Test, LinearAlgebra, Statistics

const LONGER_TESTS = false

const CACHE_TEST_ALGS = [Euler(),Midpoint(),RK4(),SSPRK22(),SSPRK33(),SSPRK53(),
  SSPRK63(),SSPRK73(),SSPRK83(),SSPRK432(),SSPRK932(),SSPRK54(),SSPRK104(),CarpenterKennedy2N54(),
  BS3(),BS5(),DP5(),DP5Threaded(),DP8(),Feagin10(),Feagin12(),Feagin14(),TanYam7(),
  Tsit5(),TsitPap8(),Vern6(),Vern7(),Vern8(),Vern9(),OwrenZen3(),OwrenZen4(),OwrenZen5()]

if haskey(ENV,"GROUP")
    group = ENV["GROUP"]
else
    group = "All"
end

is_APPVEYOR = ( Sys.iswindows() && haskey(ENV,"APPVEYOR") )

#Start Test Script

@time begin
if group == "All" || group == "Interface"
  @time include("discrete_algorithm_test.jl")
  @time include("ode/ode_tstops_tests.jl")
  @time include("ode/ode_backwards_test.jl")
  @time include("ode/ode_initdt_tests.jl")
  @time include("mass_matrix_tests.jl")
  @time include("differentiation_traits_tests.jl")
  @time include("ode/ode_saveat_tests.jl")
  @time include("ode/ode_saveidxs_tests.jl")
  @time include("static_array_tests.jl")
  @time include("data_array_test.jl")
  @time include("umodified_test.jl")
  @time include("composite_algorithm_test.jl")
  @time include("complex_tests.jl")
  @time include("stiffness_detection_test.jl")
  @time include("export_tests.jl")
  @time include("utility_tests.jl")
end

if group == "All" || group == "Integrators"
  @time include("reinit_test.jl")
  @time include("ode/ode_event_tests.jl")
  @time include("ode/ode_cache_tests.jl")
  @time include("iterator_tests.jl")
  @time include("integrator_interface_tests.jl")
  @time include("ode/ode_add_steps_tests.jl")
  @time include("check_error.jl")
end

if !is_APPVEYOR && ( group == "All" || group == "Regression" )
  @time include("ode/ode_twodimlinear_tests.jl")
  @time include("ode/ode_dense_tests.jl")
end

if !is_APPVEYOR && ( group == "All" || group == "AlgConvergence_I" )
  # ~ 250 s
  @time include("partitioned_methods_tests.jl")
  # ~ 400 s
  @time include("ode/ode_convergence_tests.jl")
  # ~ 2 s
  @time include("ode/adams_tests.jl")
  # ~ 50 s
  @time include("ode/nordsieck_tests.jl")
  #@time @testset "Linear Methods Tests" begin include("linear_method_tests.jl") end
  # ~ 170 s
  @time include("ode/ode_ssprk_tests.jl")
  # ~ 25 s
  @time include("owrenzen_tests.jl")
  @time include("ode/rkc_tests.jl")
end

if !is_APPVEYOR && ( group == "All" || group == "AlgConvergence_II" )
  # ~ 110 s
  @time include("split_methods_tests.jl")
  # ~ 550 s
  @time include("ode/ode_rosenbrock_tests.jl")
  # ~ 40 s
  @time include("linear_nonlinear_convergence_tests.jl")
  # ~ 140 s
  @time include("linear_nonlinear_krylov_tests.jl")
end
end # @time
