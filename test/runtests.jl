using OrdinaryDiffEq
using Test, LinearAlgebra, Statistics

const CPU_FLOPS = peakflops()
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

tic()
if group == "All" || group == "Interface"
  # Fail
  @time include("discrete_algorithm_test.jl")
  @time include("ode/ode_tstops_tests.jl")
  @time include("ode/ode_backwards_test.jl")
  # Fail DevTool
  @time include("ode/ode_initdt_tests.jl")
  @time include("mass_matrix_tests.jl")
  @time include("differentiation_traits_tests.jl")
  @time include("ode/ode_saveat_tests.jl")
  # Fail
  @time include("ode/ode_saveidxs_tests.jl")
  @time include("static_array_tests.jl")
  # Fail
  @time include("data_array_test.jl")
  @time include("umodified_test.jl")
  @time include("composite_algorithm_test.jl")
  @time include("complex_tests.jl")
  # Fail ???
  @time include("stiffness_detection_test.jl")
  @time include("export_tests.jl")
  @time include("utility_tests.jl")
end

if group == "All" || group == "Integrators"
    @time @testset "Reinit Tests" begin include("reinit_test.jl") end
    # Fail
    @time @testset "Events Tests" begin include("ode/ode_event_tests.jl") end
    @time @testset "Cache Tests" begin include("ode/ode_cache_tests.jl") end
    # TODO
    @time @testset "Iterator Tests" begin include("iterator_tests.jl") end
    @time @testset "Integrator Interface Tests" begin include("integrator_interface_tests.jl") end
    # Fail
    @time @testset "Add Steps Tests" begin include("ode/ode_add_steps_tests.jl") end
    @time @testset "Error Check Tests" begin include("check_error.jl") end
end

if !is_APPVEYOR && ( group == "All" || group == "Regression" )
    @time @testset "Linear Tests" begin include("ode/ode_twodimlinear_tests.jl") end
    @time @testset "Dense Tests" begin include("ode/ode_dense_tests.jl") end
end

if !is_APPVEYOR && ( group == "All" || group == "AlgConvergence_I" )
    # ~ 250 s
    @time @testset "Partitioned Methods Tests" begin include("partitioned_methods_tests.jl") end
    # ~ 400 s
    @time @testset "Convergence Tests" begin include("ode/ode_convergence_tests.jl") end
    # ~ 2 s
    @time @testset "Adams Variable Coefficients Tests" begin include("ode/adams_tests.jl") end
    # ~ 50 s
    @time @testset "Nordsieck Tests" begin include("ode/nordsieck_tests.jl") end
    #@time @testset "Linear Methods Tests" begin include("linear_method_tests.jl") end
    # ~ 170 s
    @time @testset "SSPRK Tests" begin include("ode/ode_ssprk_tests.jl") end
    # ~ 25 s
    @time @testset "OwrenZen Tests" begin include("owrenzen_tests.jl") end
end

if !is_APPVEYOR && ( group == "All" || group == "AlgConvergence_II" )
    # ~ 110 s
    @time @testset "Split Methods Tests" begin include("split_methods_tests.jl") end
    # ~ 550 s
    @time @testset "Rosenbrock Tests" begin include("ode/ode_rosenbrock_tests.jl") end
    # ~ 40 s
    @time @testset "Linear-Nonlinear Methods Tests" begin include("linear_nonlinear_convergence_tests.jl") end
    # ~ 140 s
    @time @testset "Linear-Nonlinear Krylov Methods Tests" begin include("linear_nonlinear_krylov_tests.jl") end
end

toc()
