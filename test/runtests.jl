const LONGER_TESTS = false

if haskey(ENV,"GROUP")
    group = ENV["GROUP"]
else
    group = "All"
end

is_APPVEYOR = ( Sys.iswindows() && haskey(ENV,"APPVEYOR") )

macro testset_module(name, expr)
  if name isa String
    mod = gensym()
    testname = name
  elseif name isa Expr && name.head == :(=) && length(name.args) == 2
    mod, testname = name.args[1]
  else
    error("""
          Use `@testset_module` like the following:
          @time @testset_module "Benchmark Tests" begin include("benchmark_tests.jl") end
          @time @testset_module BenchmarkTests = "Benchmark Tests" begin include("benchmark_tests.jl") end
          """)
  end
  quote
    @eval module $mod
      using Test
      @testset $testname $expr
    end
    nothing
  end
end

#Start Test Script

@time begin
if group == "All" || group == "Interface"
  @time @testset_module "Discrete Algorithm Tests" begin include("discrete_algorithm_test.jl") end
  @time @testset_module "Tstops Tests" begin include("ode/ode_tstops_tests.jl") end
  @time @testset_module "Backwards Tests" begin include("ode/ode_backwards_test.jl") end
  @time @testset_module "Initdt Tests" begin include("ode/ode_initdt_tests.jl") end
  @time @testset_module "Mass Matrix Tests" begin include("mass_matrix_tests.jl") end
  @time @testset_module "Differentiation Trait Tests" begin include("differentiation_traits_tests.jl") end
  @time @testset_module "saveat Tests" begin include("ode/ode_saveat_tests.jl") end
  @time @testset_module "save_idxs Tests" begin include("ode/ode_saveidxs_tests.jl") end
  @time @testset_module "Static Array Tests" begin include("static_array_tests.jl") end
  @time @testset_module "Data Array Tests" begin include("data_array_test.jl") end
  @time @testset_module "u_modifed Tests" begin include("umodified_test.jl") end
  @time @testset_module "Composite Algorithm Tests" begin include("composite_algorithm_test.jl") end
  @time @testset_module "Complex Tests" begin include("complex_tests.jl") end
  @time @testset_module "Stiffness Detection Tests" begin include("stiffness_detection_test.jl") end
  @time @testset_module "Export tests" begin include("export_tests.jl") end
  @time @testset_module "Derivative Utilities Tests" begin include("utility_tests.jl") end
end

if group == "All" || group == "Integrators"
  @time @testset_module "Reinit Tests" begin include("reinit_test.jl") end
  @time @testset_module "Events Tests" begin include("ode/ode_event_tests.jl") end
  @time @testset_module "Cache Tests" begin include("ode/ode_cache_tests.jl") end
  @time @testset_module "Iterator Tests" begin include("iterator_tests.jl") end
  @time @testset_module "Integrator Interface Tests" begin include("integrator_interface_tests.jl") end
  @time @testset_module "Add Steps Tests" begin include("ode/ode_add_steps_tests.jl") end
  @time @testset_module "Error Check Tests" begin include("check_error.jl") end
end

if !is_APPVEYOR && ( group == "All" || group == "Regression" )
  @time @testset_module "Linear Tests" begin include("ode/ode_twodimlinear_tests.jl") end
  @time @testset_module "Dense Tests" begin include("ode/ode_dense_tests.jl") end
end

if !is_APPVEYOR && ( group == "All" || group == "AlgConvergence_I" )
  # ~ 250 s
  @time @testset_module "Partitioned Methods Tests" begin include("partitioned_methods_tests.jl") end
  # ~ 400 s
  @time @testset_module "Convergence Tests" begin include("ode/ode_convergence_tests.jl") end
  # ~ 2 s
  @time @testset_module "Adams Variable Coefficients Tests" begin include("ode/adams_tests.jl") end
  # ~ 50 s
  @time @testset_module "Nordsieck Tests" begin include("ode/nordsieck_tests.jl") end
  #@time @testset "Linear Methods Tests" begin include("linear_method_tests.jl") end
  # ~ 170 s
  @time @testset_module "SSPRK Tests" begin include("ode/ode_ssprk_tests.jl") end
  # ~ 25 s
  @time @testset_module "OwrenZen Tests" begin include("owrenzen_tests.jl") end
  @time @testset_module "Runge-Kutta-Chebyshev Tests" begin include("ode/rkc_tests.jl") end
end

if !is_APPVEYOR && ( group == "All" || group == "AlgConvergence_II" )
  # ~ 110 s
  @time @testset_module "Split Methods Tests" begin include("split_methods_tests.jl") end
  # ~ 550 s
  @time @testset_module "Rosenbrock Tests" begin include("ode/ode_rosenbrock_tests.jl") end
  # ~ 40 s
  @time @testset_module "Linear-Nonlinear Methods Tests" begin include("linear_nonlinear_convergence_tests.jl") end
  # ~ 140 s
  @time @testset_module "Linear-Nonlinear Krylov Methods Tests" begin include("linear_nonlinear_krylov_tests.jl") end
end
end # @time
