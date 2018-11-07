using SafeTestsets
const LONGER_TESTS = false

if haskey(ENV,"GROUP")
    group = ENV["GROUP"]
else
    group = "All"
end

is_APPVEYOR = ( Sys.iswindows() && haskey(ENV,"APPVEYOR") )

#Start Test Script

@time begin
if group == "All" || group == "Interface"
  @time @safetestset "Discrete Algorithm Tests" begin include("discrete_algorithm_test.jl") end
  @time @safetestset "Tstops Tests" begin include("ode/ode_tstops_tests.jl") end
  @time @safetestset "Backwards Tests" begin include("ode/ode_backwards_test.jl") end
  @time @safetestset "Initdt Tests" begin include("ode/ode_initdt_tests.jl") end
  @time @safetestset "Mass Matrix Tests" begin include("mass_matrix_tests.jl") end
  @time @safetestset "Differentiation Trait Tests" begin include("differentiation_traits_tests.jl") end
  @time @safetestset "Inf Tests" begin include("inf_handling.jl") end
  @time @safetestset "saveat Tests" begin include("ode/ode_saveat_tests.jl") end
  @time @safetestset "save_idxs Tests" begin include("ode/ode_saveidxs_tests.jl") end
  @time @safetestset "Static Array Tests" begin include("static_array_tests.jl") end
  @time @safetestset "Data Array Tests" begin include("data_array_test.jl") end
  @time @safetestset "u_modifed Tests" begin include("umodified_test.jl") end
  @time @safetestset "Composite Algorithm Tests" begin include("composite_algorithm_test.jl") end
  @time @safetestset "Complex Tests" begin include("complex_tests.jl") end
  @time @safetestset "Stiffness Detection Tests" begin include("stiffness_detection_test.jl") end
  @time @safetestset "Composite Interpolation Tests" begin include("composite_interpolation.jl") end
  @time @safetestset "Export tests" begin include("export_tests.jl") end
  @time @safetestset "Derivative Utilities Tests" begin include("utility_tests.jl") end
end

if group == "All" || group == "Integrators"
  @time @safetestset "Reinit Tests" begin include("reinit_test.jl") end
  @time @safetestset "Events Tests" begin include("ode/ode_event_tests.jl") end
  @time @safetestset "Cache Tests" begin include("ode/ode_cache_tests.jl") end
  @time @safetestset "Iterator Tests" begin include("iterator_tests.jl") end
  @time @safetestset "Integrator Interface Tests" begin include("integrator_interface_tests.jl") end
  @time @safetestset "Add Steps Tests" begin include("ode/ode_add_steps_tests.jl") end
  @time @safetestset "Error Check Tests" begin include("check_error.jl") end
end

if !is_APPVEYOR && ( group == "All" || group == "Regression" )
  @time @safetestset "Linear Tests" begin include("ode/ode_twodimlinear_tests.jl") end
  @time @safetestset "Dense Tests" begin include("ode/ode_dense_tests.jl") end
end

if !is_APPVEYOR && ( group == "All" || group == "AlgConvergence_I" )
  # ~ 250 s
  @time @safetestset "Partitioned Methods Tests" begin include("partitioned_methods_tests.jl") end
  # ~ 400 s
  @time @safetestset "Convergence Tests" begin include("ode/ode_convergence_tests.jl") end
  # ~ 2 s
  @time @safetestset "Adams Variable Coefficients Tests" begin include("ode/adams_tests.jl") end
  # ~ 50 s
  @time @safetestset "Nordsieck Tests" begin include("ode/nordsieck_tests.jl") end
  @time @safetestset "Linear Methods Tests" begin include("linear_method_tests.jl") end
  # ~ 170 s
  @time @safetestset "SSPRK Tests" begin include("ode/ode_ssprk_tests.jl") end
  # ~ 25 s
  @time @safetestset "OwrenZen Tests" begin include("owrenzen_tests.jl") end
  @time @safetestset "Runge-Kutta-Chebyshev Tests" begin include("ode/rkc_tests.jl") end
end

if !is_APPVEYOR && ( group == "All" || group == "AlgConvergence_II" )
  # ~ 110 s
  @time @safetestset "Split Methods Tests" begin include("split_methods_tests.jl") end
  # ~ 550 s
  @time @safetestset "Rosenbrock Tests" begin include("ode/ode_rosenbrock_tests.jl") end
  # ~ 40 s
  @time @safetestset "Linear-Nonlinear Methods Tests" begin include("linear_nonlinear_convergence_tests.jl") end
  # ~ 140 s
  @time @safetestset "Linear-Nonlinear Krylov Methods Tests" begin include("linear_nonlinear_krylov_tests.jl") end
end
end # @time
