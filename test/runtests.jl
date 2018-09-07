if haskey(ENV,"GROUP")
    group = ENV["GROUP"]
else
    group = "All"
end

is_APPVEYOR = ( Sys.iswindows() && haskey(ENV,"APPVEYOR") )

#Start Test Script

@time begin
if group == "All" || group == "Interface"
  for file in ["discrete_algorithm_test.jl", "ode/ode_tstops_tests.jl",
               "ode/ode_backwards_test.jl", "ode/ode_initdt_tests.jl",
               "mass_matrix_tests.jl", "differentiation_traits_tests.jl",
               "ode/ode_saveat_tests.jl", "ode/ode_saveidxs_tests.jl",
               "static_array_tests.jl", "data_array_test.jl",
               "umodified_test.jl", "composite_algorithm_test.jl",
               "complex_tests.jl", "stiffness_detection_test.jl",
               "export_tests.jl", "utility_tests.jl"]
    @eval module $(gensym())
      using Test
      @time include($file)
    end
  end
end

if group == "All" || group == "Integrators"
  for file in ["reinit_test.jl", "ode/ode_event_tests.jl",
               "ode/ode_cache_tests.jl", "iterator_tests.jl",
               "integrator_interface_tests.jl", "ode/ode_add_steps_tests.jl",
               "check_error.jl"]
    @eval module $(gensym())
      using Test
      @time include($file)
    end
  end
end

if !is_APPVEYOR && ( group == "All" || group == "Regression" )
  for file in ["ode/ode_twodimlinear_tests.jl", "ode/ode_dense_tests.jl"]
    @eval module $(gensym())
      using Test
      @time include($file)
    end
  end
end

if !is_APPVEYOR && ( group == "All" || group == "AlgConvergence_I" )
  for file in ["partitioned_methods_tests.jl", "ode/ode_convergence_tests.jl",
               "ode/adams_tests.jl", "ode/nordsieck_tests.jl",
               #@time @testset "Linear Methods Tests" begin include("linear_method_tests.jl") end
               "ode/ode_ssprk_tests.jl", "owrenzen_tests.jl",
               "ode/rkc_tests.jl"]
    @eval module $(gensym())
      using Test
      @time include($file)
    end
  end
end

if !is_APPVEYOR && ( group == "All" || group == "AlgConvergence_II" )
  for file ["split_methods_tests.jl", "ode/ode_rosenbrock_tests.jl",
            "linear_nonlinear_convergence_tests.jl", "linear_nonlinear_krylov_tests.jl"]
    @eval module $(gensym())
      using Test
      @time include($file)
    end
  end
end
end # @time
