using OrdinaryDiffEq, DiffEqProblemLibrary, DiffEqDevTools

if VERSION <= v"0.5+"
  using ODEInterfaceDiffEq
end
using Base.Test

const CPU_FLOPS = peakflops()
const LONGER_TESTS = false
#Start Test Script

tic()
@time @testset "Linear Tests" begin include("ode/ode_twodimlinear_tests.jl") end
@time @testset "Convergence Tests" begin include("ode/ode_convergence_tests.jl") end
@time @testset "Adaptive Tests" begin include("ode/ode_adaptive_tests.jl") end
@time @testset "Tstops Tests" begin include("ode/ode_tstops_tests.jl") end
@time @testset "Backwards Tests" begin include("ode/ode_backwards_test.jl") end
(LONGER_TESTS) && @time @testset "Unrolled Tests" begin include("ode/ode_unrolled_comparison_tests.jl") end
@time @testset "Initial Dt Tests" begin include("ode/ode_initdt_tests.jl") end
@time @testset "Rosenbrock Tests" begin include("ode/ode_rosenbrock_tests.jl") end
@time @testset "Dense Tests" begin include("ode/ode_dense_tests.jl") end
@time @testset "In-Place Tests" begin include("ode/ode_inplace_tests.jl") end
@time @testset "Events Tests" begin include("ode/ode_event_tests.jl") end
#@time @testset "Cache Tests" begin include("ode/ode_cache_tests.jl") end
@time @testset "saveat Tests" begin include("ode/ode_saveat_tests.jl") end
(LONGER_TESTS) && @time @testset "Feagin Tests" begin include("ode/ode_feagin_tests.jl") end
@time @testset "Number Type Tests" begin include("ode/ode_numbertype_tests.jl") end
@time @testset "Ndim Complex Tests" begin include("ode/ode_ndim_complex_tests.jl") end

(LONGER_TESTS) && (@time @testset "Units Tests" begin include("units_tests.jl") end) # Too long for AppVeyor

toc()
