using OrdinaryDiffEq, DiffEqProblemLibrary, DiffEqDevTools
using Base.Test

const CPU_FLOPS = peakflops()
const TEST_PLOT = false
const LONGER_TESTS = true #Requires JLD
const TEST_CONDITIONAL_DEPS = false
const FILEIO_ENABLE = false
#Start Test Script

tic()

#ODE
println("Linear ODE Tests")
@time @test include("ode/ode_twodimlinear_tests.jl")
println("ODE Convergence Tests")
@time @test include("ode/ode_convergence_tests.jl")
println("ODE Adaptive Tests")
@time @test include("ode/ode_adaptive_tests.jl")
println("ODE Tstops Tests")
@time @test include("ode/ode_tstops_tests.jl")
println("ODE Unrolled Tests")
(LONGER_TESTS) && !is_windows() && @time @test include("ode/ode_unrolled_comparison_tests.jl")
println("ODE Initial Dt Tests")
@time @test include("ode/ode_initdt_tests.jl")
println("ODE Rosenbrock Tests")
@time @test include("ode/ode_rosenbrock_tests.jl")
println("ODE Initial Dt Tests")
!is_windows() && @time @test include("ode/ode_dense_tests.jl") # Windows 32-bit Overflow
println("ODE In-Place Tests")
@time @test include("ode/ode_inplace_tests.jl")
println("ODE Events Tests")
@time @test include("ode/ode_event_tests.jl")
println("ODE Cache Tests")
@time @test include("ode/ode_cache_tests.jl")
println("ODE saveat Tests")
@time @test include("ode/ode_saveat_tests.jl")
println("ODE Feagin Tests")
(LONGER_TESTS) && @time @test include("ode/ode_feagin_tests.jl")
println("ODE Number Type Tests")
@time @test include("ode/ode_numbertype_tests.jl")
println("ODE Ndim Complex Tests")
@time @test include("ode/ode_ndim_complex_tests.jl")

#Optional Items
println("Units Tests")
(LONGER_TESTS) && !is_windows() && (@time @test include("units_tests.jl")) # Too long for AppVeyor
println("ODEInterface Tests")
(TEST_CONDITIONAL_DEPS) && !is_windows() && (@time @test include("ode/ODEInterface_tests.jl"))
println("ODE.jl Tests")
(TEST_CONDITIONAL_DEPS) && @time @test include("ode/ODEJL_tests.jl")

toc()
