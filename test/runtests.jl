using OrdinaryDiffEq
using Base.Test

const CPU_FLOPS = peakflops()
const TEST_PLOT = false
const LONGER_TESTS = true #Requires JLD
const TEST_CONDITIONAL_DEPS = true
const FILEIO_ENABLE = false
#Start Test Script
using OrdinaryDiffEq, Compat
using Base.Test

tic()

#ODE
println("Linear ODE Tests")
@time @test include("ode/ode_twodimlinear_tests.jl")
println("ODE Convergence Tests")
@time @test include("ode/ode_convergence_tests.jl")
println("ODE Tableau Convergence Tests")
@compat !is_windows() && @time @test include("ode/ode_tableau_convergence_tests.jl") ## Windows 32-bit fails on Butcher62 convergence test
println("ODE Adaptive Tests")
@time @test include("ode/ode_adaptive_tests.jl")
println("ODE Tspan Tests")
@time @test include("ode/ode_tspan_tests.jl")
println("ODE Unrolled Tests")
(LONGER_TESTS) && @compat !is_windows() && @time @test include("ode/ode_unrolled_comparison_tests.jl")
println("ODE Initial Dt Tests")
@time @test include("ode/ode_initdt_tests.jl")
println("ODE AppxTrue Tests")
@time @test include("ode/ode_appxtrue_tests.jl")
println("ODE Rosenbrock Tests")
@time @test include("ode/ode_rosenbrock_tests.jl")
println("ODE Initial Dt Tests")
@compat !is_windows() && @time @test include("ode/ode_dense_tests.jl") # Windows 32-bit Overflow
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

#Optional Items
println("Units Tests")
(LONGER_TESTS) && @compat !is_windows() && (@time @test include("units_tests.jl")) # Too long for AppVeyor
println("ODEInterface Tests")
(TEST_CONDITIONAL_DEPS) && @compat !is_windows() && (@time @test include("ode/ODEInterface_tests.jl"))
println("ODE.jl Tests")
(TEST_CONDITIONAL_DEPS) && @time @test include("ode/ODEJL_tests.jl")
println("Sundials.jl Tests")
(TEST_CONDITIONAL_DEPS) && @time @test include("ode/Sundials_tests.jl") # Not until tags

toc()
