using SafeTestsets

@time @safetestset "Rosenbrock Convergence Tests" include("ode_rosenbrock_tests.jl")
@time @safetestset "DAE Rosenbrock AD Tests" include("dae_rosenbrock_ad_tests.jl")
