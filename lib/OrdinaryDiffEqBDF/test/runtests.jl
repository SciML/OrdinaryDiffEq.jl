using SafeTestsets

@time @safetestset "BDF Convergence Tests" include("bdf_convergence_tests.jl")
@time @safetestset "BDF Regression Tests" include("bdf_regression_tests.jl")

@time @safetestset "DAE Convergence Tests" include("dae_convergence_tests.jl")
@time @safetestset "DAE AD Tests" include("dae_ad_tests.jl")
@time @safetestset "DAE Event Tests" include("dae_event.jl")
@time @safetestset "DAE Initialization Tests" include("dae_initialization_tests.jl")
