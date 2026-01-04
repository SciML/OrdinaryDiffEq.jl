using SafeTestsets

@time @safetestset "DAE Convergence Tests" include("dae_convergence_tests.jl")
@time @safetestset "DAE AD Tests" include("dae_ad_tests.jl")
@time @safetestset "DAE Event Tests" include("dae_event.jl")
@time @safetestset "DAE Initialization Tests" include("dae_initialization_tests.jl")

@time @safetestset "BDF Inference Tests" include("inference_tests.jl")
@time @safetestset "BDF Convergence Tests" include("bdf_convergence_tests.jl")
@time @safetestset "BDF Regression Tests" include("bdf_regression_tests.jl")

# Only run QA and allocation tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end
