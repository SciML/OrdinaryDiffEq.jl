using SafeTestsets

@time @safetestset "Low Storage RK Tests" include("ode_low_storage_rk_tests.jl")

# Only run QA tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end