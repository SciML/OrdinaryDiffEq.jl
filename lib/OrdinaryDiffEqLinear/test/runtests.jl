using SafeTestsets

@time @safetestset "Linear Methods Tests" include("linear_method_tests.jl")

# Only run QA tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end