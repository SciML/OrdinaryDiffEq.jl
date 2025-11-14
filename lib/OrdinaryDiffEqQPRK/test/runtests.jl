using SafeTestsets

@time @safetestset "Quadruple Precision Tests" include("ode_quadruple_precision_tests.jl")

# Only run QA tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end