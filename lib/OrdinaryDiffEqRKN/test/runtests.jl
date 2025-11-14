using SafeTestsets

@time @safetestset "Nystrom Convergence Tests" include("nystrom_convergence_tests.jl")

# Only run QA tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end