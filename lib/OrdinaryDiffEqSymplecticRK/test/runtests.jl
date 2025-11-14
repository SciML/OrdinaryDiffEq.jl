using SafeTestsets

@time @safetestset "Synplectic Convergence Tests" include("symplectic_convergence.jl")
@time @safetestset "Synplectic Tests" include("symplectic_tests.jl")

# Only run QA tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end