using SafeTestsets

@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")

# Only run allocation and type stability tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
    @time @safetestset "JET Type Stability Tests" include("jet_tests.jl")
end