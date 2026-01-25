using SafeTestsets

# Run interpolation tests on all Julia versions
@time @safetestset "Generic RK Interpolation Tests" include("interpolation_tests.jl")

# Only run QA and allocation tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end
