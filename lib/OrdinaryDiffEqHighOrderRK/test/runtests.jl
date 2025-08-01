using SafeTestsets

@time @safetestset "High Order ERK Convergence Tests" include("high_order_erk_convergence_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
# Only run allocation tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end