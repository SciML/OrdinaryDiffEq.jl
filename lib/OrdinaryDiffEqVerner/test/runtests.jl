using SafeTestsets

# Only run QA and allocation tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
    @time @safetestset "RKV76IIa Tests" include("ode_verner_tests.jl")
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end
