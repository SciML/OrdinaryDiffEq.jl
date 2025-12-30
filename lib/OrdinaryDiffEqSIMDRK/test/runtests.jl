using SafeTestsets

@time @safetestset "Convergence Tests" begin
    include("convergence_tests.jl")
end
@time @safetestset "Adaptivity Tests" begin
    include("adaptivity_tests.jl")
end

# Only run JET tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
end
