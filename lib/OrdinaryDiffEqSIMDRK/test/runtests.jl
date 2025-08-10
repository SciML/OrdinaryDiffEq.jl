using SafeTestsets

@time @safetestset "Convergence Tests" begin
    include("convergence_tests.jl")
end
@time @safetestset "Adaptivity Tests" begin
    include("adaptivity_tests.jl")
end
