using SafeTestsets

@time @safetestset "Synplectic Convergence Tests" include("symplectic_convergence.jl")
@time @safetestset "Synplectic Tests" include("symplectic_tests.jl")
