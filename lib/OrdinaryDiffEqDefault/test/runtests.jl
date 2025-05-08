using SafeTestsets
@time @safetestset "Default Solver Tests" include("default_solver_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")