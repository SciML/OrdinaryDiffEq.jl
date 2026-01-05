using SafeTestsets

@time @safetestset "Nystrom Convergence Tests" include("nystrom_convergence_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
