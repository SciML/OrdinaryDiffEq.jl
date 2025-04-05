using SafeTestsets

@time @safetestset "Newton Tests" include("newton_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")