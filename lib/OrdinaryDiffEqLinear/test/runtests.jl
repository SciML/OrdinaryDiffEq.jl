using SafeTestsets

@time @safetestset "Linear Methods Tests" include("linear_method_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
