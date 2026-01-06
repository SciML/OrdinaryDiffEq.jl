using SafeTestsets

@time @safetestset "Nordsieck Tests" include("nordsieck_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
