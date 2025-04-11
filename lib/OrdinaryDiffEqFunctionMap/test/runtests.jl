using SafeTestsets

@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa_tests.jl")