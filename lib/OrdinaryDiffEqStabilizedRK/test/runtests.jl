using SafeTestsets

@time @safetestset "RKC Tests" include("rkc_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
