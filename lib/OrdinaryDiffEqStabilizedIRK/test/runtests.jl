using SafeTestsets

@time @safetestset "IRKC Tests" include("irkc_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")