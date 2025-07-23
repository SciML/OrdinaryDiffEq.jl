using SafeTestsets

@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
@time @safetestset "RKV76IIa Tests" include("rkv76iia_tests.jl")