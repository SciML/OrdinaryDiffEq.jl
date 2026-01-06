using SafeTestsets

@time @safetestset "Feagin Tests" include("ode_feagin_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
