using SafeTestsets

@time @safetestset "FIRK Tests" include("ode_firk_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
