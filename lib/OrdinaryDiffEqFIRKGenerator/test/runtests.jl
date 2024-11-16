using SafeTestsets

@time @safetestset "Generated FIRK Tests" include("ode_firk_tests.jl")
