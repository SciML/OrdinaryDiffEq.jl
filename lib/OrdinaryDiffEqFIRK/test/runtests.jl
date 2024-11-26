using SafeTestsets

@time @safetestset "FIRK Tests" include("ode_firk_tests.jl")
@time @safetestset "High Order FIRK Tests" include("ode_high_order_firk_tests.jl")
