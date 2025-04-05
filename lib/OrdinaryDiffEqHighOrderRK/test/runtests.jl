using SafeTestsets

@time @safetestset "High Order ERK Convergence Tests" include("high_order_erk_convergence_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")