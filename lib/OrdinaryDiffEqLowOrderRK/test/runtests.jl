using SafeTestsets

@time @safetestset "Low Order ERK Convergence Tests" include("low_order_erk_convergence_tests.jl")
@time @safetestset "OwrenZen Tests" include("owrenzen_tests.jl")
@time @safetestset "Euler SSP Tests" include("euler_ssp.jl")
@time @safetestset "JET Tests" include("jet.jl")