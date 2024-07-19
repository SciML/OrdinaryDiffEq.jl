using SafeTestsets

@time @safetestset "SSPRK Tests" include("ode_ssprk_tests.jl")
