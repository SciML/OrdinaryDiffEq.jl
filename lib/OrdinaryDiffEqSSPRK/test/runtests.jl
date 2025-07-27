using SafeTestsets

@time @safetestset "SSPRK Tests" include("ode_ssprk_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
