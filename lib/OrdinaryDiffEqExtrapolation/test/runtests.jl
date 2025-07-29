using SafeTestsets

@time @safetestset "Extrapolation Tests" include("ode_extrapolation_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
