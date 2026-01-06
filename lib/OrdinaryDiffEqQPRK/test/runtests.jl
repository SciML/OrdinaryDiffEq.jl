using SafeTestsets

@time @safetestset "Quadruple Precision Tests" include("ode_quadruple_precision_tests.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
