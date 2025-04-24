using SafeTestsets

@time @safetestset include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")