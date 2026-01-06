using SafeTestsets

@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa.jl")
@time @safetestset "DiscreteProblem Defaults" include("discrete_problem_defaults.jl")
