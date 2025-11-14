using SafeTestsets


# Only run QA tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
@time @safetestset "DiscreteProblem Defaults" include("discrete_problem_defaults.jl")