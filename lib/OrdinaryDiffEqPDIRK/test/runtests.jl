using SafeTestsets


# Only run QA tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
end