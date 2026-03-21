using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run QA tests (JET, Aqua)
if TEST_GROUP != "Core" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end

@time @safetestset "Convergance" include("sdirk_convergence_tests.jl")
@time @safetestset "DAE tests" include("dae_esdirk_test.jl")
