using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run QA tests (JET, Aqua)
if TEST_GROUP != "FUNCTIONAL" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end

# Run functional tests
if TEST_GROUP != "QA"
    @time @safetestset "OOP J_t Tracking" include("oop_jt_tracking_test.jl")
end
