using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @time @safetestset "Extrapolation Tests" include("ode_extrapolation_tests.jl")
end

if TEST_GROUP == "Multithreading"
    @time @safetestset "Multithreaded Extrapolation Tests" include("multithreading/ode_extrapolation_tests.jl")
end

# Run QA tests (JET, Aqua)
if TEST_GROUP != "Core" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
