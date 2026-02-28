using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @time @safetestset "Newton Tests" include("newton_tests.jl")
    @time @safetestset "Sparse Algebraic Detection" include("sparse_algebraic_detection_tests.jl")
    @time @safetestset "Sparse DAE Initialization" include("sparse_dae_initialization_tests.jl")
end

# Run QA tests (JET, Aqua)
if TEST_GROUP != "FUNCTIONAL" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
