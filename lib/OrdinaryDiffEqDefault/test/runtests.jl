using Pkg
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_gpu_env()
    Pkg.activate(joinpath(@__DIR__, "gpu"))
    return Pkg.instantiate()
end

function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

# Run GPU tests
if TEST_GROUP == "GPU"
    activate_gpu_env()
    @time @safetestset "Autoswitch GPU" include("gpu/autoswitch.jl")
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "Default Solver Tests" include("default_solver_tests.jl")
end

# Run QA tests (JET, Aqua)
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end
