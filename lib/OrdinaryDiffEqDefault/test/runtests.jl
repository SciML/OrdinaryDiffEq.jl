using Pkg
using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_gpu_env()
    Pkg.activate(joinpath(@__DIR__, "gpu"))
    return Pkg.instantiate()
end

# Run GPU tests
if TEST_GROUP == "GPU"
    activate_gpu_env()
    @time @safetestset "Autoswitch GPU" include("gpu/autoswitch.jl")
end

# Run functional tests
if TEST_GROUP != "QA" && TEST_GROUP != "GPU"
    @time @safetestset "Default Solver Tests" include("default_solver_tests.jl")
end

# Run QA tests (JET, Aqua)
if TEST_GROUP != "Core" && TEST_GROUP != "GPU" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
