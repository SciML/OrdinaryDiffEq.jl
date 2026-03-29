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
    @time @safetestset "Simple GPU" include("gpu/simple_gpu.jl")
    @time @safetestset "Hermite Interpolation GPU" include("gpu/hermite_test.jl")
end

# Run QA tests (JET, Aqua)
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end

# Functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "Sparse isdiag Performance" include("sparse_isdiag_tests.jl")
    @time @safetestset "DynamicQuantities + Measurements" include("dynamicquantities_measurements.jl")
end
