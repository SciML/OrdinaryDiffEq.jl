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
    @time @safetestset "RKIP Semilinear PDE GPU" include("gpu/rkip_semilinear_pde.jl")
end

# Run functional tests
if TEST_GROUP != "QA" && TEST_GROUP != "GPU"
    @safetestset "Type Safety Tests" include("type_test.jl")
    @safetestset "Cache Test" include("cache_recycling_test.jl")
    @safetestset "Fourier Semilinear PDE Tests" include("semilinear_pde_test_cpu.jl")
end

# Run QA tests (JET)
if TEST_GROUP != "Core" && TEST_GROUP != "GPU" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
end
