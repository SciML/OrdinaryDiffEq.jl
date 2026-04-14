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
    @time @safetestset "Linear LSRK GPU" include("gpu/linear_lsrk.jl")
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "Low Storage RK Tests" include("ode_low_storage_rk_tests.jl")
    @time @safetestset "RDPK3 init regression (#3384)" include("rdpk3_init_tests.jl")
end

# Run QA tests (JET, Aqua)
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
