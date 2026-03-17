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
    @time @safetestset "Simple DAE GPU" include("gpu/simple_dae.jl")
end

# Run functional tests
if TEST_GROUP == "Core" || TEST_GROUP == "ALL"
    @time @safetestset "DAE Rosenbrock AD Tests" include("dae_rosenbrock_ad_tests.jl")
    @time @safetestset "Rosenbrock Convergence Tests" include("ode_rosenbrock_tests.jl")
    @time @safetestset "IMEXRKR Tests" include("imexrkr_tests.jl")
end

# Run QA tests (JET, Aqua, AllocCheck)
if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end
