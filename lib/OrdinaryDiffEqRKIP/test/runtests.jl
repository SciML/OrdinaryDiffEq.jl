using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @safetestset "Type Safety Tests" include("type_test.jl")
    @safetestset "Cache Test" include("cache_recycling_test.jl")
    @safetestset "Fourier Semilinear PDE Tests" include("semilinear_pde_test_cpu.jl")
    # GPU tests moved to test/gpu/rkip_semilinear_pde.jl for Buildkite GPU CI
end

# Run QA tests (JET)
if TEST_GROUP != "FUNCTIONAL"
    @time @safetestset "JET Tests" include("jet.jl")
end
