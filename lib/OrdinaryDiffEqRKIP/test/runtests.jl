using SafeTestsets

@safetestset "Type Safety Tests" include("type_test.jl")
@safetestset "Cache Test" include("cache_recycling_test.jl")
@safetestset "Fourier Semilinear PDE Tests" include("semilinear_pde_test_cpu.jl")
# GPU tests moved to test/gpu/rkip_semilinear_pde.jl for Buildkite GPU CI

# Only run JET tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
end
