using SafeTestsets

@safetestset "Type Safety Tests" include("type_test.jl")
@safetestset "Cache Test" include("cache_recycling_test.jl")
@safetestset "Fourier Semilinear PDE Tests" include("semilinear_pde_test_cpu.jl")
@safetestset "CUDA Test" include("semilinear_pde_test_cuda.jl")
