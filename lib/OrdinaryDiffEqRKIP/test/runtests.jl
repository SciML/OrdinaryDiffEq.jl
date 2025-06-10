using SafeTestsets

@safetestset "Type Safety Tests" begin include("type_test.jl") end
@safetestset "Fourier Semilinear PDE Tests" begin include("semilinear_pde_test_cpu.jl") end
@safetestset "CUDA Test" begin include("semilinear_pde_test_cuda.jl") end