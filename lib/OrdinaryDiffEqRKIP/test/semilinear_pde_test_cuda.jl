using Test
using LinearAlgebra: norm
using CUDA

# Skip all CUDA tests if CUDA is not functional
if !CUDA.functional()
    @info "CUDA not functional, skipping CUDA tests"
    @testset "CUDA Tests (skipped - CUDA not functional)" begin
        @test_skip true
    end
else
    include("semilinear_fft.jl")

    CUDA.allowscalar(false) # ensure no scalar indexing is used (for performance)

    @testset "Lugiato-Lefever equation test on GPU with CUDA - Float64" begin
        @test !isnothing(lle_scan_test(Val(true)))
    end

    @testset "Lugiato-Lefever equation test on GPU with CUDA - Float32" begin
        @test !isnothing(lle_scan_test(Float32, Val(true)))
    end

    @testset "NLSE test on GPU with CUDA - Float64" begin
        sol = nlse_test(Val(true))

        # NLSE Soliton solution, only the phase change
        p1 = abs2.(sol.u[1])
        p2 = abs2.(sol.u[end])

        @test norm(p1 .- p2) / norm(p1) < 1.5e-2 # small difference with analytical sech solution due to boundary periodic condition and PDE discretization
    end

    @testset "NLSE test on GPU with CUDA - Float32" begin
        sol = nlse_test(Float32, Val(true))

        # NLSE Soliton solution, only the phase change
        p1 = abs2.(sol.u[1])
        p2 = abs2.(sol.u[end])

        @test norm(p1 .- p2) / norm(p1) < 1.5e-2 # small difference with analytical sech solution due to boundary periodic condition and PDE discretization
    end
end
