using Test
using LinearAlgebra: norm

include("semilinear_fft.jl")

@testset "Lugiato-Lefever equation test - Float64" begin
    @test !isnothing(lle_scan_test())
end

@testset "Lugiato-Lefever equation test - Float32" begin
    @test !isnothing(lle_scan_test(Float32))
end

@testset "NLSE test - Float64" begin
    sol = nlse_test()

    # NLSE Soliton solution, only the phase change
    p1 = abs2.(sol.u[1])
    p2 = abs2.(sol.u[end])

    @test norm(p1 .- p2) / norm(p1) < 1.0e-2 # small difference with analytical solution due to boundary periodic condition
end

@testset "NLSE test - Float32" begin
    sol = nlse_test(Float32)

    # NLSE Soliton solution, only the phase change
    p1 = abs2.(sol.u[1])
    p2 = abs2.(sol.u[end])

    @test norm(p1 .- p2) / norm(p1) < 1.0e-2 # small difference with analytical solution due to boundary periodic condition
end
