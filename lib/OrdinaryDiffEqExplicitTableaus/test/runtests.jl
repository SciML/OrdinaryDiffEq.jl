using OrdinaryDiffEqExplicitTableaus
using DiffEqBase
using Test

@testset "OrdinaryDiffEqExplicitTableaus" begin
    @testset "Basic tableau construction" begin
        tab = constructDormandPrince()
        @test tab isa DiffEqBase.ExplicitRKTableau
        @test tab.order == 5

        tab = constructBogakiShampine3()
        @test tab isa DiffEqBase.ExplicitRKTableau
        @test tab.order == 3

        tab = constructRKF5()
        @test tab isa DiffEqBase.ExplicitRKTableau
        @test tab.order == 5

        tab = constructEuler()
        @test tab isa DiffEqBase.ExplicitRKTableau
        @test tab.order == 1

        tab = constructRK4()
        @test tab isa DiffEqBase.ExplicitRKTableau
        @test tab.order == 4
    end

    @testset "Two-type construction" begin
        # Test that T_time controls c type and T controls A, α types
        tab = constructBogakiShampine3(BigFloat, Float64)
        @test eltype(tab.c) == Float64
        @test eltype(tab.A) == BigFloat
        @test eltype(tab.α) == BigFloat

        tab = constructDormandPrince(BigFloat, Float64)
        @test eltype(tab.c) == Float64
        @test eltype(tab.A) == BigFloat
        @test eltype(tab.α) == BigFloat

        tab = constructRK4(Float64, Float32)
        @test eltype(tab.c) == Float32
        @test eltype(tab.A) == Float64

        # Single type argument should use same type for both
        tab = constructCashKarp(Float32)
        @test eltype(tab.c) == Float32
        @test eltype(tab.A) == Float32
    end

    @testset "Tsit5 ExplicitRK tableau" begin
        tab = constructTsit5ExplicitRK(Float64)
        @test tab isa DiffEqBase.ExplicitRKTableau
        @test tab.order == 5

        B = construct_tsit5_interp_matrix(Float64)
        @test size(B) == (7, 5)
        @test eltype(B) == Float64

        B32 = construct_tsit5_interp_matrix(Float32)
        @test eltype(B32) == Float32
    end
end
