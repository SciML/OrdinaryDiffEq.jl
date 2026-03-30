using OrdinaryDiffEqImplicitTableaus
using DiffEqBase
using Test

@testset "OrdinaryDiffEqImplicitTableaus" begin
    @testset "Basic tableau construction" begin
        tab = constructImplicitEuler()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 1

        tab = constructGL2()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 2

        tab = constructGL4()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 4

        tab = constructGL6()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 6

        tab = constructRadauIIA5()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 5

        tab = constructLobattoIIIA4()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 4

        tab = constructTrapezoidalRule()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 2

        tab = constructMidpointRule()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 2
    end

    @testset "Two-type construction" begin
        # Test that T_time controls c type and T controls A, α types
        tab = constructGL2(BigFloat, Float64)
        @test eltype(tab.c) == Float64
        @test eltype(tab.A) == BigFloat
        @test eltype(tab.α) == BigFloat

        tab = constructRadauIIA5(BigFloat, Float64)
        @test eltype(tab.c) == Float64
        @test eltype(tab.A) == BigFloat
        @test eltype(tab.α) == BigFloat

        tab = constructLobattoIIIC4(Float64, Float32)
        @test eltype(tab.c) == Float32
        @test eltype(tab.A) == Float64

        # Single type argument should use same type for both
        tab = constructImplicitEuler(Float32)
        @test eltype(tab.c) == Float32
        @test eltype(tab.A) == Float32
    end

    @testset "Default construction" begin
        # Default should use Float64 for everything
        tab = constructGL4()
        @test eltype(tab.c) == Float64
        @test eltype(tab.A) == Float64
        @test eltype(tab.α) == Float64
    end
end
