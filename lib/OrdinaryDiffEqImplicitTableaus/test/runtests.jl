using OrdinaryDiffEqImplicitTableaus
using DiffEqBase
using DiffEqDevTools
using Test

const IT = OrdinaryDiffEqImplicitTableaus

@testset "OrdinaryDiffEqImplicitTableaus" begin
    @testset "Basic tableau construction" begin
        tab = IT.ImplicitEuler()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 1

        tab = IT.GL2()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 2

        tab = IT.GL4()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 4

        tab = IT.GL6()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 6

        tab = IT.RadauIIA5()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 5

        tab = IT.LobattoIIIA4()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 4

        tab = IT.TrapezoidalRule()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 2

        tab = IT.MidpointRule()
        @test tab isa DiffEqBase.ImplicitRKTableau
        @test tab.order == 2
    end

    @testset "Two-type construction" begin
        tab = IT.GL2(BigFloat, Float64)
        @test eltype(tab.c) == Float64
        @test eltype(tab.A) == BigFloat
        @test eltype(tab.α) == BigFloat

        tab = IT.RadauIIA5(BigFloat, Float64)
        @test eltype(tab.c) == Float64
        @test eltype(tab.A) == BigFloat
        @test eltype(tab.α) == BigFloat

        tab = IT.LobattoIIIC4(Float64, Float32)
        @test eltype(tab.c) == Float32
        @test eltype(tab.A) == Float64

        tab = IT.ImplicitEuler(Float32)
        @test eltype(tab.c) == Float32
        @test eltype(tab.A) == Float32
    end

    @testset "Default construction" begin
        tab = IT.GL4()
        @test eltype(tab.c) == Float64
        @test eltype(tab.A) == Float64
        @test eltype(tab.α) == Float64
    end

    @testset "Tableau order conditions (check_tableau)" begin
        setprecision(400)
        for name in names(IT; all = false)
            fn = getproperty(IT, name)
            fn isa Function || continue
            tab = try
                fn(BigFloat)
            catch
                continue
            end
            tab isa DiffEqBase.ImplicitRKTableau || continue
            @info "Testing $name..."
            @test check_tableau(tab)
        end
    end

    @testset "Stability regions" begin
        @test stability_region(IT.RadauIIA5(), initial_guess = 12.0)≈11.84 rtol=1e-2
    end
end
