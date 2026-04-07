using OrdinaryDiffEqExplicitTableaus
using DiffEqBase
using DiffEqDevTools
using Test

const ET = OrdinaryDiffEqExplicitTableaus

@testset "OrdinaryDiffEqExplicitTableaus" begin
    @testset "Basic tableau construction" begin
        tab = ET.DormandPrince()
        @test tab isa DiffEqBase.ExplicitRKTableau
        @test tab.order == 5

        tab = ET.BogakiShampine3()
        @test tab isa DiffEqBase.ExplicitRKTableau
        @test tab.order == 3

        tab = ET.RKF5()
        @test tab isa DiffEqBase.ExplicitRKTableau
        @test tab.order == 5

        tab = ET.Euler()
        @test tab isa DiffEqBase.ExplicitRKTableau
        @test tab.order == 1

        tab = ET.RK4()
        @test tab isa DiffEqBase.ExplicitRKTableau
        @test tab.order == 4
    end

    @testset "Two-type construction" begin
        tab = ET.BogakiShampine3(BigFloat, Float64)
        @test eltype(tab.c) == Float64
        @test eltype(tab.A) == BigFloat
        @test eltype(tab.α) == BigFloat

        tab = ET.DormandPrince(BigFloat, Float64)
        @test eltype(tab.c) == Float64
        @test eltype(tab.A) == BigFloat
        @test eltype(tab.α) == BigFloat

        tab = ET.RK4(Float64, Float32)
        @test eltype(tab.c) == Float32
        @test eltype(tab.A) == Float64

        tab = ET.CashKarp(Float32)
        @test eltype(tab.c) == Float32
        @test eltype(tab.A) == Float32
    end

    @testset "Tableau order conditions (check_tableau)" begin
        setprecision(400)
        for name in names(ET; all = false)
            fn = getproperty(ET, name)
            fn isa Function || continue
            tab = try
                fn(BigFloat)
            catch
                continue
            end
            tab isa DiffEqBase.ExplicitRKTableau || continue
            if tab.order < 12
                if name in (:Tsitouras9, :Tsitouras92)
                    @info "Known broken: $name"
                    @test_broken check_tableau(tab)
                else
                    @info "Testing $name..."
                    @test check_tableau(tab)
                end
            end
        end
    end

    @testset "High-order convergence tests" begin
        using OrdinaryDiffEq, Random
        using ODEProblemLibrary: prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear

        setprecision(400)
        Random.seed!(100)
        dts = 1 .// 2 .^ (8:-1:4)
        testTol = 0.3

        for bigprob in [prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear]
            for tab in [
                ET.Feagin12(BigFloat),
                ET.Feagin14(BigFloat),
                ET.Ono12(BigFloat)
            ]
                @info "Testing high-order convergence..."
                tabalg = ExplicitRK(tableau = tab)
                sim = test_convergence(dts, bigprob, tabalg)
                @test sim.𝒪est[:l∞] >= tab.order ||
                      abs(sim.𝒪est[:l∞] - tab.order) < testTol
            end
        end
    end

    @testset "Stability regions" begin
        @test stability_region(ET.DormandPrince6(), initial_guess = -3.5)≈-3.95413 rtol=1e-3
        @test stability_region(ET.TsitourasPapakostas6(), initial_guess = -3.5)≈-3.95413 rtol=1e-3

        @test imaginary_stability_interval(ET.SSPRK33())≈sqrt(3)
        @test imaginary_stability_interval(ET.SSPRK33(Float32))≈sqrt(3.0f0)
        @test imaginary_stability_interval(ET.Kutta3())≈sqrt(3)
        @test imaginary_stability_interval(ET.Kutta3(Float32))≈sqrt(3.0f0)
        @test imaginary_stability_interval(ET.RK4())≈2.8284271247
    end

    @testset "Deduce Butcher tableau" begin
        using OrdinaryDiffEq

        function coefficients_as_in_tableau(A, b, c, tab)
            if size(A) == size(tab.A)
                AA = A
                bb = b
                cc = c
            else
                AA = A[1:(end - 1), 1:(end - 1)]
                bb = b[1:(end - 1)]
                cc = c[1:(end - 1)]
            end
            AA, bb, cc
        end

        let erk = Heun()
            A, b, c = deduce_Butcher_tableau(erk)
            tab = ET.Heun()
            AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
            @test AA ≈ tab.A
            @test bb ≈ tab.α
            @test cc ≈ tab.c
        end

        let erk = OrdinaryDiffEq.Euler()
            A, b, c = deduce_Butcher_tableau(erk)
            tab = ET.Euler()
            AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
            @test AA ≈ tab.A
            @test bb ≈ tab.α
            @test cc ≈ tab.c
        end

        let erk = OrdinaryDiffEq.RK4()
            A, b, c = deduce_Butcher_tableau(erk)
            tab = ET.RK4()
            AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
            @test AA ≈ tab.A
            @test bb ≈ tab.α
            @test cc ≈ tab.c
        end

        let erk = Tsit5()
            A, b, c = deduce_Butcher_tableau(erk)
            tab = ET.Tsitouras5()
            AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
            @test AA ≈ tab.A
            @test bb ≈ tab.α
            @test cc ≈ tab.c
        end

        let erk = BS3()
            A, b, c = deduce_Butcher_tableau(erk)
            tab = ET.BogakiShampine3()
            AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
            @test AA ≈ tab.A
            @test bb ≈ tab.α
            @test cc ≈ tab.c
        end

        let erk = BS5()
            A, b, c = deduce_Butcher_tableau(erk)
            tab = ET.BogakiShampine5()
            AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
            @test AA ≈ tab.A
            @test bb ≈ tab.α
            @test cc ≈ tab.c
        end

        let erk = Vern7()
            A, b, c = deduce_Butcher_tableau(erk)
            tab = ET.Verner7()
            AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
            @test AA ≈ tab.A
            @test bb ≈ tab.α
            @test cc ≈ tab.c
        end

        let erk = Vern8()
            A, b, c = deduce_Butcher_tableau(erk)
            tab = ET.Verner8()
            AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
            @test AA ≈ tab.A
            @test bb ≈ tab.α
            @test cc ≈ tab.c
        end

        let erk = RKO65()
            A, b, c = deduce_Butcher_tableau(erk)
            tab = ET.RKO65()
            AA = A[2:end, 2:end]
            bb = b[2:end]
            cc = c[2:end]
            @test AA ≈ tab.A
            @test bb ≈ tab.α
            @test cc ≈ tab.c
        end
    end
end
