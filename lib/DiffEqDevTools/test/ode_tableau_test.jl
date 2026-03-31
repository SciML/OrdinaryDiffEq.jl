using Test
using LinearAlgebra
using OrdinaryDiffEq, DiffEqDevTools

# Since non-FSAL methods are FSAL'ed artificially in OrdinaryDiffEq.jl,
# we need to exclude this last FSAL stage for the comparison with
# tableaus.
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

function ode_tableau_tests(T)
    let erk = Heun()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructHeun()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = Euler()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructEuler()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = RK4()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructRK4()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = SSPRK22()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructSSPRK22()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = SSPRK33()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructSSPRK33()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = SSPRK432()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructSSPRK43()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = SSPRK104()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructSSPRK104()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = Tsit5()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructTsitouras5()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = BS3()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructBogakiShampine3()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = BS5()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructBogakiShampine5()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = Vern7()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructVerner7()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = Vern8()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructVerner8()
        AA, bb, cc = coefficients_as_in_tableau(A, b, c, tab)
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end

    let erk = RKO65()
        A, b, c = deduce_Butcher_tableau(erk)
        tab = constructRKO65()
        AA = A[2:end, 2:end]
        bb = b[2:end]
        cc = c[2:end]
        @test AA ≈ tab.A
        @test bb ≈ tab.α
        @test cc ≈ tab.c
    end
end

ode_tableau_tests(Float32)
ode_tableau_tests(Float64)
