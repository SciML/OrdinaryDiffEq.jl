using LinearAlgebra: diag
using OrdinaryDiffEqSDIRK
using SciMLBase: successful_retcode
using StaticArrays: @SMatrix, @SVector
using Test

function stratreac_rates(u, t)
    O1D, O, O3, O2, NO, NO2 = u
    Tr = 4.5
    Ts = 19.5
    T = mod(t / 3600, 24)
    if Tr <= T <= Ts
        Tfrac = (2 * T - Tr - Ts) / (Ts - Tr)
        sigma = 0.5 + 0.5 * cos(pi * abs(Tfrac) * Tfrac)
    else
        sigma = zero(t)
    end

    M = 8.12e16
    k1 = 2.643e-10 * sigma^3
    k2 = 8.018e-17
    k3 = 6.12e-4 * sigma
    k4 = 1.567e-15
    k5 = 1.07e-3 * sigma^2
    k6 = 7.11e-11
    k7 = 1.2e-10
    k8 = 6.062e-15
    k9 = 1.069e-11
    k10 = 1.289e-2 * sigma
    k11 = 1.0e-8

    return (
        k1 * O2,
        k2 * O * O2,
        k3 * O3,
        k4 * O3 * O,
        k5 * O3,
        k6 * M * O1D,
        k7 * O1D * O3,
        k8 * O3 * NO,
        k9 * NO2 * O,
        k10 * NO2,
        k11 * NO * O,
    )
end

function stratreac_production(u, p, t)
    r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11 = stratreac_rates(u, t)
    return @SMatrix [
        0.0 0.0 r5 0.0 0.0 0.0;
        r6 r1 + r10 r3 r1 0.0 0.0;
        0.0 r2 0.0 0.0 0.0 0.0;
        r7 r4 + r9 r4 + r7 + r8 r3 + r5 0.0 0.0;
        0.0 0.0 0.0 0.0 0.0 r9 + r10;
        0.0 0.0 0.0 0.0 r8 + r11 0.0
    ]
end

function stratreac_destruction(u, p, t)
    _, r2, _, _, _, _, _, _, _, _, r11 = stratreac_rates(u, t)
    return @SVector [0.0, r11, 0.0, r2, 0.0, 0.0]
end

function stratreac_rhs(u, p, t)
    r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11 = stratreac_rates(u, t)
    return @SVector [
        r5 - r6 - r7;
        2 * r1 - r2 + r3 - r4 + r6 - r9 + r10 - r11;
        r2 - r3 - r4 - r5 - r7 - r8;
        -r1 - r2 + r3 + 2 * r4 + r5 + 2 * r7 + r8 + r9;
        -r8 + r9 + r10 - r11;
        r8 - r9 - r10 + r11
    ]
end

function stratreac_pds_rhs(u, p, t)
    P = stratreac_production(u, p, t)
    D = stratreac_destruction(u, p, t)
    return diag(P) + vec(sum(P; dims = 2)) - vec(sum(P; dims = 1)) - vec(D)
end

@testset "ImplicitEuler static-array production-destruction RHS" begin
    u0 = @SVector [9.906e1, 6.624e8, 5.326e11, 1.697e16, 4.0e6, 1.093e9]
    tspan = (4.32e4, 3.024e5)
    solve_kwargs = (; dt = 1.0, abstol = 1.0e-5, reltol = 1.0e-5)

    sol = solve(ODEProblem(stratreac_rhs, u0, tspan), ImplicitEuler(); solve_kwargs...)
    sol_pds = solve(
        ODEProblem(stratreac_pds_rhs, u0, tspan), ImplicitEuler();
        solve_kwargs...
    )

    @test successful_retcode(sol)
    @test successful_retcode(sol_pds)
    @test isapprox(sol(2.0e5), sol_pds(2.0e5))
end
