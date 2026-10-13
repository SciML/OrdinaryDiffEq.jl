using StochasticDiffEqHighOrder
import StochasticDiffEqHighOrder.SciMLBase
using LinearAlgebra
using Test

# Non-diagonal (matrix) additive noise: OOP SRA1 must use matrix-vector
# products (gpdt * chi2), matching the IIP is_diagonal_noise branch.

@testset "SRA1 OOP non-diagonal additive noise" begin
    B = [0.1 0.05; 0.05 0.1; 0.02 0.03; 0.1 0.0]
    A = [-1.0 0.1 0.0 0.0; 0.1 -1.0 0.1 0.0; 0.0 0.1 -1.0 0.1; 0.0 0.0 0.1 -1.0]
    u0 = [1.0, 0.5, 0.25, 0.125]
    tspan = (0.0, 1.0)
    nrp = zeros(4, 2)
    p = (A = A, B = B)

    f_oop(u, p, t) = p.A * u
    g_oop(u, p, t) = p.B
    f_iip(du, u, p, t) = mul!(du, p.A, u)
    g_iip(du, u, p, t) = (du .= p.B; nothing)

    prob_oop = SDEProblem(f_oop, g_oop, u0, tspan, p; noise_rate_prototype = nrp)
    prob_iip = SDEProblem(f_iip, g_iip, copy(u0), tspan, p; noise_rate_prototype = copy(nrp))

    sol_iip = solve(prob_iip, SRA1(); seed = 1, adaptive = false, dt = 1 / 2^6)
    sol_oop = solve(prob_oop, SRA1(); seed = 1, adaptive = false, dt = 1 / 2^6)
    @test SciMLBase.successful_retcode(sol_oop)
    @test SciMLBase.successful_retcode(sol_iip)
    @test sol_oop.u ≈ sol_iip.u
    @test all(isfinite, sol_oop.u[end])

    sol_iip_adapt = solve(prob_iip, SRA1(); seed = 2)
    sol_oop_adapt = solve(prob_oop, SRA1(); seed = 2)
    @test SciMLBase.successful_retcode(sol_oop_adapt)
    @test sol_oop_adapt.u[end] ≈ sol_iip_adapt.u[end] rtol = 1.0e-10
end

@testset "SRA1 OOP diagonal additive noise still solves" begin
    f_oop(u, p, t) = -u
    g_oop(u, p, t) = [0.1, 0.2, 0.3, 0.4]
    u0 = [1.0, 0.5, 0.25, 0.125]
    prob_oop = SDEProblem(f_oop, g_oop, u0, (0.0, 1.0))
    sol_oop = solve(prob_oop, SRA1(); seed = 3, adaptive = false, dt = 1 / 2^6)
    @test SciMLBase.successful_retcode(sol_oop)
    @test all(isfinite, sol_oop.u[end])
end
