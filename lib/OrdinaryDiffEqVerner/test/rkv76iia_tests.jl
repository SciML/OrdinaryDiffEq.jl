using OrdinaryDiffEqVerner, OrdinaryDiffEqCore, DiffEqBase, Test
using LinearAlgebra

# Simple convergence test for RKV76IIa
@testset "RKV76IIa Convergence Tests" begin
    # Test problem: y' = -y, y(0) = 1
    # Exact solution: y(t) = exp(-t)
    function f!(du, u, p, t)
        du[1] = -u[1]
    end
    
    function f(u, p, t)
        -u
    end
    
    # Out-of-place test
    prob_oop = ODEProblem(f, 1.0, (0.0, 1.0))
    sol_oop = solve(prob_oop, RKV76IIa(), abstol=1e-12, reltol=1e-12)
    @test sol_oop[end] ≈ exp(-1.0) atol=1e-10
    
    # In-place test
    prob_ip = ODEProblem(f!, [1.0], (0.0, 1.0))
    sol_ip = solve(prob_ip, RKV76IIa(), abstol=1e-12, reltol=1e-12)
    @test sol_ip[end][1] ≈ exp(-1.0) atol=1e-10
    
    # Test that it's order 7
    dts = (1/2) .^ (6:10)
    errors = zeros(length(dts))
    
    for (i, dt) in enumerate(dts)
        sol = solve(prob_oop, RKV76IIa(), dt=dt, adaptive=false)
        errors[i] = abs(sol[end] - exp(-1.0))
    end
    
    # Check convergence order
    for i in 2:length(errors)
        order = log(errors[i-1]/errors[i]) / log(2)
        @test order ≈ 7 atol=0.3  # Allow some tolerance
    end
end

# Test with a system of ODEs
@testset "RKV76IIa System Test" begin
    # Lorenz system
    function lorenz!(du, u, p, t)
        σ, ρ, β = p
        du[1] = σ * (u[2] - u[1])
        du[2] = u[1] * (ρ - u[3]) - u[2]
        du[3] = u[1] * u[2] - β * u[3]
    end
    
    u0 = [1.0, 0.0, 0.0]
    p = [10.0, 28.0, 8/3]
    tspan = (0.0, 0.1)
    prob = ODEProblem(lorenz!, u0, tspan, p)
    
    sol = solve(prob, RKV76IIa(), abstol=1e-10, reltol=1e-10)
    @test length(sol.t) > 2  # Should have taken multiple steps
    @test all(isfinite.(sol[end]))  # Solution should be finite
end

# Test lazy vs non-lazy
@testset "RKV76IIa Lazy Tests" begin
    prob = ODEProblem((u,p,t) -> -u, 1.0, (0.0, 1.0))
    
    sol_lazy = solve(prob, RKV76IIa(lazy=true), abstol=1e-10, reltol=1e-10)
    sol_not_lazy = solve(prob, RKV76IIa(lazy=false), abstol=1e-10, reltol=1e-10)
    
    # Both should give the same result
    @test sol_lazy[end] ≈ sol_not_lazy[end] atol=1e-12
end