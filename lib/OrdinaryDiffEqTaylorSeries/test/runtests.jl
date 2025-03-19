using OrdinaryDiffEqTaylorSeries
using Test
using Symbolics
import OrdinaryDiffEqCore: alg_order, alg_stability_size, explicit_rk_docstring,
                           OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqMutableCache,
                           OrdinaryDiffEqConstantCache, @fold, trivial_limiter!,
                           constvalue, @unpack, perform_step!, calculate_residuals, @cache,
                           calculate_residuals!, _ode_interpolant, _ode_interpolant!,
                           CompiledFloats, @OnDemandTableauExtract, initialize!,
                           perform_step!, OrdinaryDiffEqAlgorithm,
                           DAEAlgorithm,
                           CompositeAlgorithm, _ode_addsteps!, copyat_or_push!,
                           AutoAlgSwitch, get_fsalfirstlast,
                           full_cache, DerivativeOrderNotPossibleError

using LinearAlgebra
using Plots
# using DifferentialEquations
# f(u, p, t) = cos(u)
# f!(du, u, p, t) = du .= cos.(u)

# u0 = 0.0
# u0! = [0.0]
# prob = ODEProblem{false, SciMLBase.NoSpecialize}(f, u0, (0.0, 10.0))
# prob! = ODEProblem{true, SciMLBase.NoSpecialize}(f!, u0!, (0.0, 10.0))
# sol = solve(prob, ExplicitTaylor2(), dt=0.01)
# sol! = solve(prob!, ExplicitTaylor2(), dt=0.01)
# # sol = solve(prob, DAETS(), dt=0.01)
# sol = solve(prob, ExplicitTaylor(order=Val(2)), dt=0.01)
# sol! = solve(prob!, ExplicitTaylor(order=Val(2)), dt=0.01)

# println("DONE with ODE tests")
include(joinpath(@__DIR__, "../src/DAETS_symbolics.jl"))
include(joinpath(@__DIR__, "../src/TaylorSeries_caches.jl"))
# println("Starting tests on DAETS")

# # --- Test Cases --- #
# @testset "Signature Matrix & Derivative Order Tests" begin
#     @syms t x(t) y(t)
#     # # Test max_derivative_order
#     # The variable itself.
#     order1 = max_derivative_order(x(t), x, t)
#     @test order1 == 0

#     # A multiplication expression: 7*x(t)
#     order2 = max_derivative_order(7*x(t), x, t)
#     @test order2 == 0

#     # First derivative: x'(t)
#     dx = Differential(t)(x(t))
#     # println(typeof(dx))
#     order3 = max_derivative_order(dx, x, t)
#     @test order3 == 1

#     # Second derivative: x''(t)
#     d2x = Differential(t)(dx)
#     order4 = max_derivative_order(d2x, x, t)
#     @test order4 == 2

#     # Expression that does not contain x(t): y(t)
#     order5 = max_derivative_order(y(t), x, t)
#     @test order5 == -Inf

#     # Test signature_matrix:
#     # Equation 1: f₁ = x'(t) + sin(x(t)) = 0
#     # Equation 2: f₂ = y''(t) - x(t) = 0
#     eq1 = Differential(t)(x(t)) + sin(x(t))
#     eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
#     eqs = [eq1, eq2]
#     # pass functions x and y (which take input t)
#     vars = [x, y]

#     Σ = signature_matrix(eqs, vars, t)
#     # First equation:
#     # Column 1: x(t) appears as x'(t) ⟹ order 1.
#     # Column 2: no y(t) ⟹ -Inf.
#     @test Σ[1, 1] == 1
#     @test Σ[1, 2] == -Inf

#     # Second equation:
#     # Column 1: x(t) appears as x(t) ⟹ order 0.
#     # Column 2: y(t) appears as y''(t) ⟹ order 2.
#     @test Σ[2, 1] == 0
#     @test Σ[2, 2] == 2
# end

# println("DONE Signature Matrix Tests")

# @testset "Highest Value Transversal Tests" begin
#     @syms t x(t) y(t)

#     # Same Equation
#     eq1 = Differential(t)(x(t)) + sin(x(t))
#     eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
#     eqs = [eq1, eq2]
#     vars = [x, y]
#     Σ = signature_matrix(eqs, vars, t)

#     # Expected signature matrix:
#     # [ 1  -Inf
#     #   0     2 ]
#     @test Σ == [1 -Inf; 0 2]
#     # Test HVT
#     transversal, value = highest_value_transversal(Σ)

#     # Expected transversal: [(1, 1), (2, 2)] (x'(t) in eq1 and y''(t) in eq2)
#     # Expected value: 1 + 2 = 3 TODO: Should probably add a check in case the value is infinite. In this case we say that the system is "ill-posed".
#     @test transversal == [(1, 1), (2, 2)]
#     @test value == 3
# end
# @testset "Highest Value Transversal Tests for 3x3 System" begin
#     @syms t x(t) y(t) z(t)

#     # Equation 1: f₁ = x''(t) + y(t) = 0
#     # Equation 2: f₂ = y'(t) + z(t) = 0
#     # Equation 3: f₃ = z''(t) + x(t) = 0
#     eq1 = Differential(t)(Differential(t)(x(t))) + y(t)
#     eq2 = Differential(t)(y(t)) + z(t)
#     eq3 = Differential(t)(Differential(t)(z(t))) + x(t)
#     eqs = [eq1, eq2, eq3]
#     vars = [x, y, z]
#     Σ = signature_matrix(eqs, vars, t)

#     # Expected signature matrix:
#     # [ 2   0  -Inf
#     #  -Inf  1     0
#     #   0  -Inf    2 ]
#     @test Σ == [2 0 -Inf; -Inf 1 0; 0 -Inf 2]
#     # Test HVT
#     transversal, value = highest_value_transversal(Σ)

#     # Expected transversal: [(1, 1), (2, 2), (3, 3)] (x''(t) in eq1, y'(t) in eq2, z''(t) in eq3)
#     # Expected value: 2 + 1 + 2 = 5
#     @test transversal == [(1, 1), (2, 2), (3, 3)]
#     @test value == 5
# end
# println("DONE Highest Value Transversal Tests")

# @testset "Find Offsets Tests for 2x2 System" begin
#     @syms t x(t) y(t)

#     # Same Equation as 2x2 System
#     eq1 = Differential(t)(x(t)) + sin(x(t))
#     eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
#     eqs = [eq1, eq2]
#     vars = [x, y]
#     Σ = signature_matrix(eqs, vars, t)

#     # Find the highest value transversal
#     transversal, value = highest_value_transversal(Σ)

#     # Find the offsets
#     c, d = find_offsets(Σ, transversal)

#     # Expected offsets (canonical):
#     # c = [0, 0]
#     # d = [1, 2]
#     @test c == [0, 0]
#     @test d == [1, 2]
# end
# @testset "Find Offsets Tests for 3x3 System" begin
#     @syms t x(t) y(t) z(t)

#     # Same Equation as 3x3 System
#     eq1 = Differential(t)(Differential(t)(x(t))) + y(t)
#     eq2 = Differential(t)(y(t)) + z(t)
#     eq3 = Differential(t)(Differential(t)(z(t))) + x(t)
#     eqs = [eq1, eq2, eq3]
#     vars = [x, y, z]
#     Σ = signature_matrix(eqs, vars, t)
#     transversal, value = highest_value_transversal(Σ)

#     # Test Offsets
#     c, d = find_offsets(Σ, transversal)

#     # Expected offsets (canonical):
#     # c = [0, 0, 0]
#     # d = [2, 1, 2]
#     @test c == [0, 0, 0]
#     @test d == [2, 1, 2]
# end

# println("DONE Find Offsets Tests")

# @testset "System Jacobian Tests for 2x2 System" begin
#     @syms t x(t) y(t)

#     # Same 2x2
#     eq1 = Differential(t)(x(t)) + sin(x(t))
#     eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
#     eqs = [eq1, eq2]
#     vars = [x, y]
#     Σ = signature_matrix(eqs, vars, t)
#     transversal, value = highest_value_transversal(Σ)
#     c, d = find_offsets(Σ, transversal)

#     # Convert c and d to Vector{Int}
#     c = Int.(c)
#     d = Int.(d)

#     # Test Jacobian
#     J = system_jacobian(eqs, vars, t, c, d, Σ)

#     # Expected Jacobian:
#     # [ 1   0
#     #   0   1 ]
#     @test isequal(J[1, 1], 1)
#     @test isequal(J[1, 2], 0)
#     @test isequal(J[2, 1], 0)
#     @test isequal(J[2, 2], 1)
# end
# @testset "System Jacobian Tests for 3x3 System" begin
#     @syms t x(t) y(t) z(t)

#     # Same 3x3
#     eq1 = Differential(t)(Differential(t)(x(t))) + y(t)
#     eq2 = Differential(t)(y(t)) + z(t)
#     eq3 = Differential(t)(Differential(t)(z(t))) + x(t)
#     eqs = [eq1, eq2, eq3]
#     vars = [x, y, z]
#     Σ = signature_matrix(eqs, vars, t)
#     transversal, value = highest_value_transversal(Σ)
#     c, d = find_offsets(Σ, transversal)

#     # Convert c and d to Vector{Int}
#     c = Int.(c)
#     d = Int.(d)

#     # Test Jacobian
#     J = system_jacobian(eqs, vars, t, c, d, Σ)

#     # Expected Jacobian:
#     # [1   0   0
#     #  0   1   1
#     #  0   0   1]
#     @test isequal(J[1, 1], 1)
#     @test isequal(J[1, 2], 0)
#     @test isequal(J[1, 3], 0)
#     @test isequal(J[2, 1], 0)
#     @test isequal(J[2, 2], 1)
#     @test isequal(J[2, 3], 0)
#     @test isequal(J[3, 1], 0)
#     @test isequal(J[3, 2], 0)
#     @test isequal(J[3, 3], 1)
# end
# println("DONE System Jacobian Tests")

# @testset "System Jacobian Tests for Simple Pendulum" begin
#     @syms t x(t) y(t) λ(t) G L

#     # Equations
#     f = Differential(t)(Differential(t)(x(t))) + x(t) * λ(t)
#     g = Differential(t)(Differential(t)(y(t))) + y(t) * λ(t) - G
#     h = x(t)^2 + y(t)^2 - L^2
#     eqs = [f, g, h]
#     vars = [x, y, λ]

#     # Construct the signature matrix
#     Σ = signature_matrix(eqs, vars, t)

#     # Find the highest value transversal
#     transversal, value = highest_value_transversal(Σ)

#     # Find the offsets
#     c, d = find_offsets(Σ, transversal)

#     # Convert c and d to Vector{Int}
#     c = Int.(c)
#     d = Int.(d)

#     # Construct the system Jacobian
#     J = system_jacobian(eqs, vars, t, c, d, Σ)

#     # Expected Jacobian:
#     # [1   0   x(t)
#     #  0   1   y(t)
#     #  2   2   0]
#     @test isequal(J[1, 1], 1)
#     @test isequal(J[1, 2], 0)
#     @test isequal(J[1, 3], x(t))
#     @test isequal(J[2, 1], 0)
#     @test isequal(J[2, 2], 1)
#     @test isequal(J[2, 3], y(t))
#     @test isequal(J[3, 1], 2x(t))
#     @test isequal(J[3, 2], 2y(t))
#     @test isequal(J[3, 3], 0)
# # end

@testset "Basic DAE Test" begin
    # Simple index-1 DAE:
    # x'(t) = y(t)
    # 0 = x(t) - sin(t)
    # 
    # Analytical solution:
    # x(t) = sin(t)
    # y(t) = cos(t)
    
    function dae_simple!(res, du, u, p, t)
        x, y = u
        dx, dy = du
        
        # Diffeq
        res[1] = dx - y
        # Alg eq
        res[2] = x - sin(t)
        
        return nothing
    end
    
    # Initial conditions
    u0 = [0.0, 1.0]    # [x(0), y(0)] = [sin(0), cos(0)]
    du0 = [1.0, 0.0]   # [x'(0), y'(0)] = [cos(0), -sin(0)]
    tspan = (0.0, 2π)
    prob = DAEProblem(dae_simple!, du0, u0, tspan)
    sol = solve(prob, DAETS(), dt=0.1)
    
    # Generate analytical solution
    t_analytical = 0:0.01:2π
    x_analytical = sin.(t_analytical)
    y_analytical = cos.(t_analytical)
    # numeric solutions
    t_numerical = sol.t
    x_numerical = [sol[1, i] for i in 1:length(sol.t)]
    y_numerical = [sol[2, i] for i in 1:length(sol.t)]
    
    # Plot
    p1 = plot(t_analytical, x_analytical, label="x analytici", title="DAE")
    plot!(p1, t_numerical, x_numerical, label="x numerical", marker=:circle, markersize=3)
    
    p2 = plot(t_analytical, y_analytical, label="y analytic", title="")
    plot!(p2, t_numerical, y_numerical, label="y numerical", marker=:circle, markersize=3)
    
    plot(p1, p2, layout=(2,1), size=(800, 600))
    savefig("dae_test_1_comparison.png")
    
    # test that solution matches analytical solution
    for i in eachindex(t_numerical)
        t = t_numerical[i]
        @test isapprox(x_numerical[i], sin(t), rtol=1e-2)
        @test isapprox(y_numerical[i], cos(t), rtol=1e-2)
    end
end

