using OrdinaryDiffEqTaylorSeries, ODEProblemLibrary, DiffEqDevTools
using Test
using Symbolics

@testset "Taylor2 Convergence Tests" begin
    # Test convergence
    dts = 2. .^ (-8:-4)
    testTol = 0.2
    sim = test_convergence(dts, prob_ode_linear, ExplicitTaylor2())
    @test sim.𝒪est[:final]≈2 atol=testTol
    sim = test_convergence(dts, prob_ode_2Dlinear, ExplicitTaylor2())
    @test sim.𝒪est[:final]≈2 atol=testTol
end

@testset "TaylorN Convergence Tests" begin
    # Test convergence
    dts = 2. .^ (-8:-4)
    testTol = 0.2
    for N in 3:4
        alg = ExplicitTaylor(order=Val(N))
        sim = test_convergence(dts, prob_ode_linear, alg)
        @test sim.𝒪est[:final]≈N atol=testTol
        sim = test_convergence(dts, prob_ode_2Dlinear, alg)
        @test sim.𝒪est[:final]≈N atol=testTol
    end
end

println("DONE with ODE tests")
include(joinpath(@__DIR__, "../src/DAETS_utils.jl"))
# include(joinpath(@__DIR__, "../src/TaylorSeries_caches.jl"))
println("Starting tests on DAETS")

# --- Test Cases --- #
@testset "Signature Matrix & Derivative Order Tests" begin
    @syms t x(t) y(t)
    # # Test max_derivative_order
    # The variable itself.
    order1 = max_derivative_order(x(t), x, t)
    @test order1 == 0

    # A multiplication expression: 7*x(t)
    order2 = max_derivative_order(7*x(t), x, t)
    @test order2 == 0

    # First derivative: x'(t)
    dx = Differential(t)(x(t))
    # println(typeof(dx))
    order3 = max_derivative_order(dx, x, t)
    @test order3 == 1

    # Second derivative: x''(t)
    d2x = Differential(t)(dx)
    order4 = max_derivative_order(d2x, x, t)
    @test order4 == 2

    # Expression that does not contain x(t): y(t)
    order5 = max_derivative_order(y(t), x, t)
    @test order5 == -Inf

    # Test signature_matrix:
    # Equation 1: f₁ = x'(t) + sin(x(t)) = 0
    # Equation 2: f₂ = y''(t) - x(t) = 0
    eq1 = Differential(t)(x(t)) + sin(x(t))
    eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
    eqs = [eq1, eq2]
    # pass functions x and y (which take input t)
    vars = [x, y]

    Σ = signature_matrix(eqs, vars, t)
    # First equation:
    # Column 1: x(t) appears as x'(t) ⟹ order 1.
    # Column 2: no y(t) ⟹ -Inf.
    @test Σ[1, 1] == 1
    @test Σ[1, 2] == -Inf

    # Second equation:
    # Column 1: x(t) appears as x(t) ⟹ order 0.
    # Column 2: y(t) appears as y''(t) ⟹ order 2.
    @test Σ[2, 1] == 0
    @test Σ[2, 2] == 2
end

println("DONE Signature Matrix Tests")

@testset "Highest Value Transversal Tests" begin
    @syms t x(t) y(t)

    # Same Equation
    eq1 = Differential(t)(x(t)) + sin(x(t))
    eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
    eqs = [eq1, eq2]
    vars = [x, y]
    Σ = signature_matrix(eqs, vars, t)

    # Expected signature matrix:
    # [ 1  -Inf
    #   0     2 ]
    @test Σ == [1 -Inf; 0 2]
    # Test HVT
    transversal, value = highest_value_transversal(Σ)

    # Expected transversal: [(1, 1), (2, 2)] (x'(t) in eq1 and y''(t) in eq2)
    # Expected value: 1 + 2 = 3 TODO: Should probably add a check in case the value is infinite. In this case we say that the system is "ill-posed".
    @test transversal == [(1, 1), (2, 2)]
    @test value == 3
end
@testset "Highest Value Transversal Tests for 3x3 System" begin
    @syms t x(t) y(t) z(t)

    # Equation 1: f₁ = x''(t) + y(t) = 0
    # Equation 2: f₂ = y'(t) + z(t) = 0
    # Equation 3: f₃ = z''(t) + x(t) = 0
    eq1 = Differential(t)(Differential(t)(x(t))) + y(t)
    eq2 = Differential(t)(y(t)) + z(t)
    eq3 = Differential(t)(Differential(t)(z(t))) + x(t)
    eqs = [eq1, eq2, eq3]
    vars = [x, y, z]
    Σ = signature_matrix(eqs, vars, t)

    # Expected signature matrix:
    # [ 2   0  -Inf
    #  -Inf  1     0
    #   0  -Inf    2 ]
    @test Σ == [2 0 -Inf; -Inf 1 0; 0 -Inf 2]
    # Test HVT
    transversal, value = highest_value_transversal(Σ)

    # Expected transversal: [(1, 1), (2, 2), (3, 3)] (x''(t) in eq1, y'(t) in eq2, z''(t) in eq3)
    # Expected value: 2 + 1 + 2 = 5
    @test transversal == [(1, 1), (2, 2), (3, 3)]
    @test value == 5
end
println("DONE Highest Value Transversal Tests")

@testset "Find Offsets Tests for 2x2 System" begin
    @syms t x(t) y(t)

    # Same Equation as 2x2 System
    eq1 = Differential(t)(x(t)) + sin(x(t))
    eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
    eqs = [eq1, eq2]
    vars = [x, y]
    Σ = signature_matrix(eqs, vars, t)

    # Find the highest value transversal
    transversal, value = highest_value_transversal(Σ)

    # Find the offsets
    c, d = find_offsets(Σ, transversal)

    # Expected offsets (canonical):
    # c = [0, 0]
    # d = [1, 2]
    @test c == [0, 0]
    @test d == [1, 2]
end
@testset "Find Offsets Tests for 3x3 System" begin
    @syms t x(t) y(t) z(t)

    # Same Equation as 3x3 System
    eq1 = Differential(t)(Differential(t)(x(t))) + y(t)
    eq2 = Differential(t)(y(t)) + z(t)
    eq3 = Differential(t)(Differential(t)(z(t))) + x(t)
    eqs = [eq1, eq2, eq3]
    vars = [x, y, z]
    Σ = signature_matrix(eqs, vars, t)
    transversal, value = highest_value_transversal(Σ)

    # Test Offsets
    c, d = find_offsets(Σ, transversal)

    # Expected offsets (canonical):
    # c = [0, 0, 0]
    # d = [2, 1, 2]
    @test c == [0, 0, 0]
    @test d == [2, 1, 2]
end

println("DONE Find Offsets Tests")

@testset "System Jacobian Tests for 2x2 System" begin
    @syms t x(t) y(t)

    # Same 2x2
    eq1 = Differential(t)(x(t)) + sin(x(t))
    eq2 = Differential(t)(Differential(t)(y(t))) - x(t)
    eqs = [eq1, eq2]
    vars = [x, y]
    Σ = signature_matrix(eqs, vars, t)
    transversal, value = highest_value_transversal(Σ)
    c, d = find_offsets(Σ, transversal)

    # Convert c and d to Vector{Int}
    c = Int.(c)
    d = Int.(d)

    # Test Jacobian
    J = system_jacobian(eqs, vars, t, c, d, Σ)

    # Expected Jacobian:
    # [ 1   0
    #   0   1 ]
    @test isequal(J[1, 1], 1)
    @test isequal(J[1, 2], 0)
    @test isequal(J[2, 1], 0)
    @test isequal(J[2, 2], 1)
end
@testset "System Jacobian Tests for 3x3 System" begin
    @syms t x(t) y(t) z(t)

    # Same 3x3
    eq1 = Differential(t)(Differential(t)(x(t))) + y(t)
    eq2 = Differential(t)(y(t)) + z(t)
    eq3 = Differential(t)(Differential(t)(z(t))) + x(t)
    eqs = [eq1, eq2, eq3]
    vars = [x, y, z]
    Σ = signature_matrix(eqs, vars, t)
    transversal, value = highest_value_transversal(Σ)
    c, d = find_offsets(Σ, transversal)

    # Convert c and d to Vector{Int}
    c = Int.(c)
    d = Int.(d)

    # Test Jacobian
    J = system_jacobian(eqs, vars, t, c, d, Σ)

    # Expected Jacobian:
    # [1   0   0
    #  0   1   1
    #  0   0   1]
    @test isequal(J[1, 1], 1)
    @test isequal(J[1, 2], 0)
    @test isequal(J[1, 3], 0)
    @test isequal(J[2, 1], 0)
    @test isequal(J[2, 2], 1)
    @test isequal(J[2, 3], 0)
    @test isequal(J[3, 1], 0)
    @test isequal(J[3, 2], 0)
    @test isequal(J[3, 3], 1)
end
println("DONE System Jacobian Tests")

@testset "System Jacobian Tests for Simple Pendulum" begin
    @syms t x(t) y(t) λ(t) G L

    # Equations
    f = Differential(t)(Differential(t)(x(t))) + x(t) * λ(t)
    g = Differential(t)(Differential(t)(y(t))) + y(t) * λ(t) - G
    h = x(t)^2 + y(t)^2 - L^2
    eqs = [f, g, h]
    vars = [x, y, λ]

    # Construct the signature matrix
    Σ = signature_matrix(eqs, vars, t)

    # Find the highest value transversal
    transversal, value = highest_value_transversal(Σ)

    # Find the offsets
    c, d = find_offsets(Σ, transversal)

    # Convert c and d to Vector{Int}
    c = Int.(c)
    d = Int.(d)

    # Construct the system Jacobian
    J = system_jacobian(eqs, vars, t, c, d, Σ)

    # Expected Jacobian:
    # [1   0   x(t)
    #  0   1   y(t)
    #  2   2   0]
    @test isequal(J[1, 1], 1)
    @test isequal(J[1, 2], 0)
    @test isequal(J[1, 3], x(t))
    @test isequal(J[2, 1], 0)
    @test isequal(J[2, 2], 1)
    @test isequal(J[2, 3], y(t))
    @test isequal(J[3, 1], 2x(t))
    @test isequal(J[3, 2], 2y(t))
    @test isequal(J[3, 3], 0)
end
# @testset "compute_taylor_coefficients! Tests" begin
#     @syms t x(t) y(t)
#     @syms G L

#     # Define a simple ODE system: x'(t) = y(t), y'(t) = -x(t)
#     f = [Differential(t)(x(t)) - y(t), Differential(t)(y(t)) + x(t)]
#     vars = [x, y]

#     # Create a mock integrator and cache
#     integrator = (u = [1.0, 0.0], t = 0.0, dt = 0.1, f = f, p = nothing)
#     cache = (Σ = nothing, c = nothing, d = nothing, J = nothing, xTS = nothing, xtrial = nothing, htrial = 0.1, e = 0.0, tmp = nothing)

#     # Compute Taylor coefficients
#     xcur = compute_taylor_coefficients!(integrator, cache)

#     # Expected Taylor coefficients for x(t) and y(t) at t = 0:
#     # x(t) = 1 - t^2/2 + t^4/24 - ...
#     # y(t) = t - t^3/6 + t^5/120 - ...
#     @test xcur[1] ≈ 1.0  # x(0)
#     @test xcur[2] ≈ 0.0  # y(0)
# end
