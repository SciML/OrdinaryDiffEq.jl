using OrdinaryDiffEqVerner, OrdinaryDiffEqCore, DiffEqBase, Test
using LinearAlgebra
using DiffEqDevTools, Random, Plots
import OrdinaryDiffEqLowStorageRK
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear

# Define the ODE function (in-place and out-of-place forms)
function f!(du, u, p, t)
    du[1] = -u[1]
end

function f(u, p, t)
    -u
end

# Problem setup
t_end = 64.0
testTol = 1.0
setprecision(256)

# Use analytic solution for error calculation
prob_oop = ODEProblem(
    ODEFunction(f; analytic = (u0, p, t) -> exp(-t)),
    1.0,
    (0.0, t_end)
)

function test_algorithm_convergence(prob, alg, dts, t_end; testTol=1.0)
    alg_name = string(typeof(alg).name.name)
    expected_order = OrdinaryDiffEqVerner.alg_order(alg)

    # Use DiffEqDevTools to estimate overall convergence order
    sim = test_convergence(dts, prob, alg)

    # Test that the convergence order is close to expected
    @test sim.𝒪est[:final] ≈ expected_order atol=testTol
end

# Define algorithms to test
algorithms = [
    # Vern6(),
    Vern7(),
    Vern8(),
    Vern9(),
    RKV76IIa()
]

# Define timesteps for testing
dts = [1, 0.5, 0.25, 0.125]

# Test each algorithm

for alg in algorithms
    test_algorithm_convergence(prob_oop, alg, dts, t_end; testTol=testTol)
end