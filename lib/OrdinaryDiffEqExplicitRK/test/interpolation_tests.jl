"""
    ExplicitRK Tests -

"""

using OrdinaryDiffEqExplicitRK
using OrdinaryDiffEqExplicitRK: constructTsit5ExplicitRK, constructTsit5ExplicitRKSimple
using OrdinaryDiffEqCore
using DiffEqBase
using Test
using Random
import SciMLBase

Random.seed!(100)

# ============================================================================
# Test Problems
# ============================================================================

# Scalar linear ODE: y' = 1.01y
f_linear(u, p, t) = 1.01 * u
prob_ode_linear = ODEProblem(
    ODEFunction(f_linear; analytic=(u0, p, t) -> u0 * exp(1.01t)),
    1/2, (0.0, 1.0)
)

# 2D linear ODE (in-place)
function f_2Dlinear!(du, u, p, t)
    du[1] = 1.01 * u[1]
    du[2] = 1.01 * u[2]
end
prob_ode_2Dlinear = ODEProblem(
    ODEFunction(f_2Dlinear!; analytic=(u0, p, t) -> u0 .* exp(1.01t)),
    [1/2, 1/2], (0.0, 1.0)
)

# ============================================================================
# Basic Solve Tests
# ============================================================================

prob = prob_ode_linear
sol = solve(prob, ExplicitRK(tableau=constructTsit5ExplicitRK()))
@test length(sol) < 20
@test SciMLBase.successful_retcode(sol)

sol = solve(prob, ExplicitRK(tableau=constructTsit5ExplicitRKSimple()))
@test length(sol) < 20
@test SciMLBase.successful_retcode(sol)

prob = prob_ode_2Dlinear
sol = solve(prob, ExplicitRK(tableau=constructTsit5ExplicitRK()))
@test length(sol) < 20
@test SciMLBase.successful_retcode(sol)

sol = solve(prob, ExplicitRK(tableau=constructTsit5ExplicitRKSimple()))
@test length(sol) < 20
@test SciMLBase.successful_retcode(sol)

# ============================================================================
# Tableau Comparison Tests
# ============================================================================

# Compare Full vs Simple tableau - should give same results
prob = prob_ode_linear
tableau_full = constructTsit5ExplicitRK()
tableau_simple = constructTsit5ExplicitRKSimple()

sol1 = solve(prob, ExplicitRK(tableau=tableau_full); dt=1/2^6, adaptive=false, save_everystep=false)
sol2 = solve(prob, ExplicitRK(tableau=tableau_simple); dt=1/2^6, adaptive=false, save_everystep=false)
@test abs(sol1.u[end] - sol2.u[end]) < 1e-10

prob = prob_ode_2Dlinear
sol1 = solve(prob, ExplicitRK(tableau=tableau_full); dt=1/2^3, adaptive=false, save_everystep=false)
sol2 = solve(prob, ExplicitRK(tableau=tableau_simple); dt=1/2^3, adaptive=false, save_everystep=false)
@test minimum(abs.(sol1.u[end] - sol2.u[end]) .< 1e-10)

# Adaptive stepping should give same length
sol1 = solve(prob, ExplicitRK(tableau=tableau_full); dt=1/2^6)
sol2 = solve(prob, ExplicitRK(tableau=tableau_simple); dt=1/2^6)
@test length(sol1) == length(sol2)
@test SciMLBase.successful_retcode(sol1)
@test SciMLBase.successful_retcode(sol2)

# ============================================================================
# Interpolation Tests
# ============================================================================

# Test problem with known analytical solution
f(u, p, t) = -u
prob_interp = ODEProblem(ODEFunction(f; analytic=(u0, p, t) -> exp(-t)), 1.0, (0.0, 1.0))


sol = solve(prob_interp, ExplicitRK(tableau=constructTsit5ExplicitRK()), dense=true)
@test SciMLBase.successful_retcode(sol)
@test sol(0.0) ≈ 1.0 atol=1e-8
@test sol(0.5) ≈ exp(-0.5) atol=0.1
@test sol(1.0) ≈ exp(-1.0) atol=0.1


function f_vec!(du, u, p, t)
    du[1] = -u[1]
    du[2] = -2*u[2]
end
exact_vec(t) = [exp(-t), exp(-2t)]
prob_interp_vec = ODEProblem(ODEFunction(f_vec!; analytic=(u0, p, t) -> exact_vec(t)), [1.0, 1.0], (0.0, 1.0))

sol = solve(prob_interp_vec, ExplicitRK(tableau=constructTsit5ExplicitRK()), dense=true)
@test SciMLBase.successful_retcode(sol)
@test sol(0.5) ≈ exact_vec(0.5) atol=0.1

# Interpolation with idxs
@test sol(0.5, idxs=1) ≈ exp(-0.5) atol=0.1
@test sol(0.5, idxs=2) ≈ exp(-1.0) atol=0.1
@test sol(0.5, idxs=[1,2]) ≈ exact_vec(0.5) atol=0.1

# Interpolation at step boundaries
@test sol(0.0) ≈ [1.0, 1.0] atol=1e-6
@test sol(1.0) ≈ sol.u[end] atol=1e-5
