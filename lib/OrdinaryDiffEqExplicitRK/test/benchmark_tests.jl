# Large System Benchmark: ExplicitRK vs Tsit5
# Comparing performance on systems with 10,000+ ODEs

using OrdinaryDiffEq
using BenchmarkTools
using Printf
using OrdinaryDiffEqExplicitRK
using DiffEqBase
using LinearAlgebra
using Plots


include("../src/algorithms.jl")

# Import tableau construction functions
using OrdinaryDiffEqExplicitRK: constructDormandPrince
using OrdinaryDiffEqExplicitRK: constructTsit5ExplicitRKSimple

using OrdinaryDiffEqExplicitRK: constructTsit5ExplicitRK
# ============================================================================
# Problem 1: Lorenz-96 System (Atmospheric Model)
# ============================================================================
"""
    lorenz96!(du, u, p, t)

The Lorenz-96 model - a chaotic dynamical system used in atmospheric science.
Scalable to any dimension N.

du[i]/dt = (u[i+1] - u[i-2]) * u[i-1] - u[i] + F

where F is a forcing parameter (typically F = 8)
"""
function lorenz96!(du, u, p, t)
    F = p
    N = length(u)

    # Periodic boundary conditions
    @inbounds for i in 1:N
        i_m2 = mod1(i - 2, N)
        i_m1 = mod1(i - 1, N)
        i_p1 = mod1(i + 1, N)

        du[i] = (u[i_p1] - u[i_m2]) * u[i_m1] - u[i] + F
    end
    return nothing
end

function create_lorenz96_problem(N=10000; F=8.0, tspan=(0.0, 10.0))
    # Initial conditions: small random perturbations around F
    u0 = F .+ 0.01 .* randn(N)
    return ODEProblem(lorenz96!, u0, tspan, F)
end

# ============================================================================
# Problem 2: 1D Reaction-Diffusion System (PDE Discretization)
# ============================================================================
"""
    reaction_diffusion!(du, u, p, t)

1D reaction-diffusion PDE discretized with finite differences:
∂u/∂t = D * ∂²u/∂x² + R(u)

where R(u) = α*u*(1-u) is a reaction term (Fisher-KPP equation)
"""
function reaction_diffusion!(du, u, p, t)
    D, α, dx = p
    N = length(u)
    inv_dx2 = 1.0 / (dx * dx)

    # Neumann boundary conditions (zero flux)
    @inbounds du[1] = D * (u[2] - 2*u[1] + u[1]) * inv_dx2 + α * u[1] * (1 - u[1])

    @inbounds for i in 2:N-1
        du[i] = D * (u[i+1] - 2*u[i] + u[i-1]) * inv_dx2 + α * u[i] * (1 - u[i])
    end

    @inbounds du[N] = D * (u[N] - 2*u[N] + u[N-1]) * inv_dx2 + α * u[N] * (1 - u[N])

    return nothing
end

function create_reaction_diffusion_problem(N=10000; D=0.01, α=1.0, tspan=(0.0, 5.0))
    # Spatial discretization
    L = 100.0  # Domain length
    dx = L / (N - 1)

    # Initial condition: step function
    u0 = zeros(N)
    u0[N÷4:3*N÷4] .= 1.0

    return ODEProblem(reaction_diffusion!, u0, tspan, (D, α, dx))
end

# ============================================================================
# Problem 3: Coupled Oscillators (Network Dynamics)
# ============================================================================
"""
    coupled_oscillators!(du, u, p, t)

Network of coupled harmonic oscillators:
d²x[i]/dt² = -ω²*x[i] + K*Σ(x[j] - x[i]) for connected j

Converted to first-order system: [x₁, v₁, x₂, v₂, ..., xₙ, vₙ]
"""
function coupled_oscillators!(du, u, p, t)
    ω, K = p
    N = length(u) ÷ 2

    @inbounds for i in 1:N
        x_idx = 2*i - 1
        v_idx = 2*i

        # Position derivative = velocity
        du[x_idx] = u[v_idx]

        # Velocity derivative = acceleration
        # Couple to neighbors (ring topology)
        i_prev = mod1(i - 1, N)
        i_next = mod1(i + 1, N)

        x_prev_idx = 2*i_prev - 1
        x_next_idx = 2*i_next - 1

        coupling = K * ((u[x_prev_idx] - u[x_idx]) + (u[x_next_idx] - u[x_idx]))
        du[v_idx] = -ω*ω*u[x_idx] + coupling
    end

    return nothing
end

function create_coupled_oscillators_problem(N=5000; ω=1.0, K=0.1, tspan=(0.0, 50.0))
    # N oscillators → 2N equations (position + velocity)
    u0 = zeros(2*N)
    # Random initial positions and velocities
    u0[1:2:end] = randn(N)  # positions
    u0[2:2:end] = randn(N)  # velocities

    return ODEProblem(coupled_oscillators!, u0, tspan, (ω, K))
end

# ============================================================================
# Problem 4: Lotka Volterra System (prey-predator Model)
# ============================================================================
"""
    lotka_volterra_blockdiag!(du, u, p, t)

Block-diagonal system of independent Lotka-Volterra pairs.
Each pair: du[2i-1] = α*u[2i-1] - γ*u[2i-1]*u[2i]
           du[2i]   = -β*u[2i] + δ*u[2i-1]*u[2i]
"""
function lotka_volterra_blockdiag!(du, u, p, t)
    α, β, γ, δ, npairs = p
    @inbounds for i in 1:npairs
        xidx = 2i - 1
        yidx = 2i
        x, y = u[xidx], u[yidx]
        du[xidx] = α * x - γ * x * y
        du[yidx] = -β * y + δ * x * y
    end
    return nothing
end

function create_lotka_volterra_problem(N=10000; α=1.5, β=3.0, γ=1.0, δ=1.0, tspan=(0.0, 10.0))
    @assert N % 2 == 0 "N must be even (pairs of equations)"
    npairs = N ÷ 2
    u0 = repeat([1.0, 1.0], npairs)
    p = (α, β, γ, δ, npairs)
    return ODEProblem(lotka_volterra_blockdiag!, u0, tspan, p)
end

# ============================================================================
# Benchmark Utilities
# ============================================================================

function benchmark_solver(prob, alg; name="Solver")
    println("\n" * "="^70)
    println("Benchmarking: $name")
    println("="^70)

    # Warmup
    sol = solve(prob, alg, saveat=1.0, abstol=1e-6, reltol=1e-3)

    # Benchmark
    println("\nRunning benchmark (5 samples, 1 evaluation each)...")
    b = @benchmark solve($prob, $alg, saveat=1.0, abstol=1e-6, reltol=1e-3) samples=5 evals=1

    println("\nResults:")
    println("  Minimum time:  $(BenchmarkTools.prettytime(minimum(b.times)))")
    println("  Median time:   $(BenchmarkTools.prettytime(median(b.times)))")
    println("  Mean time:     $(BenchmarkTools.prettytime(mean(b.times)))")
    println("  Memory:        $(BenchmarkTools.prettymemory(b.memory))")
    println("  Allocations:   $(b.allocs)")
    println("\nSolution stats:")
    println("  Steps:         $(sol.stats.naccept)")
    println("  Function evals: $(sol.stats.nf)")
    println("  Final time:    $(sol.t[end])")

    return (benchmark=b, solution=sol)
end

function compare_solvers(prob_name, prob; use_tsit5_tableau=false)
    println("\n" * "#"^70)
    println("# Problem: $prob_name")
    println("# System size: $(length(prob.u0)) equations")
    println("#"^70)

    # Benchmark ExplicitRK with chosen tableau
    if use_tsit5_tableau
        tableau = constructTsit5ExplicitRKSimple()
        tableau_name = "Tsit5 Tableau"
    else
        tableau = constructDormandPrince()
        tableau_name = "Dormand-Prince"
    end

    result_explicit = benchmark_solver(prob, ExplicitRK(tableau=tableau), name="ExplicitRK ($tableau_name)")
    result_tsit5 = benchmark_solver(prob, Tsit5(), name="Tsit5 (Specialized)")
    b_explicit, sol_explicit = result_explicit.benchmark, result_explicit.solution
    b_tsit5, sol_tsit5 = result_tsit5.benchmark, result_tsit5.solution

    # Compare
    println("\n" * "="^70)
    println("COMPARISON")
    println("="^70)

    time_ratio = median(b_explicit.times) / median(b_tsit5.times)
    mem_ratio = b_explicit.memory / max(b_tsit5.memory, 1)  # Avoid div by 0
    alloc_ratio = b_explicit.allocs / max(b_tsit5.allocs, 1)

    @printf("  Speedup (Tsit5 vs ExplicitRK): %.2fx\n", time_ratio)
    @printf("  Memory reduction:              %.2fx\n", mem_ratio)
    @printf("  Allocation reduction:          %.2fx\n", alloc_ratio)

    return (explicit=(benchmark=b_explicit, solution=sol_explicit), tsit5=(benchmark=b_tsit5, solution=sol_tsit5))
end

# ============================================================================
# Main Benchmark Suite
# ============================================================================

function run_large_system_benchmarks(; use_tsit5_tableau=false)
    println("="^70)
    println("LARGE ODE SYSTEM BENCHMARKS")
    println("Comparing ExplicitRK vs Tsit5 on 10,000+ equation systems")
    if use_tsit5_tableau
        println("Using Tsit5 tableau for ExplicitRK")
    else
        println("Using Dormand-Prince tableau for ExplicitRK")
    end
    println("="^70)


    # Test 1: Lorenz-96 (10,000 equations)
    println("\n📊 Test 1: Lorenz-96 Atmospheric Model")
    prob1 = create_lorenz96_problem(10000, tspan=(0.0, 5.0))
    results1 = compare_solvers("Lorenz-96 (N=10,000)", prob1, use_tsit5_tableau=use_tsit5_tableau)

    # Commented out these tests based on feedback

    # Test 2: Reaction-Diffusion (10,000 equations)
    # println("\n📊 Test 2: Reaction-Diffusion PDE")
    # prob2 = create_reaction_diffusion_problem(10000, tspan=(0.0, 2.0))
    # results2 = compare_solvers("Reaction-Diffusion 1D (N=10,000)", prob2, use_tsit5_tableau=use_tsit5_tableau)

    # # Test 3: Coupled Oscillators (10,000 equations = 5,000 oscillators)
    # println("\n📊 Test 3: Coupled Oscillators Network")
    # prob3 = create_coupled_oscillators_problem(5000, tspan=(0.0, 10.0))
    # results3 = compare_solvers("Coupled Oscillators (5,000 × 2 = 10,000 eqs)", prob3, use_tsit5_tableau=use_tsit5_tableau)

    # Test 4: Lotka-Volterra (10,000 equations = 5,000 pairs)
    println("\n📊 Test 4: Lotka-Volterra Block-Diagonal System")
    prob4 = create_lotka_volterra_problem(10000, tspan=(0.0, 10.0))
    results4 = compare_solvers("Lotka-Volterra (5,000 × 2 = 10,000 eqs)", prob4, use_tsit5_tableau=use_tsit5_tableau)

    println("\n" * "="^70)
    println("BENCHMARK SUITE COMPLETE")
    println("="^70)

    # return (lorenz96=results1, reaction_diffusion=results2, oscillators=results3, lotka_volterra=results4)
    return (lorenz96=results1, lotka_volterra=results4)
end

# ============================================================================
# Quick Test (Smaller Systems for Development)
# ============================================================================

function quick_test()
    println("\n🔬 Quick Test (smaller systems for development)")
    println("="^70)

    # Small Lorenz-96
    prob = create_lorenz96_problem(100, tspan=(0.0, 1.0))

    println("\nTesting ExplicitRKTsit5Tableau...")
    # dp_tableau = constructDormandPrince()
    dp_tableau = constructTsit5ExplicitRK()
    @time sol1 = solve(prob, ExplicitRK(tableau=dp_tableau), abstol=1e-6, reltol=1e-5)
    println("  Steps: $(sol1.stats.naccept), Function evals: $(sol1.stats.nf)")

    println("\nTesting Tsit5...")
    @time sol2 = solve(prob, Tsit5(), abstol=1e-6, reltol=1e-5)
    println("  Steps: $(sol2.stats.naccept), Function evals: $(sol2.stats.nf)")

    println("\n✅ Quick test passed!")
end

function run_small_system_benchmarks(; use_tsit5_tableau=true)
    println("\n" * "="^70)
    println("SMALL ODE SYSTEM BENCHMARKS")
    println("Comparing ExplicitRK vs Tsit5 on small systems")
    println("="^70)

    # Test 1: Small Lorenz-96 (100 equations)
    println("\n📊 Test 1: Lorenz-96 Atmospheric Model (N=100)")
    prob1 = create_lorenz96_problem(100, tspan=(0.0, 5.0))
    results1 = compare_solvers("Lorenz-96 (N=100)", prob1, use_tsit5_tableau=use_tsit5_tableau)

    # Test 2: Small Reaction-Diffusion (100 equations)
    println("\n📊 Test 2: Reaction-Diffusion PDE (N=100)")
    prob2 = create_reaction_diffusion_problem(100, tspan=(0.0, 2.0))
    results2 = compare_solvers("Reaction-Diffusion 1D (N=100)", prob2, use_tsit5_tableau=use_tsit5_tableau)

    # Test 3: Small Coupled Oscillators (100 oscillators × 2 = 200 equations)
    println("\n📊 Test 3: Coupled Oscillators Network (N=100 oscillators)")
    prob3 = create_coupled_oscillators_problem(100, tspan=(0.0, 10.0))
    results3 = compare_solvers("Coupled Oscillators (100 × 2 = 200 eqs)", prob3, use_tsit5_tableau=use_tsit5_tableau)

    # Test 4: Small Lotka-Volterra (100 pairs × 2 = 200 equations)
    println("\n📊 Test 4: Lotka-Volterra Block-Diagonal System (N=200)")
    prob4 = create_lotka_volterra_problem(200, tspan=(0.0, 10.0))
    results4 = compare_solvers("Lotka-Volterra (100 × 2 = 200 eqs)", prob4, use_tsit5_tableau=use_tsit5_tableau)

    println("\n" * "="^70)
    println("BENCHMARK SUITE COMPLETE")
    println("="^70)

    return (lorenz96=results1, reaction_diffusion=results2, oscillators=results3, lotka_volterra=results4)
end

# After running benchmarks
results = run_large_system_benchmarks(use_tsit5_tableau=true)
# results = run_small_system_benchmarks()
sol_explicit = results.lotka_volterra.explicit.solution
sol_tsit5 = results.lotka_volterra.tsit5.solution
# Compare final state
final_explicit = sol_explicit.u[end]
final_tsit5 = sol_tsit5.u[end]
diff_norm = norm(final_explicit .- final_tsit5)
println("Norm of difference between ExplicitRK and Tsit5 final values: $diff_norm")
# Compare all time points
all_diffs = [norm(sol_explicit.u[i] .- sol_tsit5.u[i]) for i in 1:length(sol_explicit.u)]
println("Max difference across all time points: $(maximum(all_diffs))")

# Quick verification
plot(sol_explicit.t, sol_explicit.u, label="ExplicitRK")
plot!(sol_tsit5.t, sol_tsit5.u, label="Tsit5")
