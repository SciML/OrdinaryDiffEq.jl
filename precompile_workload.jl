using OrdinaryDiffEq
println("Precompiling OrdinaryDiffEq.jl with key solver workloads...")

# Simple scalar ODE for basic functionality
f_scalar(u, p, t) = -u
prob_scalar = ODEProblem(f_scalar, 1.0, (0.0, 1.0))

# Lorenz system - classic benchmark problem
function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8/3) * u[3]
end
u0_lorenz = [1.0, 0.0, 0.0]
tspan_lorenz = (0.0, 10.0)
prob_lorenz = ODEProblem(lorenz!, u0_lorenz, tspan_lorenz)

# Van der Pol oscillator - stiff problem
function vanderpol!(du, u, p, t)
    μ = p[1]
    du[1] = u[2]
    du[2] = μ * (1 - u[1]^2) * u[2] - u[1]
end
u0_vdp = [2.0, 0.0]
tspan_vdp = (0.0, 6.3)
p_vdp = [1000.0]  # stiff parameter
prob_vdp = ODEProblem(vanderpol!, u0_vdp, tspan_vdp, p_vdp)

# Large system for testing high-dimensional problems
function large_system!(du, u, p, t)
    n = length(u)
    for i in 1:n
        du[i] = -u[i] + 0.1 * sin(t * i)
    end
end
u0_large = rand(100)
prob_large = ODEProblem(large_system!, u0_large, (0.0, 1.0))

println("  - Precompiling Tsit5 (explicit Runge-Kutta)...")
# Tsit5 - Most popular explicit solver
solve(prob_scalar, Tsit5(), reltol=1e-8, abstol=1e-10)
solve(prob_lorenz, Tsit5(), reltol=1e-8, abstol=1e-10)
solve(prob_large, Tsit5(), reltol=1e-6, abstol=1e-8)

println("  - Precompiling FBDF (multistep method)...")
# FBDF - Backward differentiation formula for stiff problems
solve(prob_scalar, FBDF(), reltol=1e-8, abstol=1e-10)
solve(prob_vdp, FBDF(), reltol=1e-6, abstol=1e-8)
solve(prob_large, FBDF(), reltol=1e-6, abstol=1e-8)

println("  - Precompiling Rodas5P (Rosenbrock method)...")
# Rodas5P - Rosenbrock method for stiff problems
solve(prob_scalar, Rodas5P(), reltol=1e-8, abstol=1e-10)
solve(prob_vdp, Rodas5P(), reltol=1e-6, abstol=1e-8)
solve(prob_lorenz, Rodas5P(), reltol=1e-8, abstol=1e-10)

println("  - Precompiling Vern9 (high-order explicit)...")
# Vern9 - High-order explicit method for non-stiff problems
solve(prob_scalar, Vern9(), reltol=1e-12, abstol=1e-14)
solve(prob_lorenz, Vern9(), reltol=1e-10, abstol=1e-12)

println("  - Precompiling common callback patterns...")
# Common callback patterns
condition(u, t, integrator) = u[1] - 1.0
affect!(integrator) = nothing
cb = ContinuousCallback(condition, affect!)
solve(prob_lorenz, Tsit5(), callback=cb, reltol=1e-8)

# Discrete callback
discrete_condition(u, t, integrator) = t == 5.0
discrete_cb = DiscreteCallback(discrete_condition, affect!)
solve(prob_lorenz, Tsit5(), callback=discrete_cb, reltol=1e-8)

println("  - Precompiling dense output and interpolation...")
# Dense output - commonly used for plotting/analysis
sol = solve(prob_lorenz, Tsit5(), dense=true, reltol=1e-8)
sol(5.0)  # Interpolation

println("  - Precompiling saveat and save patterns...")
# Saveat - specified save times
solve(prob_lorenz, Tsit5(), saveat=0.1, reltol=1e-8)

# Save only at end
solve(prob_lorenz, Tsit5(), save_everystep=false, save_start=false, reltol=1e-8)

println("  - Precompiling in-place vs out-of-place...")
# Out-of-place version
f_oop(u, p, t) = [-u[1] + u[2], -u[2]]
prob_oop = ODEProblem(f_oop, [1.0, 0.5], (0.0, 1.0))
solve(prob_oop, Tsit5(), reltol=1e-8)

println("  - Precompiling mass matrix problems...")
# Mass matrix problems
function mm_f!(du, u, p, t)
    du[1] = u[2]
    du[2] = -u[1]
end
M = [1.0 0.0; 0.0 0.0]  # Singular mass matrix
prob_mm = ODEProblem(ODEFunction(mm_f!, mass_matrix=M), [1.0, 0.0], (0.0, 1.0))
solve(prob_mm, Rodas5P(), reltol=1e-8)

println("  - Precompiling parameter sensitivity...")
# Parameter-dependent problem (common in sensitivity analysis)
function param_f!(du, u, p, t)
    du[1] = p[1] * u[1] + p[2] * u[2]
    du[2] = p[3] * u[1] - p[4] * u[2]
end
prob_param = ODEProblem(param_f!, [1.0, 0.0], (0.0, 1.0), [0.1, 0.2, -0.3, 0.4])
solve(prob_param, Tsit5(), reltol=1e-8)
solve(prob_param, Rodas5P(), reltol=1e-8)

println("Precompilation workload completed successfully!")
