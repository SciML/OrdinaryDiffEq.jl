# Usage Guide

This page provides guidance on using StochasticDiffEq.jl effectively.

## Basic Usage

### Problem Definition

StochasticDiffEq.jl uses the standard DifferentialEquations.jl problem interface:

```julia
using StochasticDiffEq

# For scalar problems
function f(u, p, t)  # drift
    return μ * u
end

function g(u, p, t)  # diffusion  
    return σ * u
end

# For in-place systems
function f!(du, u, p, t)
    du[1] = μ * u[1]
    du[2] = -ν * u[2]
end

function g!(du, u, p, t)
    du[1] = σ₁ * u[1]
    du[2] = σ₂ * u[2]
end

# Create problem
prob = SDEProblem(f, g, u0, tspan)
```

### Solver Selection

Choose solvers based on your problem characteristics:

```julia
# Default - good for most problems
sol = solve(prob)

# Specify solver explicitly
sol = solve(prob, SOSRI())          # Recommended for diagonal noise
sol = solve(prob, SOSRA())          # Optimal for additive noise  
sol = solve(prob, SKenCarp())       # For stiff problems
sol = solve(prob, EM())             # For maximum efficiency
```

### Algorithm Parameters

Most solvers accept parameters for customization:

```julia
# Euler-Maruyama with step splitting
sol = solve(prob, EM(split = true))

# RKMilCommute with Stratonovich interpretation
sol = solve(prob, RKMilCommute(interpretation = :Stratonovich))

# Implicit methods with solver options
sol = solve(prob, SKenCarp(linsolve = KrylovJL_GMRES()))
```

## Tolerances and Adaptive Stepping

Set absolute and relative tolerances:

```julia
sol = solve(prob, SOSRI(), abstol = 1e-6, reltol = 1e-3)
```

For fixed time stepping:

```julia
sol = solve(prob, EM(), dt = 0.01, adaptive = false)
```

## Noise Types

### Diagonal Noise

Most common case - each component has independent noise:

```julia
function g!(du, u, p, t)
    du[1] = σ₁ * u[1]
    du[2] = σ₂ * u[2]
end
```

### Scalar Noise

Single noise source affects all components:

```julia
function g!(du, u, p, t)
    du[1] = σ * u[1]
    du[2] = σ * u[2]
end
```

### Non-diagonal Noise

Multiple noise sources with cross-terms:

```julia
function g!(du, u, p, t)
    du[1] = σ₁₁ * u[1] + σ₁₂ * u[2]
    du[2] = σ₂₁ * u[1] + σ₂₂ * u[2]
end
```

### Additive Noise

Noise independent of solution:

```julia
function g!(du, u, p, t)
    du[1] = σ₁
    du[2] = σ₂
end
```

## Itô vs Stratonovich

Specify interpretation when creating problems or choosing solvers:

```julia
# Itô interpretation (default)
prob = SDEProblem(f!, g!, u0, tspan, interpretation = :Ito)

# Stratonovich interpretation  
prob = SDEProblem(f!, g!, u0, tspan, interpretation = :Stratonovich)

# Or at solver level
sol = solve(prob, RKMil(interpretation = :Stratonovich))
```

## RNG Control

Each `solve` or `init` call needs a random number generator (RNG) for
constructing noise processes. The RNG is resolved from the keyword arguments
in the following priority order:

1. **`rng` provided** — use the given `AbstractRNG` directly. Any `seed` kwarg
   is ignored. If a `TaskLocalRNG` is passed (i.e. `Random.default_rng()`), it
   is converted to a concrete `Xoshiro` seeded from one draw of the task-local
   stream, so the integrator never shares the global random stream.
2. **`seed` provided (nonzero)** — construct `Xoshiro(seed)`.
3. **Problem seed** (`prob.seed != 0`) — construct `Xoshiro(prob.seed)`.
4. **Neither** — generate a random seed and construct `Xoshiro` from it.

### Examples

```julia
using Random

# Reproducible results with an explicit RNG
rng = Xoshiro(42)
sol = solve(prob, EM(); dt = 0.01, rng)

# Same seed produces identical trajectories
rng2 = Xoshiro(42)
sol2 = solve(prob, EM(); dt = 0.01, rng = rng2)
sol.u == sol2.u  # true

# The older seed keyword still works (constructs Xoshiro(seed) internally,
# so reproducibility depends on the internal RNG type; prefer `rng` for
# guaranteed reproducibility across library versions)
sol = solve(prob, EM(); dt = 0.01, seed = UInt64(42))
```

### RNG ownership

The `rng` keyword controls **framework-constructed** randomness (noise processes
created internally by the solver, and the integrator's own RNG). If you supply
your own noise process via `SDEProblem(...; noise = my_W)`, that noise object's
internal RNG remains under your control and is **not** modified by the `rng`
keyword or by `set_rng!`.

### Integrator RNG interface

The integrator implements the SciMLBase RNG interface:

```julia
integ = init(prob, EM(); dt = 0.01, rng = Xoshiro(42))

SciMLBase.has_rng(integ)       # true
SciMLBase.get_rng(integ)       # the Xoshiro RNG
SciMLBase.set_rng!(integ, Xoshiro(99))  # replace with same-type RNG
```

`set_rng!` requires the new RNG to be the same concrete type as the current one.
For framework-constructed noise, it also syncs the noise process RNG. For
user-provided noise, only `integrator.rng` is updated.

The `reinit!` function also accepts an `rng` keyword:

```julia
reinit!(integ, u0; rng = Xoshiro(99))
```

## Performance Tips

 1. **Use appropriate solvers**: Match solver to problem type
 2. **In-place functions**: Use `f!(du,u,p,t)` for better performance
 3. **Tolerances**: Don't make tolerances unnecessarily strict
 4. **Static arrays**: Use `StaticArrays.jl` for small systems
 5. **GPU**: Use `CuArrays.jl` for large problems

## Common Pitfalls

 1. **Wrong noise type**: Ensure solver supports your noise structure
 2. **Stiffness**: Use appropriate stiff solvers for stiff problems
 3. **Commuting noise**: Use specialized solvers for better efficiency
 4. **High dimensions**: Consider weak convergence methods for Monte Carlo

## Integration with DifferentialEquations.jl

StochasticDiffEq.jl integrates with the broader ecosystem:

```julia
using DifferentialEquations

# Callbacks
condition(u, t, integrator) = u[1] - 0.5
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

sol = solve(prob, SOSRI(), callback = cb)

# Ensemble simulations
monte_prob = EnsembleProblem(prob)
sim = solve(monte_prob, SOSRI(), EnsembleThreads(), trajectories = 1000)
```
