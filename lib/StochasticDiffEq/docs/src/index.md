# StochasticDiffEq.jl

StochasticDiffEq.jl is a component package in the DifferentialEquations ecosystem for solving stochastic differential equations (SDEs) and random ordinary differential equations (RODEs). It provides a comprehensive suite of high-performance numerical methods for stochastic problems.

## Installation

To install StochasticDiffEq.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("StochasticDiffEq")
```

## Quick Start

Here's a simple example of solving an SDE:

```julia
using StochasticDiffEq

# Define the drift and diffusion functions
function drift!(du, u, p, t)
    du[1] = 1.01 * u[1]
end

function diffusion!(du, u, p, t)
    du[1] = 0.87 * u[1]
end

# Initial condition and time span
u0 = [0.5]
tspan = (0.0, 1.0)

# Define the SDE problem
prob = SDEProblem(drift!, diffusion!, u0, tspan)

# Solve using the default algorithm
sol = solve(prob)
```

## Solver Categories

StochasticDiffEq.jl provides several categories of solvers optimized for different types of problems:

### Nonstiff Solvers

  - **Basic Methods**: Euler-Maruyama, Heun methods
  - **SRA/SRI Methods**: High-order adaptive methods (SOSRI, SOSRA)
  - **High Weak Order**: Methods optimized for weak convergence (DRI1)
  - **Commutative Noise**: Specialized methods for commuting noise terms

### Stiff Solvers

  - **Implicit Methods**: Drift-implicit methods for stiff problems
  - **Split-Step Methods**: Methods handling stiffness in diffusion
  - **Stabilized Methods**: SROCK-type methods for parabolic PDEs

### Jump-Diffusion

  - **Tau-Leaping**: Methods for jump-diffusion processes

## Recommended Methods

For most users, we recommend starting with these methods:

  - **General Purpose**: `SOSRI()` - Excellent for diagonal/scalar Itô SDEs
  - **Additive Noise**: `SOSRA()` - Optimal for problems with additive noise
  - **Stiff Problems**: `SKenCarp()` - Best for stiff problems with additive noise
  - **Commutative Noise**: `RKMilCommute()` - For multi-dimensional commutative noise
  - **High Efficiency**: `EM()` - When computational speed is most important

## Advanced Features

  - Adaptive time stepping with sophisticated error control
  - Support for all noise types (diagonal, non-diagonal, additive, scalar)
  - Both Itô and Stratonovich interpretations
  - Integration with the broader DifferentialEquations.jl ecosystem
  - GPU compatibility for high-performance computing
  - Extensive callback and event handling capabilities

See the individual solver pages for detailed information about each method's properties, when to use them, and their theoretical foundations.
