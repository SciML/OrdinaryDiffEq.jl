# Implicit Methods for Stiff SDEs

When SDEs have stiff drift terms, explicit methods may require impractically small time steps for stability. Implicit methods treat the drift term implicitly while keeping the diffusion explicit, providing excellent stability properties.

## Recommended Stiff Methods

### SKenCarp - Stochastic KenCarp (Highly Recommended for Stiff Additive Noise)

```@docs
SKenCarp
```

## Basic Implicit Methods

### ImplicitEM - Implicit Euler-Maruyama

```@docs
ImplicitEM
```

### ImplicitEulerHeun - Implicit Euler-Heun (Stratonovich)

```@docs
ImplicitEulerHeun
```

### ImplicitRKMil - Implicit Runge-Kutta Milstein

```@docs
ImplicitRKMil
```

## Split-Step Implicit Methods

### ISSEM - Implicit Split-Step Euler-Maruyama

```@docs
ISSEM
```

### ISSEulerHeun - Implicit Split-Step Euler-Heun

```@docs
ISSEulerHeun
```

## Method Selection Guide

### Problem Classification:

 1. **Mildly stiff drift**: ImplicitEM, ImplicitRKMil
 2. **Stiff additive noise**: SKenCarp (highly recommended)
 3. **Fully stiff (including diffusion)**: ISSEM, ISSEulerHeun
 4. **Stratonovich problems**: ImplicitEulerHeun, ISSEulerHeun

### Performance Ranking for Stiff Problems:

 1. **SKenCarp** - Best for stiff problems with additive noise
 2. **ISSEM/ISSEulerHeun** - For fully stiff problems
 3. **ImplicitRKMil** - Higher order for mildly stiff problems
 4. **ImplicitEM** - Robust fallback option

## Understanding Stiffness in SDEs

**Drift Stiffness**: Large negative eigenvalues in the drift term f(u,t) requiring small time steps for explicit stability.

**Diffusion Stiffness**: Large coefficients in the diffusion term g(u,t) causing stability issues.

**Detection**: If explicit methods require very small dt or produce unstable solutions, try implicit methods.

## Configuration Options

All implicit methods share common configuration parameters:

```julia
# Linear solver options
ImplicitEM(linsolve = KrylovJL_GMRES())

# Nonlinear solver options  
ImplicitEM(nlsolve = NLNewton(max_iter = 20))

# Jacobian computation
ImplicitEM(autodiff = true, chunk_size = 4)

# Theta method parameter
ImplicitEM(theta = 0.5)  # Trapezoidal rule
```

## Symplectic Methods

For Hamiltonian SDEs, use symplectic variants:

```julia
SImplicitMidpoint()  # Symplectic implicit midpoint
STrapezoid()        # Symplectic trapezoidal rule
```

## Performance Tips

 1. **Jacobian**: Provide analytical Jacobian when possible
 2. **Linear solver**: Choose appropriate solver for problem structure
 3. **Preconditioning**: Use preconditioners for large systems
 4. **Theta parameter**: θ=0.5 often provides good accuracy/stability balance

## References

  - Standard implicit ODE methods adapted to SDEs
  - Milstein, G.N., "Numerical Integration of Stochastic Differential Equations"
