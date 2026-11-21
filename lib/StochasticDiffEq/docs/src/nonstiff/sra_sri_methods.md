# SRA/SRI Methods - Stochastic Runge-Kutta

The SRA (Stochastic Runge-Kutta for Additive noise) and SRI (Stochastic Runge-Kutta for Itô) methods provide high-order adaptive solvers for different noise structures. These are among the most effective methods for their respective problem classes.

## Recommended Methods

### SOSRI - Stability-Optimized SRI (Recommended)

```@docs
SOSRI
```

### SOSRA - Stability-Optimized SRA (Optimal for Additive Noise)

```@docs
SOSRA
```

## Alternative SRI Methods

### SRIW1 - SRI Weak Order 2

```@docs
SRIW1
```

### SRIW2 - SRI Weak Order 3

```@docs
SRIW2
```

### SOSRI2 - Alternative Stability-Optimized SRI

```@docs
SOSRI2
```

## Alternative SRA Methods

### SRA1 - Original SRA Method

```@docs
SRA1
```

### SRA2 - SRA Method Version 2

```@docs
SRA2
```

### SRA3 - SRA Method with Weak Order 3

```@docs
SRA3
```

### SOSRA2 - Alternative Stability-Optimized SRA

```@docs
SOSRA2
```

## Configurable Methods

### SRA - Configurable SRA with Custom Tableaux

```@docs
SRA
```

### SRI - Configurable SRI with Custom Tableaux

```@docs
SRI
```

## Method Selection Guide

### For Diagonal/Scalar Noise:

 1. **First choice: SOSRI** - Best overall performance and stability
 2. **Alternative: SRIW1** - Standard SRI method
 3. **High weak order: SRIW2** - When weak order 3 is needed

### For Additive Noise:

 1. **First choice: SOSRA** - Optimal for additive noise structure
 2. **Alternative: SRA1** - Standard SRA method
 3. **High weak order: SRA3** - When weak order 3 is needed

### Performance Characteristics:

  - **SOSRI/SOSRA**: Stability-optimized, robust to high tolerances
  - **SRIWx/SRAx**: Standard methods with proven theoretical properties
  - **SRA/SRI**: Allow custom tableaux for specialized applications

## Theoretical Foundation

The SRA and SRI methods are based on stochastic Runge-Kutta theory:

**SRA Methods** exploit the additive noise structure:

```
du = f(u,t)dt + σ(t)dW
```

Where the diffusion σ doesn't depend on the solution u.

**SRI Methods** handle the general diagonal case:

```
du = f(u,t)dt + g(u,t)dW
```

Where each component has independent noise.

Both method families achieve:

  - Strong order 1.5 convergence
  - Weak order 2.0 or higher
  - Adaptive time stepping with embedded error estimation
  - A-stable or L-stable properties (for optimized versions)

## References

  - Rößler A., "Runge–Kutta Methods for the Strong Approximation of Solutions of Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952
