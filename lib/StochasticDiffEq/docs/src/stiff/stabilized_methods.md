# Stabilized Methods (SROCK Family)

Stabilized Runge-Kutta Chebyshev (SROCK) methods provide stability for mildly stiff problems through extended stability regions rather than implicit treatment. These methods are particularly effective for parabolic PDEs discretized by method of lines.

## SROCK Methods

### SROCK1 - First Order Stabilized Method

```@docs
SROCK1
```

### SROCK2, KomBurSROCK2, SROCKC2 - Second Order Methods

```@docs
SROCK2
```

```@docs
KomBurSROCK2
```

```@docs
SROCKC2
```

### SROCKEM - Stabilized Euler-Maruyama

```@docs
SROCKEM
```

### SKSROCK - Stabilized Method with Post-Processing

```@docs
SKSROCK
```

### TangXiaoSROCK2 - Alternative Second Order Method

```@docs
TangXiaoSROCK2
```

## When to Use Stabilized Methods

**Ideal for:**

  - Parabolic PDEs discretized by method of lines
  - Problems with moderate stiffness (not extremely stiff)
  - Large systems where implicit methods are expensive
  - When eigenvalue spectrum is primarily negative real

**Advantages over implicit methods:**

  - No linear system solves required
  - Better for large systems
  - Parallelizable
  - No Jacobian computation needed

**Disadvantages:**

  - Limited to moderate stiffness
  - May require eigenvalue estimation
  - Not effective for highly oscillatory problems

## Stability Regions

SROCK methods extend stability along the negative real axis:

  - **Standard explicit**: Stability region ~ [-2, 0]
  - **SROCK methods**: Stability region ~ [-s², 0] where s is the number of stages

The number of stages s is chosen based on estimated eigenvalues.

## Eigenvalue Estimation

Most SROCK methods accept an `eigen_est` parameter:

```julia
# Automatic estimation (default)
SROCK1()

# Manual estimation
SROCK1(eigen_est = -100.0)  # Largest eigenvalue magnitude

# Custom estimation function
SROCK1(eigen_est = my_estimator)
```

## Method Selection Guide

 1. **SROCK1**: Basic first-order method, most robust
 2. **SROCK2**: Higher accuracy, good general choice
 3. **SROCKEM**: When Euler-Maruyama structure is preferred
 4. **SKSROCK**: Advanced features, post-processing options
 5. **SROCKC2**: Conservative second-order variant

## Problem Suitability

**Well-suited:**

  - Reaction-diffusion equations
  - Heat equations with stochastic terms
  - Large sparse systems
  - Method of lines discretizations

**Not well-suited:**

  - Highly stiff problems (use implicit methods)
  - Problems with complex eigenvalue spectra
  - Small dense systems (overhead not justified)

## Configuration Tips

```julia
# For PDE problems
SROCK2(eigen_est = estimate_spectral_radius(A))

# For uncertain problems, start conservatively
SROCK1()  # Most robust

# For higher accuracy
SROCK2()  # Good balance
```

## Performance Considerations

  - Stage count increases with stiffness
  - Eigenvalue estimation cost
  - Memory requirements for internal stages
  - Better scalability than implicit methods

## References

  - Chebyshev methods for parabolic problems
  - ROCK methods for stiff ODEs
  - Stabilized explicit methods for PDEs
