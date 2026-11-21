# Basic Nonstiff Methods

This page covers the fundamental explicit methods for solving SDEs. These methods are suitable for non-stiff problems and provide the foundation for more advanced algorithms.

## Euler-Maruyama Methods

### EM - Euler-Maruyama

```@docs
EM
```

### EulerHeun - Euler-Heun

```@docs
EulerHeun
```

### LambaEM - Adaptive Euler-Maruyama

```@docs
LambaEM
```

### LambaEulerHeun - Adaptive Euler-Heun

```@docs
LambaEulerHeun
```

## Milstein Methods

### RKMil - Runge-Kutta Milstein

```@docs
RKMil
```

## Split Methods

### SplitEM - Split Euler-Maruyama

```@docs
SplitEM
```

## When to Use Basic Methods

**Use EM when:**

  - Computational efficiency is most important
  - Problem is not stiff
  - Any noise type (most flexible)
  - Simple implementation needed

**Use EulerHeun when:**

  - Working in Stratonovich interpretation
  - Need slightly better accuracy than EM
  - Problem has non-commutative noise

**Use LambaEM/LambaEulerHeun when:**

  - Want adaptive time stepping with basic methods
  - Need error control but not high accuracy
  - Good balance of simplicity and adaptivity

**Use RKMil when:**

  - Higher accuracy than Euler methods
  - Problem has diagonal or scalar noise
  - Strong order 1.0 convergence required

These methods form the foundation of stochastic numerical analysis. While higher-order methods often provide better performance, the basic methods are essential for:

  - Initial testing and prototyping
  - Problems where simplicity is preferred
  - Educational purposes
  - Fallback options when advanced methods fail
