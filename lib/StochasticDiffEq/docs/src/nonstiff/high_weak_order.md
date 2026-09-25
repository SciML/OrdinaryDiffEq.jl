# High Weak Order Methods

These methods are specifically designed for problems where weak convergence is more important than strong convergence. They are optimal for Monte Carlo simulations, computing expectations, moments, and other statistical properties of solutions.

## Recommended High Weak Order Methods

### DRI1 - Debrabant-Rößler Method (Weak Order 2)

```@docs
DRI1
```

### DRI1NM - Debrabant-Rößler for Non-mixing Diagonal Problems

```@docs
DRI1NM
```

## Other Weak Order 2 Methods

### RI1, RI3, RI5, RI6 - Rößler Methods

```@docs
RI1
```

```@docs
RI3
```

```@docs
RI5
```

```@docs
RI6
```

### RDI Methods - Alternative Weak Order 2

```@docs
RDI2WM
```

```@docs
RDI3WM
```

```@docs
RDI4WM
```

### W2Ito1 - Efficient Weak Order 2

```@docs
W2Ito1
```

## Fixed Step Methods

### PL1WM, PL1WMA - Platen Methods

```@docs
PL1WM
```

```@docs
PL1WMA
```

### Stratonovich Methods

#### RS1, RS2 - Rößler Stratonovich Methods

```@docs
RS1
```

```@docs
RS2
```

#### NON, NON2 - Non-commutative Stratonovich

```@docs
NON
```

```@docs
NON2
```

#### COM - Commutative Stratonovich

```@docs
COM
```

## Specialized Methods

### SIEA, SMEA, SIEB, SMEB - Tocino-Vigo-Aguiar Methods

```@docs
SIEA
```

```@docs
SMEA
```

```@docs
SIEB
```

```@docs
SMEB
```

## Weak vs Strong Convergence

**Strong Convergence**: Measures pathwise error E[|X(T) - X_h(T)|^p]
**Weak Convergence**: Measures error in expectations E[f(X(T))] - E[f(X_h(T))]

### When to Use Weak Order Methods:

  - Monte Carlo simulations
  - Computing expectations and moments
  - Statistical analysis of SDEs
  - When pathwise accuracy is not critical
  - Large ensemble simulations

### Advantages:

  - Often more efficient for statistical quantities
  - Can use larger time steps while maintaining weak accuracy
  - Optimized error constants for better practical performance

## Method Selection Guide

 1. **General purpose weak order 2**: DRI1
 2. **Non-mixing diagonal**: DRI1NM
 3. **Fixed step**: PL1WM, RS1/RS2
 4. **Stratonovich**: RS1, RS2, NON, NON2
 5. **Specialized applications**: RI methods, RDI methods

## References

  - Debrabant, K. and Rößler A., "Families of efficient second order Runge–Kutta methods for the weak approximation of Itô stochastic differential equations"
  - Rößler A., "Second Order Runge–Kutta Methods for Itô Stochastic Differential Equations"
