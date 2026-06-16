# Commutative Noise Methods

When multiple noise sources satisfy commutativity conditions, specialized methods can achieve higher accuracy and efficiency compared to general methods. These methods avoid expensive Lévy area computations while maintaining high order.

## Recommended Commutative Noise Methods

### RKMilCommute - Runge-Kutta Milstein for Commutative Noise

```@docs
RKMilCommute
```

### RKMilGeneral - General Milstein for Non-commutative Noise

```@docs
RKMilGeneral
```

## Three-Stage Milstein Methods

### WangLi3SMil Family - Fixed Step Milstein Methods

```@docs
WangLi3SMil_A
```

```@docs
WangLi3SMil_B
```

```@docs
WangLi3SMil_C
```

```@docs
WangLi3SMil_D
```

```@docs
WangLi3SMil_E
```

```@docs
WangLi3SMil_F
```

## Understanding Commutative Noise

### Commutativity Condition

Noise terms g₁, g₂, ..., gₘ are **commutative** if:

```
[gᵢ, gⱼ] = gᵢ(∂gⱼ/∂x) - gⱼ(∂gᵢ/∂x) = 0
```

for all pairs (i,j).

### Examples of Commutative Noise:

 1. **Additive noise**: g(u,t) = σ(t) (independent of u)
 2. **Scalar multiplicative**: g(u,t) = σ(t)u (same u dependence)
 3. **Diagonal with same function**: gᵢ(u,t) = σᵢ(t)h(u)

### Examples of Non-commutative Noise:

 1. **Different multiplicative terms**: g₁ = σ₁u₁, g₂ = σ₂u₂
 2. **Cross-coupling**: g₁ = σ₁₁u₁ + σ₁₂u₂, g₂ = σ₂₁u₁ + σ₂₂u₂

## Method Selection Guide

### For Commutative Noise:

 1. **RKMilCommute** - Adaptive, excellent general choice
 2. **WangLi3SMil methods** - Fixed step, when dt is predetermined

### For Non-commutative Noise:

 1. **RKMilGeneral** - Handles general case with Lévy area
 2. **Fall back to SRI/SRA methods** - More robust for complex noise

### For Uncertain Commutativity:

  - Test with RKMilCommute first
  - If results seem incorrect, switch to RKMilGeneral or SRI methods

## Computational Advantages

**Commutative case:**

  - No Lévy area computation needed
  - Simpler stochastic integrals
  - Higher efficiency per step
  - Better scalability to high dimensions

**Non-commutative case:**

  - Requires Lévy area approximation
  - More expensive per step
  - Uses specialized algorithms (LevyArea.jl integration)

## Performance Tips

 1. **Verify commutativity** before using specialized methods
 2. **Use appropriate interpretation** (Itô vs Stratonovich)
 3. **Consider problem dimension** - benefits increase with system size
 4. **Test accuracy** - commutative methods may be more sensitive

## Iterated Integrals and Lévy Area

### For Commutative Noise (IICommutative):

Only simple stochastic integrals ∫₀ᵗ dWₛ are needed.

### For Non-commutative Noise (IILevyArea):

Requires double integrals ∫₀ᵗ ∫₀ˢ dWᵤdWₛ (Lévy area).

RKMilGeneral automatically chooses optimal Lévy area algorithms via LevyArea.jl.

## References

  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential Equations"
  - Wang and Li, "Three-stage stochastic Runge-Kutta methods for stochastic differential equations"
  - Kastner, F. and Rößler, A., "LevyArea.jl" for Lévy area computation
