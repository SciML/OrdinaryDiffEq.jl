# Split-Step Methods for Fully Stiff Problems

When both drift and diffusion terms are stiff, split-step methods treat both parts implicitly. These methods can handle the most challenging stiff SDE problems.

## Split-Step Implicit Methods

### ISSEM - Implicit Split-Step Euler-Maruyama

```@docs
ISSEM
```

### ISSEulerHeun - Implicit Split-Step Euler-Heun

```@docs
ISSEulerHeun
```

## When to Use Split-Step Methods

**Use ISSEM when:**

  - Both drift and diffusion are stiff
  - Explicit methods fail even with small time steps
  - Standard implicit methods (ImplicitEM) are insufficient
  - Working with Itô interpretation

**Use ISSEulerHeun when:**

  - Both drift and diffusion are stiff
  - Working with Stratonovich interpretation
  - Need implicit treatment of diffusion term

## Understanding Full Stiffness

**Drift Stiffness**: Large negative eigenvalues in f(u,t)
**Diffusion Stiffness**: Large coefficients in g(u,t) causing instability

**Detection Signs:**

  - ImplicitEM still requires very small time steps
  - Solutions become unstable despite implicit drift treatment
  - Large diffusion coefficients cause numerical artifacts

## Algorithm Structure

Split-step methods solve:

```
du = f(u,t)dt + g(u,t)dW
```

By treating both terms implicitly through operator splitting or fully implicit schemes.

## Performance Considerations

  - More expensive per step than drift-only implicit methods
  - May require smaller time steps than expected
  - Jacobian computations for both drift and diffusion
  - Nonlinear solve complexity increases

## Configuration

Same options as other implicit methods:

```julia
ISSEM(
    linsolve = KrylovJL_GMRES(),
    nlsolve = NLNewton(),
    theta = 1.0,
    autodiff = true
)
```

## Alternative Approaches

If split-step methods are too expensive:

 1. Try stabilized explicit methods (SROCK family)
 2. Consider method of lines for PDE problems
 3. Use shorter time intervals with restarts
 4. Reformulate problem to reduce stiffness

## References

  - Implicit-explicit splitting schemes for SDEs
  - Stochastic operator splitting methods
