# Type Traits

Many of the DiffEqBase abstract types have associated traits. These can be used
to check compatibility and apply separate code paths. For example, a parameter
estimation algorithm can set the default for using autodifferentiation by checking
if the algorithm is compatible with autodifferentiation.

Below are the abstract types along with the associated trait functions. These
are listed as:

```julia
f(x)
```

where `f` is the trait function and `x` is the any type which subtypes the abstract
type.

### AbstractODEProblem

  - `isinplace` : Returns true if the problem uses in-place functions

### AbstractRODEProblem

  - `is_diagonal_noise` : Returns true if the noise is diagonal.

### DEAlgorithm

  - `isautodifferentiable` : Returns true if the algorithm is autodifferentiable.
