# Controlling Solver Verbosity

OrdinaryDiffEq.jl provides fine-grained control over diagnostic messages, warnings, and errors
through the `verbose` keyword argument for `solve`. The verbosity system allows you to control what
information is displayed during the solve process. See [SciMLLogging.jl](https://docs.sciml.ai/SciMLLogging/dev/) for more details. 

## Migrating from `verbose::Bool`

In v6, `verbose` was a `Bool`. In v7 it is a [`DEVerbosity`](@ref) object (or a
`SciMLLogging` preset, which is converted to one). Both `DEVerbosity` and the
`SciMLLogging` module are exported by OrdinaryDiffEq, so no additional package is
needed:

```julia
using OrdinaryDiffEq

prob = ODEProblem((du, u, p, t) -> (du[1] = -u[1]; nothing), [1.0], (0.0, 1.0))

solve(prob, Tsit5(); verbose = DEVerbosity(SciMLLogging.None()))  # was verbose = false
solve(prob, Tsit5(); verbose = DEVerbosity())                     # was verbose = true
```

From a solver package that does not re-export them, add `using DiffEqBase, SciMLLogging`.

```@docs
DiffEqBase.DEVerbosity
```
