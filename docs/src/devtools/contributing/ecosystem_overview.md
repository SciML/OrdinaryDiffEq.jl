# Ecosystem Overview

So you're looking to help out DifferentialEquations.jl? We'd be happy to have
your help. It is recommended you first discuss with some of the developers
[on the Gitter channel](https://gitter.im/JuliaDiffEq/Lobby)
to make sure that you're up-to-date with current developments.

### The Common Interface

The DiffEq ecosystem is built around the common interface. This is the interface
for the solvers:

```julia
__solve(prob, alg; kwargs...)
__init(prob, alg; kwargs...)
```

and the standard methods for dealing with solutions. A higher level `solve` and
`init` is given by DiffEqBase.jl for functional and distributional inputs.
Users build problem
types for solvers to act on, and add-on components which use the solution types
for higher-level analysis like parameter estimation and sensitivity analysis.

One can add components at any of these levels to improve the functionality of
the system as a whole.

### Organizational Setup

JuliaDiffEq is setup in a distributed manner to allow developers to retain authoritative
control and licensing for their own packages/algorithms, yet contribute to the
greater ecosystem. This gives a way for researchers to target a wide audience
of users, but not have to fully contribute to public packages or be restricted in
licensing. At the center of the ecosystem is DiffEqBase which holds
the `Problem`, `Solution`, and `Algorithm` types (the algorithms are defined in
DiffEqBase to be accessible by the `default_algorithm` function. One can opt out
of this). Then there's the component solvers, which includes the `*DiffEq` packages
(OrdinaryDiffEq, StochasticDiffEq, etc.) which implement different methods for
`solve`. Then there are the add-on packages, such as the `DiffEq*` packages (DiffEqParamEstim,
DiffEqDevTools) which add functionality to the `Problem`+`solve` setup. Lastly,
there's DifferentialEquations.jl which is a metapackage which holds all of these
pieces together as one cohesive unit.

If one wants their package to officially join the ecosystem, it will need to be
moved to the JuliaDiffEq organization so that maintenance can occur (but the
core algorithms will only be managed by the package owners themselves). The `Algorithm`
types can then be moved to `DiffEqBase`, and after testing the package will be added to
the list packages exported by DifferentialEquations.jl and the corresponding documentation.
