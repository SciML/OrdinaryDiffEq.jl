# GlobalDiffEq

Differential equation solvers and solver wrappers with global (accumulated)
error estimates and global error control for the
[DifferentialEquations.jl common solver interface](https://docs.sciml.ai/DiffEqDocs/stable/).
Formerly the standalone [GlobalDiffEq.jl](https://github.com/SciML/GlobalDiffEq.jl)
repository, now maintained as a sublibrary of OrdinaryDiffEq.jl.

Standard adaptive ODE solvers control the *local* error of each step, which
does not bound the error accumulated over the whole integration. This package
provides:

  - `GlobalRichardson`: global Richardson extrapolation of whole solves of any
    fixed-step method, interpreting `abstol`/`reltol` as global tolerances.

See the [OrdinaryDiffEq.jl documentation](https://docs.sciml.ai/OrdinaryDiffEq/stable/)
for the full API reference.
