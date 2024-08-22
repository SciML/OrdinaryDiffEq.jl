  - This repository follows the [SciMLStyle](https://github.com/SciML/SciMLStyle) and the SciML [ColPrac](https://github.com/SciML/ColPrac).
  - Please run `using JuliaFormatter, OrdinaryDiffEq; format(joinpath(dirname(pathof(OrdinaryDiffEq)), ".."))` before committing.
  - Add tests for any new features.
  - Additional help on contributing to the numerical solvers can be found at https://docs.sciml.ai/DiffEqDevDocs/stable/

## Developing Locally

OrdinaryDiffEq is a very large package and thus it uses a sublibrary approach to keep down
the total number of dependencies per solver. As a consequence, it requires a bit of special
handling compared to some other Julia packages. When running the full test suite, it's
recommended that one has dev'd all of the relevant packages. This can be done via:

```julia
pathtolibrary = Pkg.pkgdir(OrdinaryDiffEq)
sublibs = string.((pathtolibrary,), readdir(pathtolibrary))
Pkg.develop(map(name -> Pkg.PackageSpec.(; path = name), sublibs));
```

and then running `Pkg.test("OrdinaryDiffEq")` will run the global tests locally. Each of the
per-solver subtests, which includes convergence tests and other per-library handling,
is then ran by using a specific solver set, such as
`Pkg.test("OrdinaryDiffEqAdamsBashforthMoulton")`

## Dependency Structure

There is a tree dependency structure to the sublibraries of OrdinaryDiffEq. The core stepper
lives in OrdinaryDiffEqCore.jl. Explicit methods only require the core. Then on top of the
core is OrdinaryDiffEqDifferentiation.jl which defines the utilities used for differentiation
and Jacobian construction. This is used directly by methods like Rosenbrock and exponential
integrators, which use Jacobians but not implicit solvers. Finally there's
OrdinaryDiffEqNonlinearSolve, which is the largest set of dependencies and pulls in the
NonlinearSolve.jl stack for handling implicit equations, and is thus the core dependency
of the implicit integrators. OrdinaryDiffEqNonlinearSolve relies on OrdinaryDiffEqDifferentiation
which relies on OrdinaryDiffEqCore.
