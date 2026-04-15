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
using Pkg
pathtolibrary = Pkg.pkgdir(OrdinaryDiffEq)
sublibs = joinpath.(pathtolibrary, "lib", readdir(joinpath(pathtolibrary, "lib")))
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

## Test Group Structure

Each sublibrary declares its CI jobs through `test/test_groups.toml`. Groups listed there
(e.g. `Core`, `QA`, `GPU`, `ModelingToolkit`) become matrix entries dispatched by
`.github/scripts/compute_affected_sublibraries.jl`. The group name is passed to the
sublibrary's `test/runtests.jl` via `ENV["ODEDIFFEQ_TEST_GROUP"]`.

### Per-group test environments

A group that needs dependencies beyond the sublibrary's main `[targets].test` list should
use its own isolated environment under `test/<group>/Project.toml`. `test/runtests.jl`
activates that environment before running the group's tests:

```julia
function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
end

if (TEST_GROUP == "QA" || TEST_GROUP == "ALL") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "Allocation Tests" include("qa/allocation_tests.jl")
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end
```

This pattern is the standard for the `QA` (JET + Aqua + AllocCheck) and `GPU` groups.
Use it for any new group whose dependencies shouldn't leak into the main test env or
into reverse-dependency resolution.

Two conventions that keep this working:

 1. Heavy test-only deps (JET, Aqua, AllocCheck, CUDA, MTK, …) live **only** in
    `test/<group>/Project.toml`, not in the sublibrary's `[extras]`. Do not rely on
    `Pkg.add` from inside a test file.
 2. `test/<group>/Project.toml` uses relative `[sources]` paths (e.g.
    `path = "../.."` for the sublibrary itself, `path = "../../../OrdinaryDiffEqCore"`
    for internal deps) so that the isolated env resolves against the in-tree sources.
