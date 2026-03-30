# SciML Scientific Machine Learning Developer Documentation

This is the developer documentation and Contributor's Guide for the SciML ecosystem. It explains the common
interface and some the package internals to help developers contribute.

If you have any questions, or just want to chat about solvers/using the package, please feel free to use the [Gitter channel](https://gitter.im/JuliaDiffEq/Lobby). For bug reports, feature requests, etc., please submit an issue.

### Overview

The SciML ecosystem is built around the common interface. The common interface
is a type-based interface where users define problems as a type, and solvers
plug into the ecosystem by defining an algorithm to give a new dispatch to

```julia
__solve(prob, alg; kwargs...)
__init(prob, alg; kwargs...)
```

There is a top level `solve` and `init` function which is
in DiffEqBase.jl that handles distribution and function input
(along with extra warnings) before sending the problems
to the packages.

There is then an ecosystem of add-on components which use the
common solver interface to add analysis tools for differential
equations.

### Contributing to the Ecosystem

There are many ways to help the ecosystem. One way you can contribute is to give
pull requests (PRs) to existing packages. Another way to contribute is to add your
own package to the ecosystem. Adding your own package to the ecosystem allows
you to keep executive control and licensing over your methods,
but allows users of DifferentialEquations.jl to use your methods via the common
interface, and makes your package compatible with the add-on tools (sensitivity
analysis, parameter estimation, etc). Note that, in order for the method to be
used as a default, one is required to move their
package to the SciML organization so that way common maintenance (such
as fixing deprecation warnings, updating tests to newer versions, and emergency
fixes / disabling) can be allowed by SciML members. However, the lead developer
of the package maintains administrative control, and thus any change to the core
algorithms by other SciML members will only be given through PRs.

Even if you don't have the time to contribute new solver algorithms or add-on tools,
there's always ways to help! Improved plot recipes and new series recipes are
always nice to add more default plots. It is always helpful to have benchmarks
between different algorithms to see "which is best". Adding examples IJulia
notebooks to DiffEqTutorials.jl is a good way to share knowledge about DifferentialEquations.jl.
Also, please feel free to comb through the solvers and look for ways to make them
more efficient. Lastly, the documentation could always use improvements. If you
have any questions on how to help, just ask them in the Gitter!

### Code of Conduct

All contributors must adhere to the [NumFOCUS Code of Conduct](https://numfocus.org/code-of-conduct).
Treat everyone with respect. Failure to comply will result in individuals being banned from the community.

### Contributor Guide

```@contents
Pages = [
  "contributing/ecosystem_overview.md",
  "contributing/adding_packages.md",
  "contributing/adding_algorithms.md",
  "contributing/defining_problems.md",
  "contributing/diffeq_internals.md",
  "contributing/type_traits.md"
]
Depth = 2
```

### Algorithm Development Tools

The following algorithm development tools are provided by DiffEqDevTools.jl

```@contents
Pages = [
  "alg_dev/test_problems.md",
  "alg_dev/convergence.md",
  "alg_dev/benchmarks.md"
]
Depth = 2
```

### Internal Documentation

```@contents
Pages = [
  "internals/fem_tools.md",
  "internals/notes_on_algorithms.md",
  "internals/tableaus.md"
]
Depth = 2
```

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@raw html
You can also download the 
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
       "/assets/Manifest.toml"
```

```@raw html
">manifest</a> file and the
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link = "https://github.com/SciML/" * name * ".jl/tree/gh-pages/v" * version *
       "/assets/Project.toml"
```

```@raw html
">project</a> file.
```
