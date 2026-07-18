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
  - `GLEE23`, `GLEE24`, `GLEE35`: explicit general linear methods with built-in
    global error estimation (Constantinescu 2016), which propagate the solution
    together with an asymptotically correct estimate of its global error.
  - `MM5GEE`: the Makazaga-Murua (2003) Dormand-Prince-based order-5 scheme
    with a cheap built-in global error estimate.
  - `GlobalErrorTransport`: global error estimation and `gtol`-based control by
    integrating the linearized error-transport equation driven by the
    dense-output defect (Shampine 1986, Berzins 1988, Lang-Verwer 2007).
  - `GlobalDefectCorrection`: global error estimation and control by solving
    for the correction (Zadunaisky 1976, Dormand-Duckers-Prince 1984/1989),
    requiring no Jacobian.
  - `GlobalAdjoint`: adjoint-based a posteriori endpoint error estimation and
    control (Cao and Petzold 2004), available as a package extension when
    SciMLSensitivity and QuadGK are loaded.

See the [OrdinaryDiffEq.jl documentation](https://docs.sciml.ai/OrdinaryDiffEq/stable/)
for the full API reference.
