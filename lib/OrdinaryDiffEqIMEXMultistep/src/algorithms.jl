# IMEX Multistep methods

"""
    CNAB2(; autodiff = AutoForwardDiff(), concrete_jac = nothing, linsolve = nothing,
        nlsolve = NLNewton(), extrapolant = :linear)

Second-order Crank-Nicolson Adams-Bashforth IMEX multistep method. The first
component of a `SplitODEProblem` is treated implicitly with Crank-Nicolson, while
the second component is treated explicitly with Adams-Bashforth.

`CNAB2` uses fixed time steps. Supply `dt` when solving and do not enable adaptive
time stepping.

# Fields

  - `linsolve`: optional linear solver configuration used by the implicit solve.
  - `nlsolve`: nonlinear solver algorithm used for the implicit equation.
  - `extrapolant`: strategy used to initialize the nonlinear solve.
  - `autodiff`: automatic-differentiation backend used to construct Jacobians.
  - `concrete_jac`: controls whether a concrete Jacobian is cached when supported.

# Keywords

  - `autodiff = AutoForwardDiff()`: ADTypes backend used for Jacobian construction.
  - `concrete_jac = nothing`: retain the default Jacobian-concretization policy. Set a
    boolean value to request or disable a concrete Jacobian where supported.
  - `linsolve = nothing`: linear solver algorithm or configuration. `nothing` selects
    the OrdinaryDiffEq default.
  - `nlsolve = NLNewton()`: nonlinear solver used for the implicit equation.
  - `extrapolant = :linear`: initial-guess strategy for the nonlinear solve.

# Returns

An algorithm object that can solve split ODE problems with fixed time steps.

# Examples

```julia
using OrdinaryDiffEqIMEXMultistep: CNAB2
using SciMLBase: SplitODEProblem, solve

f_implicit(u, p, t) = -9u
f_explicit(u, p, t) = -u
prob = SplitODEProblem(f_implicit, f_explicit, 1.0, (0.0, 1.0))
sol = solve(prob, CNAB2(), dt = 0.05)
```

# References

Y. He and J. Li, "Numerical implementation of the Crank-Nicolson/Adams-Bashforth
scheme for the time-dependent Navier-Stokes equations," *International Journal for
Numerical Methods in Fluids* 62 (2010), pp. 647-659.
"""
struct CNAB2{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end

function CNAB2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear
    )
    autodiff = _fixup_ad(autodiff)

    return CNAB2(
        linsolve,
        nlsolve,
        extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end

"""
    CNLF2(; autodiff = AutoForwardDiff(), concrete_jac = nothing, linsolve = nothing,
        nlsolve = NLNewton(), extrapolant = :linear)

Second-order Crank-Nicolson Leapfrog IMEX multistep method. The first component of
a `SplitODEProblem` is treated implicitly with Crank-Nicolson, while the second
component is advanced with an explicit leapfrog formula.

`CNLF2` uses fixed time steps. Supply `dt` when solving and do not enable adaptive
time stepping.

# Fields

  - `linsolve`: optional linear solver configuration used by the implicit solve.
  - `nlsolve`: nonlinear solver algorithm used for the implicit equation.
  - `extrapolant`: strategy used to initialize the nonlinear solve.
  - `autodiff`: automatic-differentiation backend used to construct Jacobians.
  - `concrete_jac`: controls whether a concrete Jacobian is cached when supported.

# Keywords

  - `autodiff = AutoForwardDiff()`: ADTypes backend used for Jacobian construction.
  - `concrete_jac = nothing`: retain the default Jacobian-concretization policy. Set a
    boolean value to request or disable a concrete Jacobian where supported.
  - `linsolve = nothing`: linear solver algorithm or configuration. `nothing` selects
    the OrdinaryDiffEq default.
  - `nlsolve = NLNewton()`: nonlinear solver used for the implicit equation.
  - `extrapolant = :linear`: initial-guess strategy for the nonlinear solve.

# Returns

An algorithm object that can solve split ODE problems with fixed time steps.

# Examples

```julia
using OrdinaryDiffEqIMEXMultistep: CNLF2
using SciMLBase: SplitODEProblem, solve

f_implicit(u, p, t) = -9u
f_explicit(u, p, t) = -u
prob = SplitODEProblem(f_implicit, f_explicit, 1.0, (0.0, 1.0))
sol = solve(prob, CNLF2(), dt = 0.05)
```

# References

N. Jiang, M. Kubacki, W. Layton, M. Moraiti, and H. Tran, "A Crank-Nicolson
Leapfrog stabilization: Unconditional stability and two applications," *Journal of
Computational and Applied Mathematics* 281 (2015), pp. 263-276.
"""
struct CNLF2{AD, F, F2, CJ} <:
    OrdinaryDiffEqNewtonAlgorithm
    linsolve::F
    nlsolve::F2
    extrapolant::Symbol
    autodiff::AD
    concrete_jac::CJ
end
function CNLF2(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :linear
    )
    autodiff = _fixup_ad(autodiff)

    return CNLF2(
        linsolve,
        nlsolve,
        extrapolant,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end
