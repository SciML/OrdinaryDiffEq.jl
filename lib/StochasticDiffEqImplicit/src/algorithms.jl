using ADTypes: AutoForwardDiff
using OrdinaryDiffEqCore: _fixup_ad

"""
    ImplicitEM(;
        theta = 1, symplectic = false, linsolve = nothing,
        nlsolve = NLNewton(), extrapolant = :constant,
        new_jac_conv_bound = 1.0e-3,
        autodiff = AutoForwardDiff(), concrete_jac = nothing
    )

**ImplicitEM: Drift-Implicit Euler-Maruyama (Stiff)**

Drift-implicit theta-method version of Euler-Maruyama for Itô SDEs. The drift is
integrated with the theta method while the diffusion is kept explicit, which gives
stability for problems with a stiff drift term.

## Method Properties

  - **Strong Order**: 0.5
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive
  - **Noise types**: All (diagonal, non-diagonal, scalar)
  - **SDE interpretation**: Itô

## Mathematical Formulation

```math
u_{n+1} = u_n + dt (θ f(u_{n+1}, t_{n+1}) + (1 - θ) f(u_n, t_n)) + g(u_n, t_n) ΔW_n
```

## When to Use

  - SDEs with a stiff drift term and non-stiff diffusion
  - As a robust fallback when explicit methods require impractically small steps
  - Problems where only a low-order method is needed but stability is essential

If the diffusion term is also stiff, prefer the split-step method [`ISSEM`](@ref).
For Stratonovich problems, use [`ImplicitEulerHeun`](@ref).

## Keyword Arguments

  - `theta`: Implicitness parameter of the drift discretization (default: `1`, i.e.
    drift-implicit Euler). `theta = 1/2` gives the trapezoidal rule ([`STrapezoid`](@ref)).
  - `symplectic`: If `true`, uses the symplectic midpoint discretization
    (forces `theta = 1/2`). See [`SImplicitMidpoint`](@ref).
  - `linsolve`: LinearSolve.jl algorithm used for the Newton linear systems
    (default: `nothing`, i.e. the LinearSolve.jl default).
  - `nlsolve`: Nonlinear solver used for the implicit stage (default: `NLNewton()`).
  - `extrapolant`: Initial guess method for the nonlinear solve (default: `:constant`).
  - `new_jac_conv_bound`: Convergence bound below which the Jacobian is reused
    instead of recomputed (default: `1e-3`).
  - `autodiff`: ADTypes.jl backend used for the Jacobian (default: `AutoForwardDiff()`).
  - `concrete_jac`: Whether to build a concrete Jacobian matrix even when a
    matrix-free linear solver is used (default: `nothing`, i.e. decided automatically).

## References

  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential
    Equations", Springer (1992).
  - Milstein, G.N., "Numerical Integration of Stochastic Differential Equations",
    Kluwer (1995).
"""
struct ImplicitEM{AD, F, F2, T2, T3, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::F2
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T3
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
end
function ImplicitEM(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitEM(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant, new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac),
    )
end

"""
    STrapezoid(; kwargs...)

**STrapezoid: Stochastic Trapezoidal Rule (Stiff)**

Convenience constructor for the trapezoidal-rule member of the [`ImplicitEM`](@ref)
family, i.e. `ImplicitEM(theta = 1/2)`. The drift is discretized with the trapezoidal
rule and the diffusion is kept explicit.

## Method Properties

  - **Strong Order**: 0.5
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive
  - **Noise types**: All (diagonal, non-diagonal, scalar)
  - **SDE interpretation**: Itô

## When to Use

  - Stiff Itô SDEs where the second-order accurate (in the deterministic part)
    trapezoidal drift discretization is preferred over backward Euler
  - Problems where backward Euler's numerical damping is undesirable

Note that the trapezoidal rule is only A-stable, not L-stable, so for very stiff
drifts the default `ImplicitEM()` (`theta = 1`) may be more robust.

## Keyword Arguments

All keyword arguments of [`ImplicitEM`](@ref) are accepted and forwarded, except
that `theta` is fixed to `1/2`.
"""
STrapezoid(; kwargs...) = ImplicitEM(; theta = 1 / 2, kwargs...)

"""
    SImplicitMidpoint(; kwargs...)

**SImplicitMidpoint: Stochastic Implicit Midpoint Rule (Stiff, Symplectic)**

Convenience constructor for the symplectic midpoint member of the [`ImplicitEM`](@ref)
family, i.e. `ImplicitEM(theta = 1/2, symplectic = true)`. The midpoint discretization
of the drift makes the deterministic part of the integrator symplectic, so it preserves
quadratic invariants and has good long-time energy behavior.

## Method Properties

  - **Strong Order**: 0.5
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive
  - **Noise types**: All (diagonal, non-diagonal, scalar)
  - **SDE interpretation**: Itô
  - **Symplectic**: Yes (deterministic part)

## When to Use

  - Stiff Hamiltonian or otherwise structure-preserving SDEs
  - Long-time integrations where drift in conserved quantities must be avoided
  - When the non-damping behavior of the midpoint rule is preferable to the
    dissipativity of backward Euler

## Keyword Arguments

All keyword arguments of [`ImplicitEM`](@ref) are accepted and forwarded, except
that `theta` is fixed to `1/2` and `symplectic` is fixed to `true`.

## References

  - Milstein, G.N., Repin, Yu.M., Tretyakov, M.V., "Numerical Methods for Stochastic
    Systems Preserving Symplectic Structure", SIAM J. Numer. Anal., 40 (4),
    pp. 1583–1604 (2002). DOI: 10.1137/S0036142901395588.
"""
SImplicitMidpoint(; kwargs...) = ImplicitEM(; theta = 1 / 2, symplectic = true, kwargs...)

"""
    ImplicitEulerHeun(;
        theta = 1, symplectic = false, linsolve = nothing,
        nlsolve = NLNewton(), extrapolant = :constant,
        new_jac_conv_bound = 1.0e-3,
        autodiff = AutoForwardDiff(), concrete_jac = nothing
    )

**ImplicitEulerHeun: Drift-Implicit Euler-Heun (Stiff, Stratonovich)**

Drift-implicit theta-method version of the Euler-Heun scheme, the Stratonovich
counterpart of [`ImplicitEM`](@ref). The drift is treated implicitly while the
Stratonovich diffusion is handled with the explicit Heun predictor-corrector.

## Method Properties

  - **Strong Order**: 0.5
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive
  - **Noise types**: All (diagonal, non-diagonal, scalar)
  - **SDE interpretation**: Stratonovich

## When to Use

  - Stratonovich SDEs with a stiff drift term
  - Physical models where the noise is a smooth approximation limit (Wong-Zakai),
    which is naturally Stratonovich
  - As the Stratonovich analogue of `ImplicitEM()`

If the diffusion is also stiff, prefer the split-step method [`ISSEulerHeun`](@ref).

## Keyword Arguments

  - `theta`: Implicitness parameter of the drift discretization (default: `1`, i.e.
    drift-implicit Euler). `theta = 1/2` gives the trapezoidal rule.
  - `symplectic`: If `true`, uses the symplectic midpoint discretization
    (forces `theta = 1/2`).
  - `linsolve`: LinearSolve.jl algorithm used for the Newton linear systems
    (default: `nothing`, i.e. the LinearSolve.jl default).
  - `nlsolve`: Nonlinear solver used for the implicit stage (default: `NLNewton()`).
  - `extrapolant`: Initial guess method for the nonlinear solve (default: `:constant`).
  - `new_jac_conv_bound`: Convergence bound below which the Jacobian is reused
    instead of recomputed (default: `1e-3`).
  - `autodiff`: ADTypes.jl backend used for the Jacobian (default: `AutoForwardDiff()`).
  - `concrete_jac`: Whether to build a concrete Jacobian matrix even when a
    matrix-free linear solver is used (default: `nothing`, i.e. decided automatically).

## References

  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential
    Equations", Springer (1992).
"""
struct ImplicitEulerHeun{AD, F, N, T2, T3, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T3
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
end
function ImplicitEulerHeun(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitEulerHeun(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac),
    )
end

"""
    ImplicitRKMil(;
        theta = 1, symplectic = false, linsolve = nothing,
        nlsolve = NLNewton(), extrapolant = :constant,
        new_jac_conv_bound = 1.0e-3,
        interpretation = AlgorithmInterpretation.Ito,
        autodiff = AutoForwardDiff(), concrete_jac = nothing
    )

**ImplicitRKMil: Drift-Implicit Runge-Kutta Milstein (Stiff)**

Drift-implicit theta-method version of the derivative-free Runge-Kutta Milstein
scheme. The Milstein correction gives strong order 1.0, twice the strong order of
[`ImplicitEM`](@ref), at the cost of being restricted to diagonal and scalar noise.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive
  - **Noise types**: Diagonal and scalar noise only
  - **SDE interpretation**: Itô or Stratonovich (selected by `interpretation`)

## When to Use

  - Stiff SDEs with diagonal or scalar noise where strong order 1.0 is needed
  - When [`ImplicitEM`](@ref) converges too slowly in the strong sense
  - Pathwise-accurate simulations (e.g. multilevel Monte Carlo) of stiff problems

Note that `ImplicitRKMil` is only compatible with diagonal-noise problems; for
non-diagonal noise use [`ImplicitEM`](@ref) or [`ISSEM`](@ref).

## Keyword Arguments

  - `theta`: Implicitness parameter of the drift discretization (default: `1`, i.e.
    drift-implicit Euler). `theta = 1/2` gives the trapezoidal rule.
  - `symplectic`: If `true`, uses the symplectic midpoint discretization
    (forces `theta = 1/2`).
  - `interpretation`: `AlgorithmInterpretation.Ito` (default) or
    `AlgorithmInterpretation.Stratonovich`, selecting which SDE interpretation the
    Milstein correction is built for.
  - `linsolve`: LinearSolve.jl algorithm used for the Newton linear systems
    (default: `nothing`, i.e. the LinearSolve.jl default).
  - `nlsolve`: Nonlinear solver used for the implicit stage (default: `NLNewton()`).
  - `extrapolant`: Initial guess method for the nonlinear solve (default: `:constant`).
  - `new_jac_conv_bound`: Convergence bound below which the Jacobian is reused
    instead of recomputed (default: `1e-3`).
  - `autodiff`: ADTypes.jl backend used for the Jacobian (default: `AutoForwardDiff()`).
  - `concrete_jac`: Whether to build a concrete Jacobian matrix even when a
    matrix-free linear solver is used (default: `nothing`, i.e. decided automatically).

## References

  - Milstein, G.N., "Numerical Integration of Stochastic Differential Equations",
    Kluwer (1995).
  - Kloeden, P.E., Platen, E., "Numerical Solution of Stochastic Differential
    Equations", Springer (1992).
"""
struct ImplicitRKMil{AD, F, N, T2, T3, interpretation, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T3
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
end
function ImplicitRKMil(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3, interpretation = SciMLBase.AlgorithmInterpretation.Ito
    )
    autodiff = _fixup_ad(autodiff)
    return ImplicitRKMil{
        typeof(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(symplectic ? 1 / 2 : theta), typeof(new_jac_conv_bound),
        interpretation, typeof(_unwrap_val(concrete_jac)),
    }(
        linsolve, nlsolve, symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac),
    )
end

"""
    ISSEM(;
        theta = 1, symplectic = false, linsolve = nothing,
        nlsolve = NLNewton(), extrapolant = :constant,
        new_jac_conv_bound = 1.0e-3,
        autodiff = AutoForwardDiff(), concrete_jac = nothing
    )

**ISSEM: Implicit Split-Step Euler-Maruyama (Fully Stiff)**

Implicit split-step Euler-Maruyama method for Itô SDEs. Unlike [`ImplicitEM`](@ref),
which is implicit only in the drift, the split-step formulation applies the implicit
solve to an intermediate state and then adds the diffusion, which gives stability for
problems where the *diffusion* term is stiff as well.

## Method Properties

  - **Strong Order**: 0.5
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive
  - **Noise types**: All (diagonal, non-diagonal, scalar)
  - **SDE interpretation**: Itô

## Mathematical Formulation

The step is split into a purely deterministic implicit solve followed by the
diffusion update, so the noise never enters the nonlinear system:

```math
u^* = u_n + dt (θ f(u^*, t_{n+1}) + (1 - θ) f(u_n, t_n))
```

```math
u_{n+1} = u^* + g(u_n, t_n) ΔW_n
```

## When to Use

  - SDEs that are stiff in both the drift and the diffusion
  - Problems where [`ImplicitEM`](@ref) still requires very small time steps
  - Mean-square stability sensitive problems (e.g. large multiplicative noise)

For Stratonovich problems, use [`ISSEulerHeun`](@ref).

## Keyword Arguments

  - `theta`: Implicitness parameter of the drift solve (default: `1`).
  - `symplectic`: If `true`, uses the symplectic midpoint discretization
    (forces `theta = 1/2`).
  - `linsolve`: LinearSolve.jl algorithm used for the Newton linear systems
    (default: `nothing`, i.e. the LinearSolve.jl default).
  - `nlsolve`: Nonlinear solver used for the implicit stage (default: `NLNewton()`).
  - `extrapolant`: Initial guess method for the nonlinear solve (default: `:constant`).
  - `new_jac_conv_bound`: Convergence bound below which the Jacobian is reused
    instead of recomputed (default: `1e-3`).
  - `autodiff`: ADTypes.jl backend used for the Jacobian (default: `AutoForwardDiff()`).
  - `concrete_jac`: Whether to build a concrete Jacobian matrix even when a
    matrix-free linear solver is used (default: `nothing`, i.e. decided automatically).

## References

  - Higham, D.J., Mao, X., Stuart, A.M., "Strong Convergence of Euler-Type Methods
    for Nonlinear Stochastic Differential Equations", SIAM J. Numer. Anal., 40 (3),
    pp. 1041–1063 (2002). DOI: 10.1137/S0036142901389530.
"""
struct ISSEM{AD, F, N, T2, T3, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T3
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
end
function ISSEM(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
    )
    autodiff = _fixup_ad(autodiff)
    return ISSEM(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac),
    )
end

"""
    ISSEulerHeun(;
        theta = 1, symplectic = false, linsolve = nothing,
        nlsolve = NLNewton(), extrapolant = :constant,
        new_jac_conv_bound = 1.0e-3,
        autodiff = AutoForwardDiff(), concrete_jac = nothing
    )

**ISSEulerHeun: Implicit Split-Step Euler-Heun (Fully Stiff, Stratonovich)**

Implicit split-step Euler-Heun method, the Stratonovich counterpart of
[`ISSEM`](@ref). The drift is advanced with an implicit solve that contains no noise
terms, and the Stratonovich diffusion is then applied with the Heun
predictor-corrector, giving stability for problems that are stiff in both the drift
and the diffusion.

## Method Properties

  - **Strong Order**: 0.5
  - **Weak Order**: 1.0
  - **Time stepping**: Adaptive
  - **Noise types**: All (diagonal, non-diagonal, scalar)
  - **SDE interpretation**: Stratonovich

## When to Use

  - Stratonovich SDEs that are stiff in both drift and diffusion
  - Problems where [`ImplicitEulerHeun`](@ref) still requires very small time steps
  - Physical models with large multiplicative noise in the Stratonovich sense

## Keyword Arguments

  - `theta`: Implicitness parameter of the drift solve (default: `1`).
  - `symplectic`: If `true`, uses the symplectic midpoint discretization
    (forces `theta = 1/2`).
  - `linsolve`: LinearSolve.jl algorithm used for the Newton linear systems
    (default: `nothing`, i.e. the LinearSolve.jl default).
  - `nlsolve`: Nonlinear solver used for the implicit stage (default: `NLNewton()`).
  - `extrapolant`: Initial guess method for the nonlinear solve (default: `:constant`).
  - `new_jac_conv_bound`: Convergence bound below which the Jacobian is reused
    instead of recomputed (default: `1e-3`).
  - `autodiff`: ADTypes.jl backend used for the Jacobian (default: `AutoForwardDiff()`).
  - `concrete_jac`: Whether to build a concrete Jacobian matrix even when a
    matrix-free linear solver is used (default: `nothing`, i.e. decided automatically).

## References

  - Higham, D.J., Mao, X., Stuart, A.M., "Strong Convergence of Euler-Type Methods
    for Nonlinear Stochastic Differential Equations", SIAM J. Numer. Anal., 40 (3),
    pp. 1041–1063 (2002). DOI: 10.1137/S0036142901389530.
"""
struct ISSEulerHeun{AD, F, N, T2, T3, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    theta::T2
    extrapolant::Symbol
    new_jac_conv_bound::T3
    symplectic::Bool
    autodiff::AD
    concrete_jac::CJ
end
function ISSEulerHeun(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        extrapolant = :constant,
        theta = 1, symplectic = false,
        new_jac_conv_bound = 1.0e-3,
    )
    autodiff = _fixup_ad(autodiff)
    return ISSEulerHeun(
        linsolve, nlsolve,
        symplectic ? 1 / 2 : theta,
        extrapolant,
        new_jac_conv_bound, symplectic,
        autodiff, _unwrap_val(concrete_jac),
    )
end

"""
    SKenCarp(;
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :min_correct,
        new_jac_conv_bound = 1.0e-3, ode_error_est = true,
        autodiff = AutoForwardDiff(), concrete_jac = nothing
    )

**SKenCarp: Stochastic KenCarp for Additive Noise (Stiff, Recommended)**

Stochastic extension of the KenCarp4 additive Runge-Kutta (ESDIRK) method for SDEs
with additive noise. The deterministic part is integrated with the L-stable,
stiffly-accurate KenCarp4 tableau, while the additive noise is handled with the
higher-order stochastic quadrature of Rößler's SRA schemes, so the method retains
strong order 1.5 on additive-noise problems.

## Method Properties

  - **Strong Order**: 1.5 (for additive noise)
  - **Deterministic Order**: 4 (KenCarp4)
  - **Time stepping**: Adaptive
  - **Noise types**: Additive noise (diagonal, non-diagonal, and scalar)
  - **SDE interpretation**: Itô (Itô and Stratonovich coincide for additive noise)
  - **Stability**: L-stable, stiffly accurate

## When to Use

  - The recommended default for stiff SDEs with additive noise
  - Problems where the deterministic dynamics dominate and need high accuracy
  - Chemical, mechanical, or PDE-discretization models driven by additive forcing

Additive noise means the diffusion function does not depend on the state (`g(u,p,t)`
is independent of `u`); for state-dependent noise use [`ImplicitEM`](@ref),
[`ImplicitRKMil`](@ref), or [`ISSEM`](@ref).

## Keyword Arguments

  - `smooth_est`: Whether to use the smoothed error estimate, which is more robust
    for stiff problems (default: `true`).
  - `ode_error_est`: Whether the deterministic (ODE) embedded error estimate is
    included in the adaptivity controller (default: `true`).
  - `linsolve`: LinearSolve.jl algorithm used for the Newton linear systems
    (default: `nothing`, i.e. the LinearSolve.jl default).
  - `nlsolve`: Nonlinear solver used for the implicit stages (default: `NLNewton()`).
  - `extrapolant`: Initial guess method for the nonlinear solves
    (default: `:min_correct`).
  - `new_jac_conv_bound`: Convergence bound below which the Jacobian is reused
    instead of recomputed (default: `1e-3`).
  - `autodiff`: ADTypes.jl backend used for the Jacobian (default: `AutoForwardDiff()`).
  - `concrete_jac`: Whether to build a concrete Jacobian matrix even when a
    matrix-free linear solver is used (default: `nothing`, i.e. decided automatically).

## References

  - Kennedy, C.A., Carpenter, M.H., "Additive Runge-Kutta schemes for
    convection-diffusion-reaction equations", Applied Numerical Mathematics, 44 (1-2),
    pp. 139–181 (2003). DOI: 10.1016/S0168-9274(02)00138-1.
  - Rackauckas, C., Nie, Q., "Stability-Optimized High Order Methods and Stiffness
    Detection for Pathwise Stiff Stochastic Differential Equations", 2020 IEEE High
    Performance Extreme Computing Conference (HPEC). DOI: 10.1109/HPEC43674.2020.9286178.
  - Rößler A., "Runge-Kutta Methods for the Strong Approximation of Solutions of
    Stochastic Differential Equations", SIAM J. Numer. Anal., 48 (3), pp. 922–952
    (2010). DOI: 10.1137/09076636X.
"""
struct SKenCarp{AD, F, N, T2, CJ} <:
    StochasticDiffEqNewtonAdaptiveAlgorithm
    linsolve::F
    nlsolve::N
    smooth_est::Bool
    extrapolant::Symbol
    new_jac_conv_bound::T2
    ode_error_est::Bool
    autodiff::AD
    concrete_jac::CJ
end

function SKenCarp(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :min_correct,
        new_jac_conv_bound = 1.0e-3,
        ode_error_est = true
    )
    autodiff = _fixup_ad(autodiff)
    return SKenCarp(
        linsolve, nlsolve, smooth_est, extrapolant, new_jac_conv_bound,
        ode_error_est,
        autodiff, _unwrap_val(concrete_jac),
    )
end
