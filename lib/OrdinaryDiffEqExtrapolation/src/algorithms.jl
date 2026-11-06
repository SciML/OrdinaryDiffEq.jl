abstract type OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm <:
OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqImplicitExtrapolationAlgorithm <:
OrdinaryDiffEqAdaptiveImplicitAlgorithm end
"""
    AitkenNeville(;
            max_order = 10, min_order = 1, init_order = 5,
            threading = false
        ) -> AitkenNeville

Adaptive explicit Euler extrapolation using the Aitken-Neville scheme and a Romberg
subdivision sequence. The method varies both its step size and extrapolation order and
is intended for smooth, non-stiff problems at tight tolerances.

# Keywords

- `max_order::Integer = 10`: largest extrapolation order the controller may select.
- `min_order::Integer = 1`: smallest extrapolation order the controller may select.
- `init_order::Integer = 5`: extrapolation order used for the first step. Choose values
  satisfying `1 <= min_order <= init_order <= max_order`.
- `threading = false`: enable concurrent internal extrapolation work when `true`; use
  `false` for serial execution.

# Returns

- `AitkenNeville`: configured adaptive extrapolation algorithm for `CommonSolve.solve`.

# Examples

```julia
using OrdinaryDiffEqExtrapolation: AitkenNeville
using CommonSolve: solve
using SciMLBase: ODEProblem

prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
sol = solve(prob, AitkenNeville(max_order = 8); abstol = 1.0e-10, reltol = 1.0e-8)
```

# References

C. Elrod et al., "Parallelizing explicit and implicit extrapolation methods for
ordinary differential equations," HPEC 2022, pp. 1-9.
"""
Base.@kwdef struct AitkenNeville{TO} <: OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm
    max_order::Int = 10
    min_order::Int = 1
    init_order::Int = 5
    threading::TO = false
end

"""
    ImplicitEulerExtrapolation(;
            autodiff = AutoForwardDiff(), concrete_jac = nothing,
            linsolve = nothing, max_order = 12, min_order = 3, init_order = 5,
            threading = false, sequence = :harmonic
        ) -> ImplicitEulerExtrapolation

Adaptive implicit Euler extrapolation, similar to SEULEX. It is intended for stiff
problems and varies both its step size and extrapolation order.

# Keywords

- `autodiff = AutoForwardDiff()`: ADTypes backend used to construct Jacobians.
- `concrete_jac = nothing`: Jacobian materialization policy. `nothing` lets the solver
  choose; use `true` or `false` to request or disable a concrete Jacobian.
- `linsolve = nothing`: LinearSolve algorithm for implicit stage systems. `nothing`
  selects the OrdinaryDiffEq default.
- `max_order::Integer = 12`: largest internal extrapolation order.
- `min_order::Integer = 3`: smallest internal extrapolation order; values below `3`
  are raised to `3`.
- `init_order::Integer = 5`: initial internal order. It is raised when necessary so
  that `min_order < init_order < max_order`.
- `threading = false`: enable concurrent extrapolation columns when `true`.
- `sequence::Symbol = :harmonic`: subdivision sequence. Supported values are
  `:harmonic`, `:romberg`, and `:bulirsch`; other values warn and use `:harmonic`.

# Returns

- `ImplicitEulerExtrapolation`: configured implicit extrapolation algorithm for
  `CommonSolve.solve`.

# Examples

```julia
using OrdinaryDiffEqExtrapolation: ImplicitEulerExtrapolation
using CommonSolve: solve
using SciMLBase: ODEProblem

prob = ODEProblem((u, p, t) -> -100u, 1.0, (0.0, 1.0))
sol = solve(prob, ImplicitEulerExtrapolation(); abstol = 1.0e-10, reltol = 1.0e-8)
```

# References

C. Elrod et al., "Parallelizing explicit and implicit extrapolation methods for
ordinary differential equations," HPEC 2022, pp. 1-9.
"""
struct ImplicitEulerExtrapolation{AD, F, TO, CJ} <:
    OrdinaryDiffEqImplicitExtrapolationAlgorithm
    linsolve::F
    max_order::Int
    min_order::Int
    init_order::Int
    threading::TO
    sequence::Symbol # Name of the subdividing sequence
    autodiff::AD
    concrete_jac::CJ
end

function ImplicitEulerExtrapolation(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing,
        max_order = 12, min_order = 3, init_order = 5,
        threading = false, sequence = :harmonic
    )
    autodiff = _fixup_ad(autodiff)

    linsolve = (
            linsolve === nothing &&
            (threading == true || threading isa PolyesterThreads)
        ) ?
        GenericLUFactorization() : linsolve

    min_order = max(3, min_order)
    init_order = max(min_order + 1, init_order)
    max_order = max(init_order + 1, max_order)

    # Warn user if orders have been changed
    if (min_order, init_order, max_order) != (min_order, init_order, max_order)
        @warn "The range of extrapolation orders and/or the initial order given to the
          `ImplicitEulerExtrapolation` algorithm are not valid and have been changed:
          Minimal order: " * lpad(min_order, 2, " ") * " --> " * lpad(min_order, 2, " ") *
            "
Maximal order: " * lpad(max_order, 2, " ") * " --> " * lpad(max_order, 2, " ") *
            "
Initial order: " * lpad(init_order, 2, " ") * " --> " * lpad(init_order, 2, " ")
    end

    # Warn user if sequence has been changed:
    if sequence != :harmonic && sequence != :romberg && sequence != :bulirsch
        @warn "The `sequence` given to the `ImplicitEulerExtrapolation` algorithm
            is not valid: it must match `:harmonic`, `:romberg` or `:bulirsch`.
            Thus it has been changed
          :$(sequence) --> :harmonic"
        sequence = :harmonic
    end
    return ImplicitEulerExtrapolation(
        linsolve, max_order, min_order,
        init_order,
        threading, sequence, autodiff,
        _unwrap_val(concrete_jac)

    )
end

"""
    ExtrapolationMidpointDeuflhard(;
            min_order = 1, init_order = 5,
            max_order = 10, sequence = :harmonic, threading = true,
            sequence_factor = 2
        ) -> ExtrapolationMidpointDeuflhard

Adaptive explicit midpoint extrapolation with Deuflhard's order and step-size
selection. Barycentric extrapolation combines the midpoint approximations. The
independent extrapolation columns can run concurrently.

# Keywords

- `min_order::Integer = 1`: smallest internal extrapolation order; values below `1`
  are raised to `1`.
- `init_order::Integer = 5`: initial internal order. It is raised to at least
  `min_order`.
- `max_order::Integer = 10`: largest internal order. It is raised to at least
  `init_order`.
- `sequence::Symbol = :harmonic`: subdivision sequence. Supported values are
  `:harmonic`, `:romberg`, and `:bulirsch`; other values warn and use `:harmonic`.
- `threading = true`: enable concurrent extrapolation columns when `true`.
- `sequence_factor::Integer = 2`: even multiplier applied to the subdivision
  sequence. Odd values warn and are replaced by `2`.

# Returns

- `ExtrapolationMidpointDeuflhard`: configured explicit extrapolation algorithm for
  `CommonSolve.solve`.

# Examples

```julia
using OrdinaryDiffEqExtrapolation: ExtrapolationMidpointDeuflhard
using CommonSolve: solve
using SciMLBase: ODEProblem

prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
alg = ExtrapolationMidpointDeuflhard(sequence = :romberg, threading = false)
sol = solve(prob, alg; abstol = 1.0e-10, reltol = 1.0e-8)
```

# References

C. Elrod et al., "Parallelizing explicit and implicit extrapolation methods for
ordinary differential equations," HPEC 2022, pp. 1-9.
"""
struct ExtrapolationMidpointDeuflhard{TO} <:
    OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm
    min_order::Int # Minimal extrapolation order
    init_order::Int # Initial extrapolation order
    max_order::Int # Maximal extrapolation order
    sequence::Symbol # Name of the subdividing sequence
    threading::TO
    sequence_factor::Int # An even factor by which sequence is scaled for midpoint extrapolation
end
function ExtrapolationMidpointDeuflhard(;
        min_order = 1, init_order = 5, max_order = 10,
        sequence = :harmonic, threading = true,
        sequence_factor = 2
    )
    # Enforce 1 <=  min_order <= init_order <= max_order:
    min_order = max(1, min_order)
    init_order = max(min_order, init_order)
    max_order = max(init_order, max_order)

    # Warn user if orders have been changed
    if (min_order, init_order, max_order) != (min_order, init_order, max_order)
        @warn "The range of extrapolation orders and/or the initial order given to the
          `ExtrapolationMidpointDeuflhard` algorithm are not valid and have been changed:
          Minimal order: " * lpad(min_order, 2, " ") * " --> " * lpad(min_order, 2, " ") *
            "
Maximal order: " * lpad(max_order, 2, " ") * " --> " * lpad(max_order, 2, " ") *
            "
Initial order: " * lpad(init_order, 2, " ") * " --> " * lpad(init_order, 2, " ")
    end

    # Warn user if sequence_factor is not even
    if sequence_factor % 2 != 0
        @warn "A non-even number cannot be used as sequence factor.
              Thus is has been changed
              $(sequence_factor) --> 2"
        sequence_factor = 2
    end

    # Warn user if sequence has been changed:
    if sequence != :harmonic && sequence != :romberg && sequence != :bulirsch
        @warn "The `sequence` given to the `ExtrapolationMidpointDeuflhard` algorithm
           is not valid: it must match `:harmonic`, `:romberg` or `:bulirsch`.
           Thus it has been changed
          :$(sequence) --> :harmonic"
        sequence = :harmonic
    end

    # Initialize algorithm
    return ExtrapolationMidpointDeuflhard(
        min_order, init_order, max_order, sequence, threading,
        sequence_factor
    )
end

"""
    ImplicitDeuflhardExtrapolation(;
            autodiff = AutoForwardDiff(),
            concrete_jac = nothing, linsolve = nothing, min_order = 1,
            init_order = 5, max_order = 10, sequence = :harmonic,
            threading = false
        ) -> ImplicitDeuflhardExtrapolation

Adaptive implicit midpoint extrapolation with Deuflhard's order and step-size
selection. It targets stiff problems and uses barycentric extrapolation of the
implicit midpoint approximations.

# Keywords

- `autodiff = AutoForwardDiff()`: ADTypes backend used to construct Jacobians.
- `concrete_jac = nothing`: Jacobian materialization policy. `nothing` lets the solver
  choose; use `true` or `false` to request or disable a concrete Jacobian.
- `linsolve = nothing`: LinearSolve algorithm for implicit stage systems. `nothing`
  selects the OrdinaryDiffEq default.
- `min_order::Integer = 1`: smallest internal extrapolation order; values below `1`
  are raised to `1`.
- `init_order::Integer = 5`: initial internal order. It is raised to at least
  `min_order`.
- `max_order::Integer = 10`: largest internal order. It is raised to at least
  `init_order`.
- `sequence::Symbol = :harmonic`: subdivision sequence. Supported values are
  `:harmonic`, `:romberg`, and `:bulirsch`; other values warn and use `:harmonic`.
- `threading = false`: enable concurrent extrapolation columns when `true`.

# Returns

- `ImplicitDeuflhardExtrapolation`: configured implicit extrapolation algorithm for
  `CommonSolve.solve`.

# Examples

```julia
using OrdinaryDiffEqExtrapolation: ImplicitDeuflhardExtrapolation
using CommonSolve: solve
using SciMLBase: ODEProblem

prob = ODEProblem((u, p, t) -> -100u, 1.0, (0.0, 1.0))
alg = ImplicitDeuflhardExtrapolation(sequence = :romberg)
sol = solve(prob, alg; abstol = 1.0e-10, reltol = 1.0e-8)
```

# References

C. Elrod et al., "Parallelizing explicit and implicit extrapolation methods for
ordinary differential equations," HPEC 2022, pp. 1-9.
"""
struct ImplicitDeuflhardExtrapolation{AD, F, TO, CJ} <:
    OrdinaryDiffEqImplicitExtrapolationAlgorithm
    linsolve::F
    min_order::Int # Minimal extrapolation order
    init_order::Int # Initial extrapolation order
    max_order::Int # Maximal extrapolation order
    sequence::Symbol # Name of the subdividing sequence
    threading::TO
    autodiff::AD
    concrete_jac::CJ
end
function ImplicitDeuflhardExtrapolation(;
        autodiff = AutoForwardDiff(),
        concrete_jac = nothing,
        linsolve = nothing,
        min_order = 1, init_order = 5, max_order = 10,
        sequence = :harmonic, threading = false
    )
    autodiff = _fixup_ad(autodiff)

    # Enforce 1 <=  min_order <= init_order <= max_order:
    min_order = max(1, min_order)
    init_order = max(min_order, init_order)
    max_order = max(init_order, max_order)

    linsolve = (
            linsolve === nothing &&
            (threading == true || threading isa PolyesterThreads)
        ) ?
        GenericLUFactorization() : linsolve

    # Warn user if orders have been changed
    if (min_order, init_order, max_order) != (min_order, init_order, max_order)
        @warn "The range of extrapolation orders and/or the initial order given to the
          `ImplicitDeuflhardExtrapolation` algorithm are not valid and have been changed:
          Minimal order: " * lpad(min_order, 2, " ") * " --> " * lpad(min_order, 2, " ") *
            "
Maximal order: " * lpad(max_order, 2, " ") * " --> " * lpad(max_order, 2, " ") *
            "
Initial order: " * lpad(init_order, 2, " ") * " --> " * lpad(init_order, 2, " ")
    end

    # Warn user if sequence has been changed:
    if sequence != :harmonic && sequence != :romberg && sequence != :bulirsch
        @warn "The `sequence` given to the `ImplicitDeuflhardExtrapolation` algorithm
           is not valid: it must match `:harmonic`, `:romberg` or `:bulirsch`.
           Thus it has been changed
          :$(sequence) --> :harmonic"
        sequence = :harmonic
    end

    # Initialize algorithm
    return ImplicitDeuflhardExtrapolation(
        linsolve, min_order,
        init_order, max_order,
        sequence, threading, autodiff,
        _unwrap_val(concrete_jac)

    )
end

"""
    ExtrapolationMidpointHairerWanner(;
            min_order = 2, init_order = 5,
            max_order = 10, sequence = :harmonic, threading = true,
            sequence_factor = 2
        ) -> ExtrapolationMidpointHairerWanner

Adaptive explicit midpoint extrapolation with the order and step-size controller from
Hairer and Wanner's ODEX algorithm. Barycentric extrapolation combines the midpoint
approximations, whose independent columns can run concurrently.

# Keywords

- `min_order::Integer = 2`: smallest internal extrapolation order; values below `2`
  are raised to `2`.
- `init_order::Integer = 5`: initial internal order. It is raised when necessary so
  that `min_order < init_order`.
- `max_order::Integer = 10`: largest internal order. It is raised when necessary so
  that `init_order < max_order`.
- `sequence::Symbol = :harmonic`: subdivision sequence. Supported values are
  `:harmonic`, `:romberg`, and `:bulirsch`; other values warn and use `:harmonic`.
- `threading = true`: enable concurrent extrapolation columns when `true`.
- `sequence_factor::Integer = 2`: even multiplier applied to the subdivision
  sequence. Odd values warn and are replaced by `2`.

# Returns

- `ExtrapolationMidpointHairerWanner`: configured explicit extrapolation algorithm
  for `CommonSolve.solve`.

# Examples

```julia
using OrdinaryDiffEqExtrapolation: ExtrapolationMidpointHairerWanner
using CommonSolve: solve
using SciMLBase: ODEProblem

prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
alg = ExtrapolationMidpointHairerWanner(threading = false)
sol = solve(prob, alg; abstol = 1.0e-10, reltol = 1.0e-8)
```

# References

C. Elrod et al., "Parallelizing explicit and implicit extrapolation methods for
ordinary differential equations," HPEC 2022, pp. 1-9.
"""
struct ExtrapolationMidpointHairerWanner{TO} <:
    OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm
    min_order::Int # Minimal extrapolation order
    init_order::Int # Initial extrapolation order
    max_order::Int # Maximal extrapolation order
    sequence::Symbol # Name of the subdividing sequence
    threading::TO
    sequence_factor::Int # An even factor by which sequence is scaled for midpoint extrapolation
end
function ExtrapolationMidpointHairerWanner(;
        min_order = 2, init_order = 5, max_order = 10,
        sequence = :harmonic, threading = true,
        sequence_factor = 2
    )
    # Enforce 2 <=  min_order
    # and min_order + 1 <= init_order <= max_order - 1:
    min_order = max(2, min_order)
    init_order = max(min_order + 1, init_order)
    max_order = max(init_order + 1, max_order)

    # Warn user if orders have been changed
    if (min_order, init_order, max_order) != (min_order, init_order, max_order)
        @warn "The range of extrapolation orders and/or the initial order given to the
          `ExtrapolationMidpointHairerWanner` algorithm are not valid and have been changed:
          Minimal order: " * lpad(min_order, 2, " ") * " --> " * lpad(min_order, 2, " ") *
            "
Maximal order: " * lpad(max_order, 2, " ") * " --> " * lpad(max_order, 2, " ") *
            "
Initial order: " * lpad(init_order, 2, " ") * " --> " * lpad(init_order, 2, " ")
    end

    # Warn user if sequence_factor is not even
    if sequence_factor % 2 != 0
        @warn "A non-even number cannot be used as sequence factor.
              Thus is has been changed
              $(sequence_factor) --> 2"
        sequence_factor = 2
    end

    # Warn user if sequence has been changed:
    if sequence != :harmonic && sequence != :romberg && sequence != :bulirsch
        @warn "The `sequence` given to the `ExtrapolationMidpointHairerWanner` algorithm
           is not valid: it must match `:harmonic`, `:romberg` or `:bulirsch`.
           Thus it has been changed
          :$(sequence) --> :harmonic"
        sequence = :harmonic
    end

    # Initialize algorithm
    return ExtrapolationMidpointHairerWanner(
        min_order, init_order, max_order, sequence, threading,
        sequence_factor
    )
end

"""
    ImplicitHairerWannerExtrapolation(;
            autodiff = AutoForwardDiff(),
            concrete_jac = nothing, linsolve = nothing, min_order = 2,
            init_order = 5, max_order = 10, sequence = :harmonic,
            threading = false
        ) -> ImplicitHairerWannerExtrapolation

Adaptive implicit midpoint extrapolation with the order and step-size controller from
Hairer and Wanner's SODEX algorithm. It is intended for stiff problems and supports
parallel evaluation of independent extrapolation columns.

# Keywords

- `autodiff = AutoForwardDiff()`: ADTypes backend used to construct Jacobians.
- `concrete_jac = nothing`: Jacobian materialization policy. `nothing` lets the solver
  choose; use `true` or `false` to request or disable a concrete Jacobian.
- `linsolve = nothing`: LinearSolve algorithm for implicit stage systems. `nothing`
  selects the OrdinaryDiffEq default.
- `min_order::Integer = 2`: smallest internal extrapolation order; values below `2`
  are raised to `2`.
- `init_order::Integer = 5`: initial internal order. It is raised when necessary so
  that `min_order < init_order`.
- `max_order::Integer = 10`: largest internal order. It is raised when necessary so
  that `init_order < max_order`.
- `sequence::Symbol = :harmonic`: subdivision sequence. Supported values are
  `:harmonic`, `:romberg`, and `:bulirsch`; other values warn and use `:harmonic`.
- `threading = false`: enable concurrent extrapolation columns when `true`.

# Returns

- `ImplicitHairerWannerExtrapolation`: configured implicit extrapolation algorithm
  for `CommonSolve.solve`.

# Examples

```julia
using OrdinaryDiffEqExtrapolation: ImplicitHairerWannerExtrapolation
using CommonSolve: solve
using SciMLBase: ODEProblem

prob = ODEProblem((u, p, t) -> -100u, 1.0, (0.0, 1.0))
alg = ImplicitHairerWannerExtrapolation(sequence = :bulirsch)
sol = solve(prob, alg; abstol = 1.0e-10, reltol = 1.0e-8)
```

# References

C. Elrod et al., "Parallelizing explicit and implicit extrapolation methods for
ordinary differential equations," HPEC 2022, pp. 1-9.
"""
struct ImplicitHairerWannerExtrapolation{AD, F, TO, CJ} <:
    OrdinaryDiffEqImplicitExtrapolationAlgorithm
    linsolve::F
    min_order::Int # Minimal extrapolation order
    init_order::Int # Initial extrapolation order
    max_order::Int # Maximal extrapolation order
    sequence::Symbol # Name of the subdividing sequence
    threading::TO
    autodiff::AD
    concrete_jac::CJ
end

function ImplicitHairerWannerExtrapolation(;
        autodiff = AutoForwardDiff(),

        concrete_jac = nothing,
        linsolve = nothing,
        min_order = 2, init_order = 5, max_order = 10,
        sequence = :harmonic, threading = false
    )

    # Enforce 2 <=  min_order
    # and min_order + 1 <= init_order <= max_order - 1:
    min_order = max(2, min_order)
    init_order = max(min_order + 1, init_order)
    max_order = max(init_order + 1, max_order)

    linsolve = (
            linsolve === nothing &&
            (threading == true || threading isa PolyesterThreads)
        ) ?
        GenericLUFactorization() : linsolve

    # Warn user if orders have been changed
    if (min_order, init_order, max_order) != (min_order, init_order, max_order)
        @warn "The range of extrapolation orders and/or the initial order given to the
          `ImplicitHairerWannerExtrapolation` algorithm are not valid and have been changed:
          Minimal order: " * lpad(min_order, 2, " ") * " --> " * lpad(min_order, 2, " ") *
            "
Maximal order: " * lpad(max_order, 2, " ") * " --> " * lpad(max_order, 2, " ") *
            "
Initial order: " * lpad(init_order, 2, " ") * " --> " * lpad(init_order, 2, " ")
    end

    # Warn user if sequence has been changed:
    if sequence != :harmonic && sequence != :romberg && sequence != :bulirsch
        @warn "The `sequence` given to the `ImplicitHairerWannerExtrapolation` algorithm
           is not valid: it must match `:harmonic`, `:romberg` or `:bulirsch`.
           Thus it has been changed
          :$(sequence) --> :harmonic"
        sequence = :harmonic
    end

    autodiff = _fixup_ad(autodiff)
    # Initialize algorithm
    return ImplicitHairerWannerExtrapolation(
        linsolve, min_order,
        init_order,
        max_order, sequence, threading, autodiff,
        _unwrap_val(concrete_jac)

    )
end

"""
    ImplicitEulerBarycentricExtrapolation(;
            autodiff = AutoForwardDiff(),
            concrete_jac = nothing, linsolve = nothing, min_order = 3,
            init_order = 5, max_order = 12, sequence = :harmonic,
            threading = false, sequence_factor = 2
        ) -> ImplicitEulerBarycentricExtrapolation

Adaptive implicit Euler extrapolation using barycentric coordinates and the
Hairer-Wanner order and step-size strategy. It targets stiff problems and can evaluate
independent extrapolation columns concurrently.

# Keywords

- `autodiff = AutoForwardDiff()`: ADTypes backend used to construct Jacobians.
- `concrete_jac = nothing`: Jacobian materialization policy. `nothing` lets the solver
  choose; use `true` or `false` to request or disable a concrete Jacobian.
- `linsolve = nothing`: LinearSolve algorithm for implicit stage systems. `nothing`
  selects the OrdinaryDiffEq default.
- `min_order::Integer = 3`: smallest internal extrapolation order; values below `3`
  are raised to `3`.
- `init_order::Integer = 5`: initial internal order. It is raised when necessary so
  that `min_order < init_order`.
- `max_order::Integer = 12`: largest internal order. It is raised when necessary so
  that `init_order < max_order`.
- `sequence::Symbol = :harmonic`: subdivision sequence. Supported values are
  `:harmonic`, `:romberg`, and `:bulirsch`; other values warn and use `:harmonic`.
- `threading = false`: enable concurrent extrapolation columns when `true`.
- `sequence_factor::Integer = 2`: positive multiplier applied to the subdivision
  sequence used by the implicit Euler base method.

# Returns

- `ImplicitEulerBarycentricExtrapolation`: configured implicit extrapolation algorithm
  for `CommonSolve.solve`.

# Examples

```julia
using OrdinaryDiffEqExtrapolation: ImplicitEulerBarycentricExtrapolation
using CommonSolve: solve
using SciMLBase: ODEProblem

prob = ODEProblem((u, p, t) -> -100u, 1.0, (0.0, 1.0))
alg = ImplicitEulerBarycentricExtrapolation(sequence = :romberg)
sol = solve(prob, alg; abstol = 1.0e-10, reltol = 1.0e-8)
```

# References

C. Elrod et al., "Parallelizing explicit and implicit extrapolation methods for
ordinary differential equations," HPEC 2022, pp. 1-9.
"""
struct ImplicitEulerBarycentricExtrapolation{AD, F, TO, CJ} <:
    OrdinaryDiffEqImplicitExtrapolationAlgorithm
    linsolve::F
    min_order::Int # Minimal extrapolation order
    init_order::Int # Initial extrapolation order
    max_order::Int # Maximal extrapolation order
    sequence::Symbol # Name of the subdividing sequence
    threading::TO
    sequence_factor::Int
    autodiff::AD
    concrete_jac::CJ
end

function ImplicitEulerBarycentricExtrapolation(;
        autodiff = AutoForwardDiff(),

        concrete_jac = nothing,
        linsolve = nothing,
        min_order = 3, init_order = 5,
        max_order = 12, sequence = :harmonic,
        threading = false, sequence_factor = 2
    )
    # Enforce 2 <=  min_order
    # and min_order + 1 <= init_order <= max_order - 1:
    min_order = max(3, min_order)
    init_order = max(min_order + 1, init_order)
    max_order = max(init_order + 1, max_order)

    linsolve = (
            linsolve === nothing &&
            (threading == true || threading isa PolyesterThreads)
        ) ?
        GenericLUFactorization() : linsolve

    # Warn user if orders have been changed
    if (min_order, init_order, max_order) != (min_order, init_order, max_order)
        @warn "The range of extrapolation orders and/or the initial order given to the
          `ImplicitEulerBarycentricExtrapolation` algorithm are not valid and have been changed:
          Minimal order: " * lpad(min_order, 2, " ") * " --> " * lpad(min_order, 2, " ") *
            "
Maximal order: " * lpad(max_order, 2, " ") * " --> " * lpad(max_order, 2, " ") *
            "
Initial order: " * lpad(init_order, 2, " ") * " --> " * lpad(init_order, 2, " ")
    end

    # Warn user if sequence has been changed:
    if sequence != :harmonic && sequence != :romberg && sequence != :bulirsch
        @warn "The `sequence` given to the `ImplicitEulerBarycentricExtrapolation` algorithm
           is not valid: it must match `:harmonic`, `:romberg` or `:bulirsch`.
           Thus it has been changed
          :$(sequence) --> :harmonic"
        sequence = :harmonic
    end

    autodiff = _fixup_ad(autodiff)
    # Initialize algorithm
    return ImplicitEulerBarycentricExtrapolation(
        linsolve,
        min_order,
        init_order,
        max_order,
        sequence,
        threading,
        sequence_factor,
        autodiff,
        _unwrap_val(concrete_jac)

    )
end
