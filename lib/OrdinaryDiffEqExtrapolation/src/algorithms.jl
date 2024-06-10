abstract type OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm <:
              OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS, AD, FDT, ST, CJ} <:
              OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD, FDT, ST, CJ} end

"""
AitkenNeville: Parallelized Explicit Extrapolation Method
Euler extrapolation using Aitken-Neville with the Romberg Sequence.
"""
struct AitkenNeville{TO} <: OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm
    max_order::Int
    min_order::Int
    init_order::Int
    threading::TO
end
function AitkenNeville(; max_order = 10, min_order = 1, init_order = 5, threading = false)
    AitkenNeville(max_order, min_order, init_order, threading)
end
"""
ImplicitEulerExtrapolation: Parallelized Implicit Extrapolation Method
Extrapolation of implicit Euler method with Romberg sequence.
Similar to Hairer's SEULEX.
"""
struct ImplicitEulerExtrapolation{CS, AD, F, P, FDT, ST, CJ, TO} <:
       OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    max_order::Int
    min_order::Int
    init_order::Int
    threading::TO
    sequence::Symbol # Name of the subdividing sequence
end

function ImplicitEulerExtrapolation(; chunk_size = Val{0}(), autodiff = true,
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}, linsolve = nothing,
        precs = DEFAULT_PRECS,
        max_order = 12, min_order = 3, init_order = 5,
        threading = false, sequence = :harmonic)
    linsolve = (linsolve === nothing &&
                (threading == true || threading isa PolyesterThreads)) ?
               RFLUFactorization(; thread = Val(false)) : linsolve

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
    ImplicitEulerExtrapolation{_unwrap_val(chunk_size), _unwrap_val(autodiff),
        typeof(linsolve), typeof(precs), diff_type,
        _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(threading)}(linsolve, precs, max_order, min_order,
        init_order,
        threading, sequence)
end
"""
ExtrapolationMidpointDeuflhard: Parallelized Explicit Extrapolation Method
Midpoint extrapolation using Barycentric coordinates
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
function ExtrapolationMidpointDeuflhard(; min_order = 1, init_order = 5, max_order = 10,
        sequence = :harmonic, threading = true,
        sequence_factor = 2)
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
    ExtrapolationMidpointDeuflhard(min_order, init_order, max_order, sequence, threading,
        sequence_factor)
end
"""
ImplicitDeuflhardExtrapolation: Parallelized Implicit Extrapolation Method
Midpoint extrapolation using Barycentric coordinates
"""
struct ImplicitDeuflhardExtrapolation{CS, AD, F, P, FDT, ST, CJ, TO} <:
       OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    min_order::Int # Minimal extrapolation order
    init_order::Int # Initial extrapolation order
    max_order::Int # Maximal extrapolation order
    sequence::Symbol # Name of the subdividing sequence
    threading::TO
end
function ImplicitDeuflhardExtrapolation(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        linsolve = nothing, precs = DEFAULT_PRECS,
        diff_type = Val{:forward},
        min_order = 1, init_order = 5, max_order = 10,
        sequence = :harmonic, threading = false)
    # Enforce 1 <=  min_order <= init_order <= max_order:
    min_order = max(1, min_order)
    init_order = max(min_order, init_order)
    max_order = max(init_order, max_order)

    linsolve = (linsolve === nothing &&
                (threading == true || threading isa PolyesterThreads)) ?
               RFLUFactorization(; thread = Val(false)) : linsolve

    # Warn user if orders have been changed
    if (min_order, init_order, max_order) != (min_order, init_order, max_order)
        @warn "The range of extrapolation orders and/or the initial order given to the
          `ImplicitDeuflhardExtrapolation` algorithm are not valid and have been changed:
          Minimal order: " * lpad(min_order, 2, " ") * " --> " * lpad(min_order, 2, " ") *
              "
Maximal order: " * lpad(max_order, 2, " ") * " --> " * lpad(max_order, 2, " ") *
              "
Initial order: " * lpad(init_order, 2, " ") * " --> " * lpad(init_order, 2, " ")
        chunk_size
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
    ImplicitDeuflhardExtrapolation{_unwrap_val(chunk_size), _unwrap_val(autodiff),
        typeof(linsolve), typeof(precs), diff_type,
        _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(threading)}(linsolve, precs, min_order,
        init_order, max_order,
        sequence, threading)
end
"""
ExtrapolationMidpointHairerWanner: Parallelized Explicit Extrapolation Method
Midpoint extrapolation using Barycentric coordinates, following Hairer's ODEX in the adaptivity behavior.
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
function ExtrapolationMidpointHairerWanner(; min_order = 2, init_order = 5, max_order = 10,
        sequence = :harmonic, threading = true,
        sequence_factor = 2)
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
    ExtrapolationMidpointHairerWanner(
        min_order, init_order, max_order, sequence, threading,
        sequence_factor)
end
"""
ImplicitHairerWannerExtrapolation: Parallelized Implicit Extrapolation Method
Midpoint extrapolation using Barycentric coordinates, following Hairer's SODEX in the adaptivity behavior.
"""
struct ImplicitHairerWannerExtrapolation{CS, AD, F, P, FDT, ST, CJ, TO} <:
       OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    min_order::Int # Minimal extrapolation order
    init_order::Int # Initial extrapolation order
    max_order::Int # Maximal extrapolation order
    sequence::Symbol # Name of the subdividing sequence
    threading::TO
end

function ImplicitHairerWannerExtrapolation(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(),
        concrete_jac = nothing,
        linsolve = nothing, precs = DEFAULT_PRECS,
        diff_type = Val{:forward},
        min_order = 2, init_order = 5, max_order = 10,
        sequence = :harmonic, threading = false)
    # Enforce 2 <=  min_order
    # and min_order + 1 <= init_order <= max_order - 1:
    min_order = max(2, min_order)
    init_order = max(min_order + 1, init_order)
    max_order = max(init_order + 1, max_order)

    linsolve = (linsolve === nothing &&
                (threading == true || threading isa PolyesterThreads)) ?
               RFLUFactorization(; thread = Val(false)) : linsolve

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

    # Initialize algorithm
    ImplicitHairerWannerExtrapolation{_unwrap_val(chunk_size), _unwrap_val(autodiff),
        typeof(linsolve), typeof(precs), diff_type,
        _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(threading)}(linsolve, precs, min_order,
        init_order,
        max_order, sequence, threading)
end

"""
ImplicitEulerBarycentricExtrapolation: Parallelized Implicit Extrapolation Method
Euler extrapolation using Barycentric coordinates, following Hairer's SODEX in the adaptivity behavior.
"""
struct ImplicitEulerBarycentricExtrapolation{CS, AD, F, P, FDT, ST, CJ, TO} <:
       OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    min_order::Int # Minimal extrapolation order
    init_order::Int # Initial extrapolation order
    max_order::Int # Maximal extrapolation order
    sequence::Symbol # Name of the subdividing sequence
    threading::TO
    sequence_factor::Int
end

function ImplicitEulerBarycentricExtrapolation(; chunk_size = Val{0}(),
        autodiff = Val{true}(),
        standardtag = Val{true}(),
        concrete_jac = nothing,
        linsolve = nothing, precs = DEFAULT_PRECS,
        diff_type = Val{:forward},
        min_order = 3, init_order = 5,
        max_order = 12, sequence = :harmonic,
        threading = false, sequence_factor = 2)
    # Enforce 2 <=  min_order
    # and min_order + 1 <= init_order <= max_order - 1:
    min_order = max(3, min_order)
    init_order = max(min_order + 1, init_order)
    max_order = max(init_order + 1, max_order)

    linsolve = (linsolve === nothing &&
                (threading == true || threading isa PolyesterThreads)) ?
               RFLUFactorization(; thread = Val(false)) : linsolve

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

    # Initialize algorithm
    ImplicitEulerBarycentricExtrapolation{_unwrap_val(chunk_size), _unwrap_val(autodiff),
        typeof(linsolve), typeof(precs), diff_type,
        _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(threading)}(linsolve,
        precs,
        min_order,
        init_order,
        max_order,
        sequence,
        threading,
        sequence_factor)
end
