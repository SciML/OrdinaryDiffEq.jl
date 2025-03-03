abstract type OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm <:
              OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS, AD, FDT, ST, CJ} <:
              OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD, FDT, ST, CJ} end
reference = """@inproceedings{elrod2022parallelizing,
  title={Parallelizing explicit and implicit extrapolation methods for ordinary differential equations},
  author={Elrod, Chris and Ma, Yingbo and Althaus, Konstantin and Rackauckas, Christopher and others},
  booktitle={2022 IEEE High Performance Extreme Computing Conference (HPEC)},
  pages={1--9},
  year={2022},
  organization={IEEE}}
"""

@doc generic_solver_docstring(
    "Euler extrapolation using Aitken-Neville with the Romberg Sequence.",
    "AitkenNeville",
    "Parallelized Explicit Extrapolation Method.",
    reference,
    """
    - `max_order`: maximum order of the adaptive order algorithm.
    - `min_order`: minimum order of the adaptive order algorithm.
    - `init_order`: initial order of the adaptive order algorithm.
    - `thread`: determines whether internal broadcasting on appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when Julia is started with multiple threads.
    """,
    """
    max_order::Int = 10,
    min_order::Int = 1,
    init_order = 3,
    thread = OrdinaryDiffEq.False(),
    """)
Base.@kwdef struct AitkenNeville{TO} <: OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm
    max_order::Int = 10
    min_order::Int = 1
    init_order::Int = 5
    threading::TO = false
end

@doc differentiation_rk_docstring(
    "Extrapolation of implicit Euler method with Romberg sequence.
Similar to Hairer's SEULEX.",
    "ImplicitEulerExtrapolation",
    "Parallelized Explicit Extrapolation Method.",
    references = reference,
    extra_keyword_description = """
    - `max_order`: maximum order of the adaptive order algorithm.
    - `min_order`: minimum order of the adaptive order algorithm.
    - `init_order`: initial order of the adaptive order algorithm.
    - `thread`: determines whether internal broadcasting on appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when Julia is started with multiple threads.
    - `sequence`: the step-number sequences, also called the subdividing sequence. Possible values are `:harmonic`, `:romberg` or `:bulirsch`.
    """,
    extra_keyword_default = """
    max_order = 12,
    min_order = 3,
    init_order = 5,
    thread = OrdinaryDiffEq.False(),
    sequence = :harmonic
    """)
struct ImplicitEulerExtrapolation{CS, AD, F, P, FDT, ST, CJ, TO} <:
       OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    max_order::Int
    min_order::Int
    init_order::Int
    threading::TO
    sequence::Symbol # Name of the subdividing sequence
    autodiff::AD
end

function ImplicitEulerExtrapolation(; chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward}(), linsolve = nothing,
        precs = DEFAULT_PRECS,
        max_order = 12, min_order = 3, init_order = 5,
        threading = false, sequence = :harmonic)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

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
    ImplicitEulerExtrapolation{_unwrap_val(chunk_size), typeof(AD_choice),
        typeof(linsolve), typeof(precs), diff_type,
        _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(threading)}(linsolve, precs, max_order, min_order,
        init_order,
        threading, sequence, AD_choice)
end

@doc generic_solver_docstring("Midpoint extrapolation using Barycentric coordinates.",
    "ExtrapolationMidpointDeuflhard",
    "Parallelized Explicit Extrapolation Method.",
    reference,
    """
    - `max_order`: maximum order of the adaptive order algorithm.
    - `min_order`: minimum order of the adaptive order algorithm.
    - `init_order`: initial order of the adaptive order algorithm.
    - `thread`: determines whether internal broadcasting on appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when Julia is started with multiple threads.
    - `sequence`: the step-number sequences, also called the subdividing sequence. Possible values are `:harmonic`, `:romberg` or `:bulirsch`.
    - `sequence_factor`: denotes which even multiple of sequence to take while evaluating internal discretizations.
    """,
    """
    max_order = 10,
    min_order = 1,
    init_order = 5,
    thread = OrdinaryDiffEq.True(),
    sequence = :harmonic,
    sequence_factor = 2,
    """)
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

@doc differentiation_rk_docstring("Midpoint extrapolation using Barycentric coordinates.",
    "ImplicitDeuflhardExtrapolation",
    "Parallelized Explicit Extrapolation Method.",
    references = reference,
    extra_keyword_description = """
    - `max_order`: maximum order of the adaptive order algorithm.
    - `min_order`: minimum order of the adaptive order algorithm.
    - `init_order`: initial order of the adaptive order algorithm.
    - `thread`: determines whether internal broadcasting on appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when Julia is started with multiple threads.
    - `sequence`: the step-number sequences, also called the subdividing sequence. Possible values are `:harmonic`, `:romberg` or `:bulirsch`.
    """,
    extra_keyword_default = """
    max_order = 10,
    min_order = 1,
    init_order = 5,
    thread = OrdinaryDiffEq.False(),
    sequence = :harmonic,
    """)
struct ImplicitDeuflhardExtrapolation{CS, AD, F, P, FDT, ST, CJ, TO} <:
       OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    min_order::Int # Minimal extrapolation order
    init_order::Int # Initial extrapolation order
    max_order::Int # Maximal extrapolation order
    sequence::Symbol # Name of the subdividing sequence
    threading::TO
    autodiff::AD
end
function ImplicitDeuflhardExtrapolation(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(), concrete_jac = nothing,
        linsolve = nothing, precs = DEFAULT_PRECS,
        diff_type = Val{:forward}(),
        min_order = 1, init_order = 5, max_order = 10,
        sequence = :harmonic, threading = false)
    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)

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
    ImplicitDeuflhardExtrapolation{_unwrap_val(chunk_size), typeof(AD_choice),
        typeof(linsolve), typeof(precs), diff_type,
        _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(threading)}(linsolve, precs, min_order,
        init_order, max_order,
        sequence, threading, AD_choice)
end

@doc generic_solver_docstring("Midpoint extrapolation using Barycentric coordinates,
    following Hairer's ODEX in the adaptivity behavior.",
    "ExtrapolationMidpointHairerWanner",
    "Parallelized Explicit Extrapolation Method.",
    reference,
    """
    - `max_order`: maximum order of the adaptive order algorithm.
    - `min_order`: minimum order of the adaptive order algorithm.
    - `init_order`: initial order of the adaptive order algorithm.
    - `thread`: determines whether internal broadcasting on appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when Julia is started with multiple threads.
    - `sequence`: the step-number sequences, also called the subdividing sequence. Possible values are `:harmonic`, `:romberg` or `:bulirsch`.
    - `sequence_factor`: denotes which even multiple of sequence to take while evaluating internal discretizations.
    """,
    """
    max_order = 10,
    min_order = 2,
    init_order = 5,
    thread = OrdinaryDiffEq.True(),
    sequence = :harmonic,
    sequence_factor = 2,
    """)
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

@doc differentiation_rk_docstring("Midpoint extrapolation using Barycentric coordinates,
    following Hairer's SODEX in the adaptivity behavior.",
    "ImplicitHairerWannerExtrapolation",
    "Parallelized Explicit Extrapolation Method.",
    references = reference,
    extra_keyword_description = """
    - `max_order`: maximum order of the adaptive order algorithm.
    - `min_order`: minimum order of the adaptive order algorithm.
    - `init_order`: initial order of the adaptive order algorithm.
    - `thread`: determines whether internal broadcasting on appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when Julia is started with multiple threads.
    - `sequence`: the step-number sequences, also called the subdividing sequence. Possible values are `:harmonic`, `:romberg` or `:bulirsch`.
    """,
    extra_keyword_default = """
    max_order = 10,
    min_order = 2,
    init_order = 5,
    thread = OrdinaryDiffEq.False(),
    sequence = :harmonic,
    """)
struct ImplicitHairerWannerExtrapolation{CS, AD, F, P, FDT, ST, CJ, TO} <:
       OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    min_order::Int # Minimal extrapolation order
    init_order::Int # Initial extrapolation order
    max_order::Int # Maximal extrapolation order
    sequence::Symbol # Name of the subdividing sequence
    threading::TO
    autodiff::AD
end

function ImplicitHairerWannerExtrapolation(;
        chunk_size = Val{0}(), autodiff = AutoForwardDiff(),
        standardtag = Val{true}(),
        concrete_jac = nothing,
        linsolve = nothing, precs = DEFAULT_PRECS,
        diff_type = Val{:forward}(),
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

    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)
    # Initialize algorithm
    ImplicitHairerWannerExtrapolation{_unwrap_val(chunk_size), typeof(AD_choice),
        typeof(linsolve), typeof(precs), diff_type,
        _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(threading)}(linsolve, precs, min_order,
        init_order,
        max_order, sequence, threading, AD_choice)
end

@doc differentiation_rk_docstring("Euler extrapolation using Barycentric coordinates,
    following Hairer's SODEX in the adaptivity behavior.",
    "ImplicitEulerBarycentricExtrapolation",
    "Parallelized Explicit Extrapolation Method.",
    references = reference,
    extra_keyword_description = """
    - `max_order`: maximum order of the adaptive order algorithm.
    - `min_order`: minimum order of the adaptive order algorithm.
    - `init_order`: initial order of the adaptive order algorithm.
    - `thread`: determines whether internal broadcasting on appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when Julia is started with multiple threads.
    - `sequence`: the step-number sequences, also called the subdividing sequence. Possible values are `:harmonic`, `:romberg` or `:bulirsch`.
    - `sequence_factor`: denotes which even multiple of sequence to take while evaluating internal discretizations.
    """,
    extra_keyword_default = """
    max_order = 10,
    min_order = 3,
    init_order = 5,
    thread = OrdinaryDiffEq.False(),
    sequence = :harmonic,
    sequence_factor = 2,
    """)
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
    autodiff::AD
end

function ImplicitEulerBarycentricExtrapolation(; chunk_size = Val{0}(),
        autodiff = AutoForwardDiff(),
        standardtag = Val{true}(),
        concrete_jac = nothing,
        linsolve = nothing, precs = DEFAULT_PRECS,
        diff_type = Val{:forward}(),
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

    AD_choice, chunk_size, diff_type = _process_AD_choice(autodiff, chunk_size, diff_type)
    # Initialize algorithm
    ImplicitEulerBarycentricExtrapolation{_unwrap_val(chunk_size), typeof(AD_choice),
        typeof(linsolve), typeof(precs), diff_type,
        _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(threading)}(linsolve,
        precs,
        min_order,
        init_order,
        max_order,
        sequence,
        threading,
        sequence_factor,
        AD_choice)
end
