abstract type OrdinaryDiffEqAlgorithm <: DiffEqBase.AbstractODEAlgorithm end
abstract type OrdinaryDiffEqAdaptiveAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqCompositeAlgorithm <: OrdinaryDiffEqAlgorithm end

abstract type OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD, FDT, ST, CJ} <:
              OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ} <:
              OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD, FDT, ST, CJ} end
abstract type OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS, AD, FDT, ST, CJ} <:
              OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD, FDT, ST, CJ} end

abstract type OrdinaryDiffEqImplicitAlgorithm{CS, AD, FDT, ST, CJ} <:
              OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ} <:
              OrdinaryDiffEqImplicitAlgorithm{CS, AD, FDT, ST, CJ} end
abstract type OrdinaryDiffEqRosenbrockAlgorithm{CS, AD, FDT, ST, CJ} <:
              OrdinaryDiffEqImplicitAlgorithm{CS, AD, FDT, ST, CJ} end
const NewtonAlgorithm = Union{OrdinaryDiffEqNewtonAlgorithm,
    OrdinaryDiffEqNewtonAdaptiveAlgorithm}
const RosenbrockAlgorithm = Union{OrdinaryDiffEqRosenbrockAlgorithm,
    OrdinaryDiffEqRosenbrockAdaptiveAlgorithm}

abstract type OrdinaryDiffEqExponentialAlgorithm{CS, AD, FDT, ST, CJ} <:
              OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD, FDT, ST, CJ} <:
              OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqLinearExponentialAlgorithm <:
              OrdinaryDiffEqExponentialAlgorithm{
    0,
    false,
    Val{:forward},
    Val{true},
    nothing
} end
const ExponentialAlgorithm = Union{OrdinaryDiffEqExponentialAlgorithm,
    OrdinaryDiffEqAdaptiveExponentialAlgorithm}

abstract type OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm <:
              OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS, AD, FDT, ST, CJ} <:
              OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD, FDT, ST, CJ} end

# DAE Specific Algorithms
abstract type DAEAlgorithm{CS, AD, FDT, ST, CJ} <: DiffEqBase.AbstractDAEAlgorithm end

# Partitioned ODE Specific Algorithms
abstract type OrdinaryDiffEqPartitionedAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqAdaptivePartitionedAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
const PartitionedAlgorithm = Union{OrdinaryDiffEqPartitionedAlgorithm,
    OrdinaryDiffEqAdaptivePartitionedAlgorithm}

struct FunctionMap{scale_by_time} <: OrdinaryDiffEqAlgorithm end
FunctionMap(; scale_by_time = false) = FunctionMap{scale_by_time}()

function DiffEqBase.remake(thing::OrdinaryDiffEqAlgorithm; kwargs...)
    T = SciMLBase.remaker_of(thing)
    T(; SciMLBase.struct_as_namedtuple(thing)..., kwargs...)
end

function DiffEqBase.remake(
        thing::Union{
            OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD, FDT,
                ST, CJ},
            OrdinaryDiffEqImplicitAlgorithm{CS, AD, FDT, ST, CJ
            },
            DAEAlgorithm{CS, AD, FDT, ST, CJ}};
        kwargs...) where {CS, AD, FDT, ST, CJ}
    T = SciMLBase.remaker_of(thing)
    T(; SciMLBase.struct_as_namedtuple(thing)...,
        chunk_size = Val{CS}(), autodiff = Val{AD}(), standardtag = Val{ST}(),
        concrete_jac = CJ === nothing ? CJ : Val{CJ}(),
        kwargs...)
end

###############################################################################

# RK methods

struct ExplicitRK{TabType} <: OrdinaryDiffEqAdaptiveAlgorithm
    tableau::TabType
end
ExplicitRK(; tableau = ODE_DEFAULT_TABLEAU) = ExplicitRK(tableau)

TruncatedStacktraces.@truncate_stacktrace ExplicitRK

@inline trivial_limiter!(u, integrator, p, t) = nothing

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

"""
    SIR54(; stage_limiter! = OrdinaryDiffEq.trivial_limiter!,
             step_limiter! = OrdinaryDiffEq.trivial_limiter!,
             thread = OrdinaryDiffEq.False())

5th order Explicit RK method suited for SIR-type epidemic models.

Like SSPRK methods, this method also takes optional arguments `stage_limiter!`
and `step_limiter!`, where `stage_limiter!` and `step_limiter!` are functions
of the form `limiter!(u, integrator, p, t)`.

The argument `thread` determines whether internal broadcasting on
appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`,
default) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when
Julia is started with multiple threads.

## Reference

@article{Kovalnogov2020RungeKuttaPS,
title={Runge–Kutta pairs suited for SIR‐type epidemic models},
author={Vladislav N. Kovalnogov and Theodore E. Simos and Ch. Tsitouras},
journal={Mathematical Methods in the Applied Sciences},
year={2020},
volume={44},
pages={5210 - 5216}
}
"""
struct SIR54{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function SIR54(; stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!,
        thread = False())
    SIR54{typeof(stage_limiter!), typeof(step_limiter!), typeof(thread)}(stage_limiter!,
        step_limiter!,
        thread)
end

# for backwards compatibility
function SIR54(stage_limiter!, step_limiter! = trivial_limiter!)
    SIR54{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

function Base.show(io::IO, alg::SIR54)
    print(io, "SIR54(stage_limiter! = ", alg.stage_limiter!,
        ", step_limiter! = ", alg.step_limiter!,
        ", thread = ", alg.thread, ")")
end

"""
    Alshina2(; stage_limiter! = OrdinaryDiffEq.trivial_limiter!,
             step_limiter! = OrdinaryDiffEq.trivial_limiter!,
             thread = OrdinaryDiffEq.False())

2nd order, 2-stage Explicit Runge-Kutta Method with optimal parameters.

Like SSPRK methods, this method also takes optional arguments `stage_limiter!`
and `step_limiter!`, where `stage_limiter!` and `step_limiter!` are functions
of the form `limiter!(u, integrator, p, t)`.

The argument `thread` determines whether internal broadcasting on
appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`,
default) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when
Julia is started with multiple threads.

## Reference

@article{Alshina2008,
doi = {10.1134/s0965542508030068},
url = {https://doi.org/10.1134/s0965542508030068},
year = {2008},
month = mar,
publisher = {Pleiades Publishing Ltd},
volume = {48},
number = {3},
pages = {395--405},
author = {E. A. Alshina and E. M. Zaks and N. N. Kalitkin},
title = {Optimal first- to sixth-order accurate Runge-Kutta schemes},
journal = {Computational Mathematics and Mathematical Physics}
}
"""
struct Alshina2{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function Alshina2(; stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!,
        thread = False())
    Alshina2{typeof(stage_limiter!), typeof(step_limiter!), typeof(thread)}(stage_limiter!,
        step_limiter!,
        thread)
end

function Alshina2(stage_limiter!, step_limiter! = trivial_limiter!)
    Alshina2{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

function Base.show(io::IO, alg::Alshina2)
    print(io, "Alshina2(stage_limiter! = ", alg.stage_limiter!,
        ", step_limiter! = ", alg.step_limiter!,
        ", thread = ", alg.thread, ")")
end

"""
    Alshina3(; stage_limiter! = OrdinaryDiffEq.trivial_limiter!,
             step_limiter! = OrdinaryDiffEq.trivial_limiter!,
             thread = OrdinaryDiffEq.False())

3rd order, 3-stage Explicit Runge-Kutta Method with optimal parameters.

Like SSPRK methods, this method also takes optional arguments `stage_limiter!`
and `step_limiter!`, where `stage_limiter!` and `step_limiter!` are functions
of the form `limiter!(u, integrator, p, t)`.

The argument `thread` determines whether internal broadcasting on
appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`,
default) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when
Julia is started with multiple threads.

## Reference

@article{Alshina2008,
doi = {10.1134/s0965542508030068},
url = {https://doi.org/10.1134/s0965542508030068},
year = {2008},
month = mar,
publisher = {Pleiades Publishing Ltd},
volume = {48},
number = {3},
pages = {395--405},
author = {E. A. Alshina and E. M. Zaks and N. N. Kalitkin},
title = {Optimal first- to sixth-order accurate Runge-Kutta schemes},
journal = {Computational Mathematics and Mathematical Physics}
}
"""
struct Alshina3{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function Alshina3(; stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!,
        thread = False())
    Alshina3{typeof(stage_limiter!), typeof(step_limiter!), typeof(thread)}(stage_limiter!,
        step_limiter!,
        thread)
end

function Alshina3(stage_limiter!, step_limiter! = trivial_limiter!)
    Alshina3{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

function Base.show(io::IO, alg::Alshina3)
    print(io, "Alshina3(stage_limiter! = ", alg.stage_limiter!,
        ", step_limiter! = ", alg.step_limiter!,
        ", thread = ", alg.thread, ")")
end

"""
    Alshina6(; stage_limiter! = OrdinaryDiffEq.trivial_limiter!,
             step_limiter! = OrdinaryDiffEq.trivial_limiter!,
             thread = OrdinaryDiffEq.False())

6th order, 7-stage Explicit Runge-Kutta Method with optimal parameters.

Like SSPRK methods, this method also takes optional arguments `stage_limiter!`
and `step_limiter!`, where `stage_limiter!` and `step_limiter!` are functions
of the form `limiter!(u, integrator, p, t)`.

The argument `thread` determines whether internal broadcasting on
appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`,
default) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when
Julia is started with multiple threads.

## Reference

@article{Alshina2008,
doi = {10.1134/s0965542508030068},
url = {https://doi.org/10.1134/s0965542508030068},
year = {2008},
month = mar,
publisher = {Pleiades Publishing Ltd},
volume = {48},
number = {3},
pages = {395--405},
author = {E. A. Alshina and E. M. Zaks and N. N. Kalitkin},
title = {Optimal first- to sixth-order accurate Runge-Kutta schemes},
journal = {Computational Mathematics and Mathematical Physics}
}
"""
struct Alshina6{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function Alshina6(; stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!,
        thread = False())
    Alshina6{typeof(stage_limiter!), typeof(step_limiter!), typeof(thread)}(stage_limiter!,
        step_limiter!,
        thread)
end

function Alshina6(stage_limiter!, step_limiter! = trivial_limiter!)
    Alshina6{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

function Base.show(io::IO, alg::Alshina6)
    print(io, "Alshina6(stage_limiter! = ", alg.stage_limiter!,
        ", step_limiter! = ", alg.step_limiter!,
        ", thread = ", alg.thread, ")")
end

################################################################################

# Symplectic methods

struct SymplecticEuler <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{verlet1967computer,
title={Computer" experiments" on classical fluids. I. Thermodynamical properties of Lennard-Jones molecules},
author={Verlet, Loup},
journal={Physical review},
volume={159},
number={1},
pages={98},
year={1967},
publisher={APS}
}
"""
struct VelocityVerlet <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{verlet1967computer,
title={Computer" experiments" on classical fluids. I. Thermodynamical properties of Lennard-Jones molecules},
author={Verlet, Loup},
journal={Physical review},
volume={159},
number={1},
pages={98},
year={1967},
publisher={APS}
}
"""
struct VerletLeapfrog <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{verlet1967computer,
title={Computer" experiments" on classical fluids. I. Thermodynamical properties of Lennard-Jones molecules},
author={Verlet, Loup},
journal={Physical review},
volume={159},
number={1},
pages={98},
year={1967},
publisher={APS}
}
"""
struct PseudoVerletLeapfrog <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{mclachlan1992accuracy,
title={The accuracy of symplectic integrators},
author={McLachlan, Robert I and Atela, Pau},
journal={Nonlinearity},
volume={5},
number={2},
pages={541},
year={1992},
publisher={IOP Publishing}
}
"""
struct McAte2 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{ruth1983canonical,
title={A canonical integration technique},
author={Ruth, Ronald D},
journal={IEEE Trans. Nucl. Sci.},
volume={30},
number={CERN-LEP-TH-83-14},
pages={2669--2671},
year={1983}
}
"""
struct Ruth3 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{mclachlan1992accuracy,
title={The accuracy of symplectic integrators},
author={McLachlan, Robert I and Atela, Pau},
journal={Nonlinearity},
volume={5},
number={2},
pages={541},
year={1992},
publisher={IOP Publishing}
}
"""
struct McAte3 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{candy1991symplectic,
title={A symplectic integration algorithm for separable Hamiltonian functions},
author={Candy, J and Rozmus, W},
journal={Journal of Computational Physics},
volume={92},
number={1},
pages={230--256},
year={1991},
publisher={Elsevier}
}
"""
struct CandyRoz4 <: OrdinaryDiffEqPartitionedAlgorithm end
struct McAte4 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{sanz1993symplectic,
title={Symplectic numerical methods for Hamiltonian problems},
author={Sanz-Serna, Jes{\'u}s Maria and Calvo, Mari-Paz},
journal={International Journal of Modern Physics C},
volume={4},
number={02},
pages={385--392},
year={1993},
publisher={World Scientific}
}
"""
struct CalvoSanz4 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{mclachlan1992accuracy,
title={The accuracy of symplectic integrators},
author={McLachlan, Robert I and Atela, Pau},
journal={Nonlinearity},
volume={5},
number={2},
pages={541},
year={1992},
publisher={IOP Publishing}
}
"""
struct McAte42 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{mclachlan1992accuracy,
title={The accuracy of symplectic integrators},
author={McLachlan, Robert I and Atela, Pau},
journal={Nonlinearity},
volume={5},
number={2},
pages={541},
year={1992},
publisher={IOP Publishing}
}
"""
struct McAte5 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{yoshida1990construction,
title={Construction of higher order symplectic integrators},
author={Yoshida, Haruo},
journal={Physics letters A},
volume={150},
number={5-7},
pages={262--268},
year={1990},
publisher={Elsevier}
}
"""
struct Yoshida6 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{kahan1997composition,
title={Composition constants for raising the orders of unconventional schemes for ordinary differential equations},
author={Kahan, William and Li, Ren-Cang},
journal={Mathematics of computation},
volume={66},
number={219},
pages={1089--1099},
year={1997}
}
"""
struct KahanLi6 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{mclachlan1995numerical,
title={On the numerical integration of ordinary differential equations by symmetric composition methods},
author={McLachlan, Robert I},
journal={SIAM Journal on Scientific Computing},
volume={16},
number={1},
pages={151--168},
year={1995},
publisher={SIAM}
}
"""
struct McAte8 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{kahan1997composition,
title={Composition constants for raising the orders of unconventional schemes for ordinary differential equations},
author={Kahan, William and Li, Ren-Cang},
journal={Mathematics of computation},
volume={66},
number={219},
pages={1089--1099},
year={1997}
}
"""
struct KahanLi8 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
@article{sofroniou2005derivation,
title={Derivation of symmetric composition constants for symmetric integrators},
author={Sofroniou, Mark and Spaletta, Giulia},
journal={Optimization Methods and Software},
volume={20},
number={4-5},
pages={597--613},
year={2005},
publisher={Taylor \\& Francis}
}
"""
struct SofSpa10 <: OrdinaryDiffEqPartitionedAlgorithm end

# Nyström methods

"""
    IRKN3

Improved Runge-Kutta-Nyström method of order three, which minimizes the amount of evaluated functions in each step. Fixed time steps only.

Second order ODE should not depend on the first derivative.

## References

@article{rabiei2012numerical,
title={Numerical Solution of Second-Order Ordinary Differential Equations by Improved Runge-Kutta Nystrom Method},
author={Rabiei, Faranak and Ismail, Fudziah and Norazak, S and Emadi, Saeid},
publisher={Citeseer}
}
"""
struct IRKN3 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
    Nystrom4

A 4th order explicit Runge-Kutta-Nyström method which can be applied directly on second order ODEs. Can only be used with fixed time steps.

In case the ODE Problem is not dependent on the first derivative consider using
[`Nystrom4VelocityIndependent`](@ref) to increase performance.

## References

E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
Springer-Verlag.
"""
struct Nystrom4 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
    FineRKN4()

A 4th order explicit Runge-Kutta-Nyström method which can be applied directly to second order ODEs.
In particular, this method allows the acceleration equation to depend on the velocity.

## References

```
@article{fine1987low,
  title={Low order practical {R}unge-{K}utta-{N}ystr{\"o}m methods},
  author={Fine, Jerry Michael},
  journal={Computing},
  volume={38},
  number={4},
  pages={281--297},
  year={1987},
  publisher={Springer}
}
```
"""
struct FineRKN4 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    FineRKN5()

A 5th order explicit Runge-Kutta-Nyström method which can be applied directly to second order ODEs.
In particular, this method allows the acceleration equation to depend on the velocity.

## References

```
@article{fine1987low,
  title={Low order practical {R}unge-{K}utta-{N}ystr{\"o}m methods},
  author={Fine, Jerry Michael},
  journal={Computing},
  volume={38},
  number={4},
  pages={281--297},
  year={1987},
  publisher={Springer}
}
```
"""
struct FineRKN5 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    Nystrom4VelocityIdependent

A 4th order explicit Runkge-Kutta-Nyström method. Used directly on second order ODEs, where the acceleration is independent from velocity (ODE Problem is not dependent on the first derivative).

More efficient then [`Nystrom4`](@ref) on velocity independent problems, since less evaluations are needed.

Fixed time steps only.

## References

E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
Springer-Verlag.
"""
struct Nystrom4VelocityIndependent <: OrdinaryDiffEqPartitionedAlgorithm end

"""
    IRKN4

Improves Runge-Kutta-Nyström method of order four, which minimizes the amount of evaluated functions in each step. Fixed time steps only.

Second order ODE should not be dependent on the first derivative.

Recommended for smooth problems with expensive functions to evaluate.

## References

@article{rabiei2012numerical,
title={Numerical Solution of Second-Order Ordinary Differential Equations by Improved Runge-Kutta Nystrom Method},
author={Rabiei, Faranak and Ismail, Fudziah and Norazak, S and Emadi, Saeid},
publisher={Citeseer}
}
"""
struct IRKN4 <: OrdinaryDiffEqPartitionedAlgorithm end

"""
    Nystrom5VelocityIndependent

A 5th order explicit Runkge-Kutta-Nyström method. Used directly on second order ODEs, where the acceleration is independent from velocity (ODE Problem is not dependent on the first derivative).
Fixed time steps only.

## References

E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
Springer-Verlag.
"""
struct Nystrom5VelocityIndependent <: OrdinaryDiffEqPartitionedAlgorithm end

"""
    DPRKN4

4th order explicit Runge-Kutta-Nyström methods. The second order ODE should not depend on the first derivative.

## References

@article{Dormand1987FamiliesOR,
title={Families of Runge-Kutta-Nystrom Formulae},
author={J. R. Dormand and Moawwad E. A. El-Mikkawy and P. J. Prince},
journal={Ima Journal of Numerical Analysis},
year={1987},
volume={7},
pages={235-250}
}
"""
struct DPRKN4 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    DPRKN5

5th order explicit Runge-Kutta-Nyström mehod. The second order ODE should not depend on the first derivative.

## References

@article{Bettis1973ARN,
title={A Runge-Kutta Nystrom algorithm},
author={Dale G. Bettis},
journal={Celestial mechanics},
year={1973},
volume={8},
pages={229-233},
publisher={Springer}
}
"""
struct DPRKN5 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    DPRKN6

6th order explicit Runge-Kutta-Nyström method. The second order ODE should not depend on the first derivative. Free 6th order interpolant.

## References

@article{dormand1987runge,
title={Runge-kutta-nystrom triples},
author={Dormand, JR and Prince, PJ},
journal={Computers \\& Mathematics with Applications},
volume={13},
number={12},
pages={937--949},
year={1987},
publisher={Elsevier}
}
"""
struct DPRKN6 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    DPRKN6FM

6th order explicit Runge-Kutta-Nyström method. The second order ODE should not depend on the first derivative.

Compared to [`DPRKN6`](@ref), this method has smaller truncation error coefficients which leads to performance gain
when only the main solution points are considered.

## References

@article{Dormand1987FamiliesOR,
title={Families of Runge-Kutta-Nystrom Formulae},
author={J. R. Dormand and Moawwad E. A. El-Mikkawy and P. J. Prince},
journal={Ima Journal of Numerical Analysis},
year={1987},
volume={7},
pages={235-250}
}
"""
struct DPRKN6FM <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    DPRKN8

8th order explicit Runge-Kutta-Nyström method. The second order ODE should not depend on the first derivative.

Not as efficient as [`DPRKN12`](@ref) when high accuracy is needed, however this solver is competitive with
[`DPRKN6`](@ref) at lax tolerances and, depending on the problem, might be a good option between performance and accuracy.

## References

@article{dormand1987high,
title={High-order embedded Runge-Kutta-Nystrom formulae},
author={Dormand, JR and El-Mikkawy, MEA and Prince, PJ},
journal={IMA Journal of Numerical Analysis},
volume={7},
number={4},
pages={423--430},
year={1987},
publisher={Oxford University Press}
}
"""
struct DPRKN8 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    DPRKN12

12th order explicit Rugne-Kutta-Nyström method. The second order ODE should not depend on the first derivative.

Most efficient when high accuracy is needed.

## References

@article{dormand1987high,
title={High-order embedded Runge-Kutta-Nystrom formulae},
author={Dormand, JR and El-Mikkawy, MEA and Prince, PJ},
journal={IMA Journal of Numerical Analysis},
volume={7},
number={4},
pages={423--430},
year={1987},
publisher={Oxford University Press}
}
"""
struct DPRKN12 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    ERKN4

Embedded 4(3) pair of explicit Runge-Kutta-Nyström methods. Integrates the periodic properties of the harmonic oscillator exactly.

The second order ODE should not depend on the first derivative.

Uses adaptive step size control. This method is extra efficient on periodic problems.

## References

@article{demba2017embedded,
title={An Embedded 4 (3) Pair of Explicit Trigonometrically-Fitted Runge-Kutta-Nystr{\"o}m Method for Solving Periodic Initial Value Problems},
author={Demba, MA and Senu, N and Ismail, F},
journal={Applied Mathematical Sciences},
volume={11},
number={17},
pages={819--838},
year={2017}
}
"""
struct ERKN4 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    ERKN5

Embedded 5(4) pair of explicit Runge-Kutta-Nyström methods. Integrates the periodic properties of the harmonic oscillator exactly.

The second order ODE should not depend on the first derivative.

Uses adaptive step size control. This method is extra efficient on periodic problems.

## References

@article{demba20165,
title={A 5 (4) Embedded Pair of Explicit Trigonometrically-Fitted Runge--Kutta--Nystr{\"o}m Methods for the Numerical Solution of Oscillatory Initial Value Problems},
author={Demba, Musa A and Senu, Norazak and Ismail, Fudziah},
journal={Mathematical and Computational Applications},
volume={21},
number={4},
pages={46},
year={2016},
publisher={Multidisciplinary Digital Publishing Institute}
}
"""
struct ERKN5 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
    ERKN7

Embedded pair of explicit Runge-Kutta-Nyström methods. Integrates the periodic properties of the harmonic oscillator exactly.

The second order ODE should not depend on the first derivative.

Uses adaptive step size control. This method is extra efficient on periodic Problems.

## References

@article{SimosOnHO,
title={On high order Runge-Kutta-Nystr{\"o}m pairs},
author={Theodore E. Simos and Ch. Tsitouras},
journal={J. Comput. Appl. Math.},
volume={400},
pages={113753}
}
"""
struct ERKN7 <: OrdinaryDiffEqAdaptivePartitionedAlgorithm end

"""
3 stage fourth order Runge-Kutta Nystrom method to solve second order linear inhomogenous IVPs.

Does not include an adaptive method. Solves for for d-dimensional differential systems of second order linear inhomogeneous equations.

!!! warn

    This method is only fourth order for these systems, the method is second order otherwise!

## References

@article{MONTIJANO2024115533,
title = {Explicit Runge–Kutta–Nyström methods for the numerical solution of second order linear inhomogeneous IVPs},
author = {J.I. Montijano and L. Rández and M. Calvo},
journal = {Journal of Computational and Applied Mathematics},
volume = {438},
pages = {115533},
year = {2024},
}
"""
struct RKN4 <: OrdinaryDiffEqAlgorithm end

################################################################################

# Adams Bashforth and Adams moulton methods

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

AB3: Adams-Bashforth Explicit Method
The 3-step third order multistep method. Ralston's Second Order Method is used to calculate starting values.
"""
struct AB3 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

AB4: Adams-Bashforth Explicit Method
The 4-step fourth order multistep method. Runge-Kutta method of order 4 is used to calculate starting values.
"""
struct AB4 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

AB5: Adams-Bashforth Explicit Method
The 3-step third order multistep method. Ralston's Second Order Method is used to calculate starting values.
"""
struct AB5 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

ABM32: Adams-Bashforth Explicit Method
It is third order method. In ABM32, AB3 works as predictor and Adams Moulton 2-steps method works as Corrector.
Ralston's Second Order Method is used to calculate starting values.
"""
struct ABM32 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

ABM43: Adams-Bashforth Explicit Method
It is fourth order method. In ABM43, AB4 works as predictor and Adams Moulton 3-steps method works as Corrector.
Runge-Kutta method of order 4 is used to calculate starting values.
"""
struct ABM43 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

ABM54: Adams-Bashforth Explicit Method
It is fifth order method. In ABM54, AB5 works as predictor and Adams Moulton 4-steps method works as Corrector.
Runge-Kutta method of order 4 is used to calculate starting values.
"""
struct ABM54 <: OrdinaryDiffEqAlgorithm end

# Variable Step Size Adams methods

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCAB3: Adaptive step size Adams explicit Method
The 3rd order Adams method. Bogacki-Shampine 3/2 method is used to calculate starting values.
"""
struct VCAB3 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCAB4: Adaptive step size Adams explicit Method
The 4th order Adams method. Runge-Kutta 4 is used to calculate starting values.
"""
struct VCAB4 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCAB5: Adaptive step size Adams explicit Method
The 5th order Adams method. Runge-Kutta 4 is used to calculate starting values.
"""
struct VCAB5 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCABM3: Adaptive step size Adams explicit Method
The 3rd order Adams-Moulton method. Bogacki-Shampine 3/2 method is used to calculate starting values.
"""
struct VCABM3 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCABM4: Adaptive step size Adams explicit Method
The 4th order Adams-Moulton method. Runge-Kutta 4 is used to calculate starting values.
"""
struct VCABM4 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCABM5: Adaptive step size Adams explicit Method
The 5th order Adams-Moulton method. Runge-Kutta 4 is used to calculate starting values.
"""
struct VCABM5 <: OrdinaryDiffEqAdaptiveAlgorithm end

# Variable Order and Variable Step Size Adams methods

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1

VCABM: Adaptive step size Adams explicit Method
An adaptive order adaptive time Adams Moulton method.
It uses an order adaptivity algorithm is derived from Shampine's DDEABM.
"""
struct VCABM <: OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm end

# IMEX Multistep methods

struct CNAB2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end

function CNAB2(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    CNAB2{
        _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(
        linsolve,
        nlsolve,
        precs,
        extrapolant)
end

struct CNLF2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end
function CNLF2(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    CNLF2{
        _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(
        linsolve,
        nlsolve,
        precs,
        extrapolant)
end

"""
QNDF1: Multistep Method
An adaptive order 1 quasi-constant timestep L-stable numerical differentiation function (NDF) method.
Optional parameter kappa defaults to Shampine's accuracy-optimal -0.1850.

See also `QNDF`.
"""
struct QNDF1{CS, AD, F, F2, P, FDT, ST, CJ, κType, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    kappa::κType
    controller::Symbol
    step_limiter!::StepLimiter
end

function QNDF1(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, kappa = -0.1850,
        controller = :Standard, step_limiter! = trivial_limiter!)
    QNDF1{
        _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(kappa), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        kappa,
        controller,
        step_limiter!)
end

"""
QBDF1: Multistep Method

An alias of `QNDF1` with κ=0.
"""
QBDF1(; kwargs...) = QNDF1(; kappa = 0, kwargs...)

"""
QNDF2: Multistep Method
An adaptive order 2 quasi-constant timestep L-stable numerical differentiation function (NDF) method.

See also `QNDF`.
"""
struct QNDF2{CS, AD, F, F2, P, FDT, ST, CJ, κType, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    kappa::κType
    controller::Symbol
    step_limiter!::StepLimiter
end

function QNDF2(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, kappa = -1 // 9,
        controller = :Standard, step_limiter! = trivial_limiter!)
    QNDF2{
        _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(kappa), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        kappa,
        controller,
        step_limiter!)
end

"""
QBDF2: Multistep Method

An alias of `QNDF2` with κ=0.
"""
QBDF2(; kwargs...) = QNDF2(; kappa = 0, kwargs...)

"""
QNDF: Multistep Method
An adaptive order quasi-constant timestep NDF method.
Utilizes Shampine's accuracy-optimal kappa values as defaults (has a keyword argument for a tuple of kappa coefficients).

@article{shampine1997matlab,
title={The matlab ode suite},
author={Shampine, Lawrence F and Reichelt, Mark W},
journal={SIAM journal on scientific computing},
volume={18},
number={1},
pages={1--22},
year={1997},
publisher={SIAM}
}
"""
struct QNDF{MO, CS, AD, F, F2, P, FDT, ST, CJ, K, T, κType, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    kappa::κType
    controller::Symbol
    step_limiter!::StepLimiter
end

function QNDF(; max_order::Val{MO} = Val{5}(), chunk_size = Val{0}(),
        autodiff = Val{true}(), standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, kappa = promote(-0.1850, -1 // 9, -0.0823, -0.0415, 0),
        controller = :Standard, step_limiter! = trivial_limiter!) where {MO}
    QNDF{MO, _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol), typeof(kappa), typeof(step_limiter!)}(
        max_order, linsolve, nlsolve, precs, κ, tol,
        extrapolant, kappa, controller, step_limiter!)
end

TruncatedStacktraces.@truncate_stacktrace QNDF

"""
QBDF: Multistep Method

An alias of `QNDF` with κ=0.
"""
QBDF(; kwargs...) = QNDF(; kappa = tuple(0 // 1, 0 // 1, 0 // 1, 0 // 1, 0 // 1), kwargs...)

"""
FBDF: Fixed leading coefficient BDF

An adaptive order quasi-constant timestep NDF method.
Utilizes Shampine's accuracy-optimal kappa values as defaults (has a keyword argument for a tuple of kappa coefficients).

@article{shampine2002solving,
title={Solving 0= F (t, y (t), y′(t)) in Matlab},
author={Shampine, Lawrence F},
year={2002},
publisher={Walter de Gruyter GmbH \\& Co. KG}
}
"""
struct FBDF{MO, CS, AD, F, F2, P, FDT, ST, CJ, K, T, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end

function FBDF(; max_order::Val{MO} = Val{5}(), chunk_size = Val{0}(),
        autodiff = Val{true}(), standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, controller = :Standard, step_limiter! = trivial_limiter!) where {MO}
    FBDF{MO, _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol), typeof(step_limiter!)}(
        max_order, linsolve, nlsolve, precs, κ, tol, extrapolant,
        controller, step_limiter!)
end

TruncatedStacktraces.@truncate_stacktrace FBDF

"""
Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton. Implicit-Explicit Methods for Time-
Dependent Partial Differential Equations. 1995 Society for Industrial and Applied Mathematics
Journal on Numerical Analysis, 32(3), pp 797-823, 1995. doi: https://doi.org/10.1137/0732037
"""
struct SBDF{CS, AD, F, F2, P, FDT, ST, CJ, K, T} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    order::Int
    ark::Bool
end

function SBDF(order; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, ark = false)
    SBDF{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol)}(linsolve,
        nlsolve,
        precs,
        κ,
        tol,
        extrapolant,
        order,
        ark)
end

# All keyword form needed for remake
function SBDF(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear,
        order, ark = false)
    SBDF{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol)}(linsolve,
        nlsolve,
        precs,
        κ,
        tol,
        extrapolant,
        order,
        ark)
end

"""
    IMEXEuler(;kwargs...)

The one-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

When applied to a `SplitODEProblem` of the form

```
u'(t) = f1(u) + f2(u)
```

The default `IMEXEuler()` method uses an update of the form

```
unew = uold + dt * (f1(unew) + f2(uold))
```

See also `SBDF`, `IMEXEulerARK`.
"""
IMEXEuler(; kwargs...) = SBDF(1; kwargs...)

"""
    IMEXEulerARK(;kwargs...)

The one-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

When applied to a `SplitODEProblem` of the form

```
u'(t) = f1(u) + f2(u)
```

A classical additive Runge-Kutta method in the sense of
[Araújo, Murua, Sanz-Serna (1997)](https://doi.org/10.1137/S0036142995292128)
consisting of the implicit and the explicit Euler method given by

```
y1   = uold + dt * f1(y1)
unew = uold + dt * (f1(unew) + f2(y1))
```

See also `SBDF`, `IMEXEuler`.
"""
IMEXEulerARK(; kwargs...) = SBDF(1; ark = true, kwargs...)

"""
    SBDF2(;kwargs...)

The two-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

See also `SBDF`.
"""
SBDF2(; kwargs...) = SBDF(2; kwargs...)

"""
    SBDF3(;kwargs...)

The three-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

See also `SBDF`.
"""
SBDF3(; kwargs...) = SBDF(3; kwargs...)

"""
    SBDF4(;kwargs...)

The four-step version of the IMEX multistep methods of

  - Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton.
    Implicit-Explicit Methods for Time-Dependent Partial Differential Equations.
    Society for Industrial and Applied Mathematics.
    Journal on Numerical Analysis, 32(3), pp 797-823, 1995.
    doi: [https://doi.org/10.1137/0732037](https://doi.org/10.1137/0732037)

See also `SBDF`.
"""
SBDF4(; kwargs...) = SBDF(4; kwargs...)

# Adams/BDF methods in Nordsieck forms
"""
AN5: Adaptive step size Adams explicit Method
An adaptive 5th order fixed-leading coefficient Adams method in Nordsieck form.
"""
struct AN5 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct JVODE{bType, aType} <: OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm
    algorithm::Symbol
    bias1::bType
    bias2::bType
    bias3::bType
    addon::aType
end

function JVODE(algorithm = :Adams; bias1 = 6, bias2 = 6, bias3 = 10,
        addon = 1 // 10^6)
    JVODE(algorithm, bias1, bias2, bias3, addon)
end
JVODE_Adams(; kwargs...) = JVODE(:Adams; kwargs...)
JVODE_BDF(; kwargs...) = JVODE(:BDF; kwargs...)

# ROCK methods

"""
Assyr Abdulle, Alexei A. Medovikov. Second Order Chebyshev Methods based on Orthogonal Polynomials.
Numerische Mathematik, 90 (1), pp 1-18, 2001. doi: https://dx.doi.org/10.1007/s002110100292

ROCK2: Stabilized Explicit Method.
Second order stabilized Runge-Kutta method.
Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.

This method takes optional keyword arguments `min_stages`, `max_stages`, and `eigen_est`.
The function `eigen_est` should be of the form

`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`,

where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix. If `eigen_est`
is not provided, `upper_bound` will be estimated using the power iteration.
"""
struct ROCK2{E} <: OrdinaryDiffEqAdaptiveAlgorithm
    min_stages::Int
    max_stages::Int
    eigen_est::E
end
function ROCK2(; min_stages = 0, max_stages = 200, eigen_est = nothing)
    ROCK2(min_stages, max_stages, eigen_est)
end

"""
    ROCK4(; min_stages = 0, max_stages = 152, eigen_est = nothing)

Assyr Abdulle. Fourth Order Chebyshev Methods With Recurrence Relation. 2002 Society for
Industrial and Applied Mathematics Journal on Scientific Computing, 23(6), pp 2041-2054, 2001.
doi: https://doi.org/10.1137/S1064827500379549

ROCK4: Stabilized Explicit Method.
Fourth order stabilized Runge-Kutta method.
Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.

This method takes optional keyword arguments `min_stages`, `max_stages`, and `eigen_est`.
The function `eigen_est` should be of the form

`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`,

where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix. If `eigen_est`
is not provided, `upper_bound` will be estimated using the power iteration.
"""
struct ROCK4{E} <: OrdinaryDiffEqAdaptiveAlgorithm
    min_stages::Int
    max_stages::Int
    eigen_est::E
end
function ROCK4(; min_stages = 0, max_stages = 152, eigen_est = nothing)
    ROCK4(min_stages, max_stages, eigen_est)
end

# SERK methods

for Alg in [:ESERK4, :ESERK5, :RKC]
    @eval begin
        struct $Alg{E} <: OrdinaryDiffEqAdaptiveAlgorithm
            eigen_est::E
        end
        $Alg(; eigen_est = nothing) = $Alg(eigen_est)
    end
end

"""
    RKC(; eigen_est = nothing)

B. P. Sommeijer, L. F. Shampine, J. G. Verwer. RKC: An Explicit Solver for Parabolic PDEs,
Journal of Computational and Applied Mathematics, 88(2), pp 315-326, 1998. doi:
https://doi.org/10.1016/S0377-0427(97)00219-7

RKC: Stabilized Explicit Method.
Second order stabilized Runge-Kutta method.
Exhibits high stability for real eigenvalues.

This method takes the keyword argument `eigen_est` of the form

`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`,

where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix. If `eigen_est`
is not provided, `upper_bound` will be estimated using the power iteration.
"""
function RKC end

"""
    ESERK4(; eigen_est = nothing)

J. Martín-Vaquero, B. Kleefeld. Extrapolated stabilized explicit Runge-Kutta methods,
Journal of Computational Physics, 326, pp 141-155, 2016. doi:
https://doi.org/10.1016/j.jcp.2016.08.042.

ESERK4: Stabilized Explicit Method.
Fourth order extrapolated stabilized Runge-Kutta method.
Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.

This method takes the keyword argument `eigen_est` of the form

`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`,

where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.
"""
function ESERK4 end

"""
    ESERK5(; eigen_est = nothing)

J. Martín-Vaquero, A. Kleefeld. ESERK5: A fifth-order extrapolated stabilized explicit Runge-Kutta method,
Journal of Computational and Applied Mathematics, 356, pp 22-36, 2019. doi:
https://doi.org/10.1016/j.cam.2019.01.040.

ESERK5: Stabilized Explicit Method.
Fifth order extrapolated stabilized Runge-Kutta method.
Exhibits high stability for real eigenvalues and is smoothened to allow for moderate sized complex eigenvalues.

This method takes the keyword argument `eigen_est` of the form

`eigen_est = (integrator) -> integrator.eigen_est = upper_bound`,

where `upper_bound` is an estimated upper bound on the spectral radius of the Jacobian matrix.
If `eigen_est` is not provided, `upper_bound` will be estimated using the power iteration.
"""
function ESERK5 end

struct SERK2{E} <: OrdinaryDiffEqAdaptiveAlgorithm
    controller::Symbol
    eigen_est::E
end
SERK2(; controller = :PI, eigen_est = nothing) = SERK2(controller, eigen_est)

struct IRKC{CS, AD, F, F2, P, FDT, ST, CJ, K, T, E} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    controller::Symbol
    eigen_est::E
end

function IRKC(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, controller = :Standard, eigen_est = nothing)
    IRKC{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol), typeof(eigen_est)}(linsolve, nlsolve, precs, κ, tol,
        extrapolant, controller, eigen_est)
end

################################################################################

# Linear Methods

for Alg in [
    :MagnusMidpoint,
    :MagnusLeapfrog,
    :LieEuler,
    :MagnusGauss4,
    :MagnusNC6,
    :MagnusGL6,
    :MagnusGL8,
    :MagnusNC8,
    :MagnusGL4,
    :RKMK2,
    :RKMK4,
    :LieRK4,
    :CG2,
    :CG3,
    :CG4a
]
    @eval struct $Alg <: OrdinaryDiffEqLinearExponentialAlgorithm
        krylov::Bool
        m::Int
        iop::Int
    end
    @eval $Alg(; krylov = false, m = 30, iop = 0) = $Alg(krylov, m, iop)
end

struct MagnusAdapt4 <: OrdinaryDiffEqAdaptiveAlgorithm end

struct LinearExponential <:
       OrdinaryDiffEqExponentialAlgorithm{1, false, Val{:forward}, Val{true}, nothing}
    krylov::Symbol
    m::Int
    iop::Int
end
LinearExponential(; krylov = :off, m = 10, iop = 0) = LinearExponential(krylov, m, iop)

struct CayleyEuler <: OrdinaryDiffEqAlgorithm end

################################################################################

# FIRK Methods

"""
@article{hairer1999stiff,
title={Stiff differential equations solved by Radau methods},
author={Hairer, Ernst and Wanner, Gerhard},
journal={Journal of Computational and Applied Mathematics},
volume={111},
number={1-2},
pages={93--111},
year={1999},
publisher={Elsevier}
}

RadauIIA3: Fully-Implicit Runge-Kutta Method
An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency.
"""
struct RadauIIA3{CS, AD, F, P, FDT, ST, CJ, Tol, C1, C2, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    extrapolant::Symbol
    κ::Tol
    maxiters::Int
    fast_convergence_cutoff::C1
    new_W_γdt_cutoff::C2
    controller::Symbol
    step_limiter!::StepLimiter
end

function RadauIIA3(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10,
        step_limiter! = trivial_limiter!)
    RadauIIA3{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(fast_convergence_cutoff),
        typeof(new_W_γdt_cutoff), typeof(step_limiter!)}(linsolve,
        precs,
        extrapolant,
        κ,
        maxiters,
        fast_convergence_cutoff,
        new_W_γdt_cutoff,
        controller,
        step_limiter!)
end

TruncatedStacktraces.@truncate_stacktrace RadauIIA3

"""
@article{hairer1999stiff,
title={Stiff differential equations solved by Radau methods},
author={Hairer, Ernst and Wanner, Gerhard},
journal={Journal of Computational and Applied Mathematics},
volume={111},
number={1-2},
pages={93--111},
year={1999},
publisher={Elsevier}
}

RadauIIA5: Fully-Implicit Runge-Kutta Method
An A-B-L stable fully implicit Runge-Kutta method with internal tableau complex basis transform for efficiency.
"""
struct RadauIIA5{CS, AD, F, P, FDT, ST, CJ, Tol, C1, C2, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    κ::Tol
    maxiters::Int
    fast_convergence_cutoff::C1
    new_W_γdt_cutoff::C2
    controller::Symbol
    step_limiter!::StepLimiter
end

function RadauIIA5(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS,
        extrapolant = :dense, fast_convergence_cutoff = 1 // 5,
        new_W_γdt_cutoff = 1 // 5,
        controller = :Predictive, κ = nothing, maxiters = 10, smooth_est = true,
        step_limiter! = trivial_limiter!)
    RadauIIA5{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(fast_convergence_cutoff),
        typeof(new_W_γdt_cutoff), typeof(step_limiter!)}(linsolve,
        precs,
        smooth_est,
        extrapolant,
        κ,
        maxiters,
        fast_convergence_cutoff,
        new_W_γdt_cutoff,
        controller,
        step_limiter!)
end
TruncatedStacktraces.@truncate_stacktrace RadauIIA5

################################################################################

# SDIRK Methods
"""
ImplicitEuler: SDIRK Method
A 1st order implicit solver. A-B-L-stable. Adaptive timestepping through a divided differences estimate via memory.
Strong-stability preserving (SSP).
"""
struct ImplicitEuler{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end

function ImplicitEuler(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant,
        controller = :PI, step_limiter! = trivial_limiter!)
    ImplicitEuler{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve, precs, extrapolant, controller, step_limiter!)
end
"""
ImplicitMidpoint: SDIRK Method
A second order A-stable symplectic and symmetric implicit solver.
Good for highly stiff equations which need symplectic integration.
"""
struct ImplicitMidpoint{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    step_limiter!::StepLimiter
end

function ImplicitMidpoint(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, step_limiter! = trivial_limiter!)
    ImplicitMidpoint{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        step_limiter!)
end

"""
Andre Vladimirescu. 1994. The Spice Book. John Wiley & Sons, Inc., New York,
NY, USA.

Trapezoid: SDIRK Method
A second order A-stable symmetric ESDIRK method.
"Almost symplectic" without numerical dampening.
Also known as Crank-Nicolson when applied to PDEs. Adaptive timestepping via divided
differences approximation to the second derivative terms in the local truncation error
estimate (the SPICE approximation strategy).
"""
struct Trapezoid{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end

function Trapezoid(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    Trapezoid{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        controller,
        step_limiter!)
end

"""
@article{hosea1996analysis,
title={Analysis and implementation of TR-BDF2},
author={Hosea, ME and Shampine, LF},
journal={Applied Numerical Mathematics},
volume={20},
number={1-2},
pages={21--37},
year={1996},
publisher={Elsevier}
}

TRBDF2: SDIRK Method
A second order A-B-L-S-stable one-step ESDIRK method.
Includes stiffness-robust error estimates for accurate adaptive timestepping, smoothed derivatives for highly stiff and oscillatory problems.
"""
struct TRBDF2{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end

function TRBDF2(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    TRBDF2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

TruncatedStacktraces.@truncate_stacktrace TRBDF2

"""
@article{hindmarsh2005sundials,
title={{SUNDIALS}: Suite of nonlinear and differential/algebraic equation solvers},
author={Hindmarsh, Alan C and Brown, Peter N and Grant, Keith E and Lee, Steven L and Serban, Radu and Shumaker, Dan E and Woodward, Carol S},
journal={ACM Transactions on Mathematical Software (TOMS)},
volume={31},
number={3},
pages={363--396},
year={2005},
publisher={ACM}
}

SDIRK2: SDIRK Method
An A-B-L stable 2nd order SDIRK method
"""
struct SDIRK2{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end

function SDIRK2(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    SDIRK2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(
        linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller,
        step_limiter!)
end

struct SDIRK22{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end

function SDIRK22(;
        chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    Trapezoid{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
        nlsolve,
        precs,
        extrapolant,
        controller,
        step_limiter!)
end

struct SSPSDIRK2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ} # Not adaptive
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
end

function SSPSDIRK2(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :constant,
        controller = :PI)
    SSPSDIRK2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller)
end

"""
@article{kvaerno2004singly,
title={Singly diagonally implicit Runge--Kutta methods with an explicit first stage},
author={Kv{\\ae}rn{\\o}, Anne},
journal={BIT Numerical Mathematics},
volume={44},
number={3},
pages={489--502},
year={2004},
publisher={Springer}
}

Kvaerno3: SDIRK Method
An A-L stable stiffly-accurate 3rd order ESDIRK method
"""
struct Kvaerno3{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function Kvaerno3(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    Kvaerno3{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

"""
@book{kennedy2001additive,
title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
author={Kennedy, Christopher Alan},
year={2001},
publisher={National Aeronautics and Space Administration, Langley Research Center}
}

KenCarp3: SDIRK Method
An A-L stable stiffly-accurate 3rd order ESDIRK method with splitting
"""
struct KenCarp3{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function KenCarp3(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    KenCarp3{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

struct CFNLIRK3{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end
function CFNLIRK3(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    CFNLIRK3{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

"""
@article{hindmarsh2005sundials,
title={{SUNDIALS}: Suite of nonlinear and differential/algebraic equation solvers},
author={Hindmarsh, Alan C and Brown, Peter N and Grant, Keith E and Lee, Steven L and Serban, Radu and Shumaker, Dan E and Woodward, Carol S},
journal={ACM Transactions on Mathematical Software (TOMS)},
volume={31},
number={3},
pages={363--396},
year={2005},
publisher={ACM}
}

Cash4: SDIRK Method
An A-L stable 4th order SDIRK method
"""
struct Cash4{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    embedding::Int
    controller::Symbol
end
function Cash4(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, embedding = 3)
    Cash4{
        _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(
        linsolve,
        nlsolve,
        precs,
        smooth_est,
        extrapolant,
        embedding,
        controller)
end

struct SFSDIRK4{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end
function SFSDIRK4(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    SFSDIRK4{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

struct SFSDIRK5{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end

function SFSDIRK5(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    SFSDIRK5{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

struct SFSDIRK6{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end

function SFSDIRK6(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    SFSDIRK6{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

struct SFSDIRK7{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end

function SFSDIRK7(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    SFSDIRK7{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

struct SFSDIRK8{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end

function SFSDIRK8(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear)
    SFSDIRK8{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

"""
E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
differential-algebraic problems. Computational mathematics (2nd revised ed.),
Springer (1996)

Hairer4: SDIRK Method
An A-L stable 4th order SDIRK method
"""
struct Hairer4{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
end
function Hairer4(;
        chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    Hairer4{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller)
end

"""
E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
differential-algebraic problems. Computational mathematics (2nd revised ed.),
Springer (1996)

Hairer42: SDIRK Method
An A-L stable 4th order SDIRK method
"""
struct Hairer42{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
end
function Hairer42(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    Hairer42{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller)
end

"""
@article{kvaerno2004singly,
title={Singly diagonally implicit Runge--Kutta methods with an explicit first stage},
author={Kv{\\ae}rn{\\o}, Anne},
journal={BIT Numerical Mathematics},
volume={44},
number={3},
pages={489--502},
year={2004},
publisher={Springer}
}

Kvaerno4: SDIRK Method
An A-L stable stiffly-accurate 4th order ESDIRK method.
"""
struct Kvaerno4{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function Kvaerno4(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    Kvaerno4{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

"""
@article{kvaerno2004singly,
title={Singly diagonally implicit Runge--Kutta methods with an explicit first stage},
author={Kv{\\ae}rn{\\o}, Anne},
journal={BIT Numerical Mathematics},
volume={44},
number={3},
pages={489--502},
year={2004},
publisher={Springer}
}

Kvaerno5: SDIRK Method
An A-L stable stiffly-accurate 5th order ESDIRK method
"""
struct Kvaerno5{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function Kvaerno5(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    Kvaerno5{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

"""
@book{kennedy2001additive,
title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
author={Kennedy, Christopher Alan},
year={2001},
publisher={National Aeronautics and Space Administration, Langley Research Center}
}

KenCarp4: SDIRK Method
An A-L stable stiffly-accurate 4th order ESDIRK method with splitting
"""
struct KenCarp4{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function KenCarp4(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    KenCarp4{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end

TruncatedStacktraces.@truncate_stacktrace KenCarp4

"""
@article{kennedy2019higher,
title={Higher-order additive Runge--Kutta schemes for ordinary differential equations},
author={Kennedy, Christopher A and Carpenter, Mark H},
journal={Applied Numerical Mathematics},
volume={136},
pages={183--205},
year={2019},
publisher={Elsevier}
}

KenCarp47: SDIRK Method
An A-L stable stiffly-accurate 4th order seven-stage ESDIRK method with splitting
"""
struct KenCarp47{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
end
function KenCarp47(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    KenCarp47{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller)
end

"""
@book{kennedy2001additive,
title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
author={Kennedy, Christopher Alan},
year={2001},
publisher={National Aeronautics and Space Administration, Langley Research Center}
}

KenCarp5: SDIRK Method
An A-L stable stiffly-accurate 5th order ESDIRK method with splitting
"""
struct KenCarp5{CS, AD, F, F2, P, FDT, ST, CJ, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function KenCarp5(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI, step_limiter! = trivial_limiter!)
    KenCarp5{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve, nlsolve, precs,
        smooth_est, extrapolant, controller, step_limiter!)
end
"""
@article{kennedy2019higher,
title={Higher-order additive Runge--Kutta schemes for ordinary differential equations},
author={Kennedy, Christopher A and Carpenter, Mark H},
journal={Applied Numerical Mathematics},
volume={136},
pages={183--205},
year={2019},
publisher={Elsevier}
}

KenCarp58: SDIRK Method
An A-L stable stiffly-accurate 5th order eight-stage ESDIRK method with splitting
"""
struct KenCarp58{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
end
function KenCarp58(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :PI)
    KenCarp58{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, smooth_est, extrapolant,
        controller)
end

# `smooth_est` is not necessary, as the embedded method is also L-stable
struct ESDIRK54I8L2SA{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function ESDIRK54I8L2SA(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    ESDIRK54I8L2SA{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller)
end

"""
@article{Kennedy2019DiagonallyIR,
title={Diagonally implicit Runge–Kutta methods for stiff ODEs},
author={Christopher A. Kennedy and Mark H. Carpenter},
journal={Applied Numerical Mathematics},
year={2019},
volume={146},
pages={221-244}
}
"""
struct ESDIRK436L2SA2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function ESDIRK436L2SA2(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    ESDIRK436L2SA2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller)
end

"""
@article{Kennedy2019DiagonallyIR,
title={Diagonally implicit Runge–Kutta methods for stiff ODEs},
author={Christopher A. Kennedy and Mark H. Carpenter},
journal={Applied Numerical Mathematics},
year={2019},
volume={146},
pages={221-244}
}
"""
struct ESDIRK437L2SA{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function ESDIRK437L2SA(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    ESDIRK437L2SA{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller)
end

"""
@article{Kennedy2019DiagonallyIR,
title={Diagonally implicit Runge–Kutta methods for stiff ODEs},
author={Christopher A. Kennedy and Mark H. Carpenter},
journal={Applied Numerical Mathematics},
year={2019},
volume={146},
pages={221-244}
}
"""
struct ESDIRK547L2SA2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function ESDIRK547L2SA2(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    ESDIRK547L2SA2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller)
end

"""
@article{Kennedy2019DiagonallyIR,
title={Diagonally implicit Runge–Kutta methods for stiff ODEs},
author={Christopher A. Kennedy and Mark H. Carpenter},
journal={Applied Numerical Mathematics},
year={2019},
volume={146},
pages={221-244}

Currently has STABILITY ISSUES, causing it to fail the adaptive tests.
Check issue https://github.com/SciML/OrdinaryDiffEq.jl/issues/1933 for more details.
}
"""
struct ESDIRK659L2SA{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function ESDIRK659L2SA(; chunk_size = Val{0}(), autodiff = Val{true}(),
        standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :linear, controller = :PI)
    ESDIRK659L2SA{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve, nlsolve, precs, extrapolant,
        controller)
end

################################################################################

# Rosenbrock Methods

#=
#### Rosenbrock23, Rosenbrock32, ode23s, ModifiedRosenbrockIntegrator

- Shampine L.F. and Reichelt M., (1997) The MATLAB ODE Suite, SIAM Journal of
Scientific Computing, 18 (1), pp. 1-22.

#### ROS2

- J. G. Verwer et al. (1999): A second-order Rosenbrock method applied to photochemical dispersion problems
  https://doi.org/10.1137/S1064827597326651

#### ROS3P

- Lang, J. & Verwer, ROS3P—An Accurate Third-Order Rosenbrock Solver Designed for
  Parabolic Problems J. BIT Numerical Mathematics (2001) 41: 731. doi:10.1023/A:1021900219772

#### ROS3, Rodas3, Ros4LStab, Rodas4, Rodas42

- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)

#### ROS2PR, ROS2S, ROS3PR, Scholz4_7
-Rang, Joachim (2014): The Prothero and Robinson example:
 Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
 https://doi.org/10.24355/dbbs.084-201408121139-0

#### RosShamp4

- L. F. Shampine, Implementation of Rosenbrock Methods, ACM Transactions on
  Mathematical Software (TOMS), 8: 2, 93-113. doi:10.1145/355993.355994

#### Veldd4, Velds4

- van Veldhuizen, D-stability and Kaps-Rentrop-methods, M. Computing (1984) 32: 229.
  doi:10.1007/BF02243574

#### GRK4T, GRK4A

- Kaps, P. & Rentrop, Generalized Runge-Kutta methods of order four with stepsize control
  for stiff ordinary differential equations. P. Numer. Math. (1979) 33: 55. doi:10.1007/BF01396495

#### Rodas23W, Rodas3P

- Steinebach G., Rodas23W / Rodas32P - a Rosenbrock-type method for DAEs with additional error estimate for dense output and Julia implementation,
 in progress

#### Rodas4P

- Steinebach G. Order-reduction of ROW-methods for DAEs and method of lines
  applications. Preprint-Nr. 1741, FB Mathematik, TH Darmstadt; 1995.

#### Rodas4P2
- Steinebach G. (2020) Improvement of Rosenbrock-Wanner Method RODASP.
  In: Reis T., Grundel S., Schoeps S. (eds) Progress in Differential-Algebraic Equations II.
  Differential-Algebraic Equations Forum. Springer, Cham. https://doi.org/10.1007/978-3-030-53905-4_6

#### Rodas5
- Di Marzo G. RODAS5(4) – Méthodes de Rosenbrock d’ordre 5(4) adaptées aux problemes
différentiels-algébriques. MSc mathematics thesis, Faculty of Science,
University of Geneva, Switzerland.

#### ROS34PRw
-Joachim Rang, Improved traditional Rosenbrock–Wanner methods for stiff ODEs and DAEs,
 Journal of Computational and Applied Mathematics,
 https://doi.org/10.1016/j.cam.2015.03.010

#### ROS3PRL, ROS3PRL2
-Rang, Joachim (2014): The Prothero and Robinson example:
 Convergence studies for Runge-Kutta and Rosenbrock-Wanner methods.
 https://doi.org/10.24355/dbbs.084-201408121139-0

#### Rodas5P
- Steinebach G.   Construction of Rosenbrock–Wanner method Rodas5P and numerical benchmarks within the Julia Differential Equations package.
 In: BIT Numerical Mathematics, 63(2), 2023

 #### Rodas23W, Rodas3P, Rodas5Pe, Rodas5Pr
- Steinebach G. Rosenbrock methods within OrdinaryDiffEq.jl - Overview, recent developments and applications -
 Preprint 2024
 https://github.com/hbrs-cse/RosenbrockMethods/blob/main/paper/JuliaPaper.pdf

=#

for Alg in [
    :ROS2,
    :ROS2PR,
    :ROS2S,
    :ROS3,
    :ROS3PR,
    :Scholz4_7,
    :ROS34PW1a,
    :ROS34PW1b,
    :ROS34PW2,
    :ROS34PW3,
    :ROS34PRw,
    :ROS3PRL,
    :ROS3PRL2,
    :RosShamp4,
    :Veldd4,
    :Velds4,
    :GRK4T,
    :GRK4A,
    :Ros4LStab]
    @eval begin
        struct $Alg{CS, AD, F, P, FDT, ST, CJ} <:
               OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
            linsolve::F
            precs::P
        end
        function $Alg(; chunk_size = Val{0}(), autodiff = Val{true}(),
                standardtag = Val{true}(), concrete_jac = nothing,
                diff_type = Val{:forward}, linsolve = nothing, precs = DEFAULT_PRECS)
            $Alg{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
                typeof(precs), diff_type, _unwrap_val(standardtag),
                _unwrap_val(concrete_jac)}(linsolve,
                precs)
        end
    end

    @eval TruncatedStacktraces.@truncate_stacktrace $Alg 1 2
end

# for Rosenbrock methods with step_limiter
for Alg in [
    :Rosenbrock23,
    :Rosenbrock32,
    :ROS3P,
    :Rodas3,
    :Rodas23W,
    :Rodas3P,
    :Rodas4,
    :Rodas42,
    :Rodas4P,
    :Rodas4P2,
    :Rodas5,
    :Rodas5P,
    :Rodas5Pe,
    :Rodas5Pr]
    @eval begin
        struct $Alg{CS, AD, F, P, FDT, ST, CJ, StepLimiter} <:
               OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
            linsolve::F
            precs::P
            step_limiter!::StepLimiter
        end
        function $Alg(; chunk_size = Val{0}(), autodiff = Val{true}(),
                standardtag = Val{true}(), concrete_jac = nothing,
                diff_type = Val{:forward}, linsolve = nothing,
                precs = DEFAULT_PRECS, step_limiter! = trivial_limiter!)
            $Alg{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
                typeof(precs), diff_type, _unwrap_val(standardtag),
                _unwrap_val(concrete_jac), typeof(step_limiter!)}(linsolve,
                precs, step_limiter!)
        end
    end

    @eval TruncatedStacktraces.@truncate_stacktrace $Alg 1 2
end
struct GeneralRosenbrock{CS, AD, F, ST, CJ, TabType} <:
       OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS, AD, Val{:forward}, ST, CJ}
    tableau::TabType
    factorization::F
end

function GeneralRosenbrock(; chunk_size = Val{0}(), autodiff = true,
        standardtag = Val{true}(), concrete_jac = nothing,
        factorization = lu!, tableau = ROSENBROCK_DEFAULT_TABLEAU)
    GeneralRosenbrock{
        _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(factorization),
        _unwrap_val(standardtag), _unwrap_val(concrete_jac), typeof(tableau)}(tableau,
        factorization)
end
"""
RosenbrockW6S4OS: Rosenbrock-W Method
A 4th order L-stable Rosenbrock-W method (fixed step only).
"""
struct RosenbrockW6S4OS{CS, AD, F, P, FDT, ST, CJ} <:
       OrdinaryDiffEqRosenbrockAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    precs::P
end
function RosenbrockW6S4OS(; chunk_size = Val{0}(), autodiff = true,
        standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:central},
        linsolve = nothing,
        precs = DEFAULT_PRECS)
    RosenbrockW6S4OS{_unwrap_val(chunk_size),
        _unwrap_val(autodiff), typeof(linsolve), typeof(precs), diff_type,
        _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(linsolve,
        precs)
end

######################################

for Alg in [:LawsonEuler, :NorsettEuler, :ETDRK2, :ETDRK3, :ETDRK4, :HochOst4]
    """
    Hochbruck, Marlis, and Alexander Ostermann. “Exponential Integrators.” Acta
      Numerica 19 (2010): 209–286. doi:10.1017/S0962492910000048.
    """
    @eval struct $Alg{CS, AD, FDT, ST, CJ} <:
                 OrdinaryDiffEqExponentialAlgorithm{CS, AD, FDT, ST, CJ}
        krylov::Bool
        m::Int
        iop::Int
    end
    @eval function $Alg(; krylov = false, m = 30, iop = 0, autodiff = true,
            standardtag = Val{true}(), concrete_jac = nothing,
            chunk_size = Val{0}(),
            diff_type = Val{:forward})
        $Alg{_unwrap_val(chunk_size), _unwrap_val(autodiff),
            diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(krylov,
            m,
            iop)
    end
end
const ETD1 = NorsettEuler # alias
for Alg in [:Exprb32, :Exprb43]
    @eval struct $Alg{CS, AD, FDT, ST, CJ} <:
                 OrdinaryDiffEqAdaptiveExponentialAlgorithm{CS, AD, FDT, ST, CJ}
        m::Int
        iop::Int
    end
    @eval function $Alg(; m = 30, iop = 0, autodiff = true, standardtag = Val{true}(),
            concrete_jac = nothing, chunk_size = Val{0}(),
            diff_type = Val{:forward})
        $Alg{_unwrap_val(chunk_size), _unwrap_val(autodiff),
            diff_type, _unwrap_val(standardtag),
            _unwrap_val(concrete_jac)}(m,
            iop)
    end
end
for Alg in [:Exp4, :EPIRK4s3A, :EPIRK4s3B, :EPIRK5s3, :EXPRB53s3, :EPIRK5P1, :EPIRK5P2]
    @eval struct $Alg{CS, AD, FDT, ST, CJ} <:
                 OrdinaryDiffEqExponentialAlgorithm{CS, AD, FDT, ST, CJ}
        adaptive_krylov::Bool
        m::Int
        iop::Int
    end
    @eval function $Alg(; adaptive_krylov = true, m = 30, iop = 0, autodiff = true,
            standardtag = Val{true}(), concrete_jac = nothing,
            chunk_size = Val{0}(), diff_type = Val{:forward})
        $Alg{_unwrap_val(chunk_size), _unwrap_val(autodiff), diff_type,
            _unwrap_val(standardtag), _unwrap_val(concrete_jac)}(adaptive_krylov,
            m,
            iop)
    end
end
struct SplitEuler <:
       OrdinaryDiffEqExponentialAlgorithm{0, false, Val{:forward}, Val{true}, nothing} end
"""
ETD2: Exponential Runge-Kutta Method
Second order Exponential Time Differencing method (in development).
"""
struct ETD2 <:
       OrdinaryDiffEqExponentialAlgorithm{0, false, Val{:forward}, Val{true}, nothing} end

#########################################

"""
E. Alberdi Celayaa, J. J. Anza Aguirrezabalab, P. Chatzipantelidisc. Implementation of
an Adaptive BDF2 Formula and Comparison with The MATLAB Ode15s. Procedia Computer Science,
29, pp 1014-1026, 2014. doi: https://doi.org/10.1016/j.procs.2014.05.091

ABDF2: Multistep Method
An adaptive order 2 L-stable fixed leading coefficient multistep BDF method.
"""
struct ABDF2{CS, AD, F, F2, P, FDT, ST, CJ, K, T, StepLimiter} <:
       OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    smooth_est::Bool
    extrapolant::Symbol
    controller::Symbol
    step_limiter!::StepLimiter
end
function ABDF2(; chunk_size = Val{0}(), autodiff = true, standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        κ = nothing, tol = nothing, linsolve = nothing, precs = DEFAULT_PRECS,
        nlsolve = NLNewton(),
        smooth_est = true, extrapolant = :linear,
        controller = :Standard, step_limiter! = trivial_limiter!)
    ABDF2{
        _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve), typeof(nlsolve),
        typeof(precs), diff_type, _unwrap_val(standardtag), _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol), typeof(step_limiter!)}(linsolve, nlsolve, precs, κ, tol,
        smooth_est, extrapolant, controller, step_limiter!)
end

#########################################

struct CompositeAlgorithm{T, F} <: OrdinaryDiffEqCompositeAlgorithm
    algs::T
    choice_function::F
end

TruncatedStacktraces.@truncate_stacktrace CompositeAlgorithm 1

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, :silence!)
    Base.Experimental.silence!(CompositeAlgorithm)
end

mutable struct AutoSwitchCache{nAlg, sAlg, tolType, T}
    count::Int
    successive_switches::Int
    nonstiffalg::nAlg
    stiffalg::sAlg
    is_stiffalg::Bool
    maxstiffstep::Int
    maxnonstiffstep::Int
    nonstifftol::tolType
    stifftol::tolType
    dtfac::T
    stiffalgfirst::Bool
    switch_max::Int
    current::Int
end

struct AutoSwitch{nAlg, sAlg, tolType, T}
    nonstiffalg::nAlg
    stiffalg::sAlg
    maxstiffstep::Int
    maxnonstiffstep::Int
    nonstifftol::tolType
    stifftol::tolType
    dtfac::T
    stiffalgfirst::Bool
    switch_max::Int
end

################################################################################
"""
MEBDF2: Multistep Method
The second order Modified Extended BDF method, which has improved stability properties over the standard BDF.
Fixed timestep only.
"""
struct MEBDF2{CS, AD, F, F2, P, FDT, ST, CJ} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
end
function MEBDF2(; chunk_size = Val{0}(), autodiff = true, standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant)
    MEBDF2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve,
        precs,
        extrapolant)
end

#################################################
"""
PDIRK44: Parallel Diagonally Implicit Runge-Kutta Method
A 2 processor 4th order diagonally non-adaptive implicit method.
"""
struct PDIRK44{CS, AD, F, F2, P, FDT, ST, CJ, TO} <:
       OrdinaryDiffEqNewtonAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    threading::TO
end
function PDIRK44(; chunk_size = Val{0}(), autodiff = true, standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant, threading = true)
    PDIRK44{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac), typeof(threading)}(linsolve, nlsolve, precs,
        extrapolant, threading)
end
### Algorithm Groups

const MultistepAlgorithms = Union{IRKN3, IRKN4,
    ABDF2,
    AB3, AB4, AB5, ABM32, ABM43, ABM54}

const SplitAlgorithms = Union{CNAB2, CNLF2, IRKC, SBDF,
    KenCarp3, KenCarp4, KenCarp47, KenCarp5, KenCarp58, CFNLIRK3}

#=
struct DBDF{CS,AD,F,F2,P,FDT,ST,CJ} <: DAEAlgorithm{CS,AD,FDT,ST,CJ}
  linsolve::F
  nlsolve::F2
  precs::P
  extrapolant::Symbol
end

DBDF(;chunk_size=Val{0}(),autodiff=Val{true}(), standardtag = Val{true}(), concrete_jac = nothing,diff_type=Val{:forward},
     linsolve=nothing,precs = DEFAULT_PRECS,nlsolve=NLNewton(),extrapolant=:linear) =
     DBDF{_unwrap_val(chunk_size),_unwrap_val(autodiff),typeof(linsolve),typeof(nlsolve),typeof(precs),diff_type,_unwrap_val(standardtag),_unwrap_val(concrete_jac)}(
     linsolve,nlsolve,precs,extrapolant)
=#

struct DImplicitEuler{CS, AD, F, F2, P, FDT, ST, CJ} <: DAEAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function DImplicitEuler(;
        chunk_size = Val{0}(), autodiff = true, standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant,
        controller = :Standard)
    DImplicitEuler{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve, precs, extrapolant, controller)
end

struct DABDF2{CS, AD, F, F2, P, FDT, ST, CJ} <: DAEAlgorithm{CS, AD, FDT, ST, CJ}
    linsolve::F
    nlsolve::F2
    precs::P
    extrapolant::Symbol
    controller::Symbol
end
function DABDF2(; chunk_size = Val{0}(), autodiff = Val{true}(), standardtag = Val{true}(),
        concrete_jac = nothing, diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(),
        extrapolant = :constant,
        controller = :Standard)
    DABDF2{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac)}(linsolve,
        nlsolve, precs, extrapolant, controller)
end

struct DFBDF{MO, CS, AD, F, F2, P, FDT, ST, CJ, K, T} <: DAEAlgorithm{CS, AD, FDT, ST, CJ}
    max_order::Val{MO}
    linsolve::F
    nlsolve::F2
    precs::P
    κ::K
    tol::T
    extrapolant::Symbol
    controller::Symbol
end
function DFBDF(; max_order::Val{MO} = Val{5}(), chunk_size = Val{0}(),
        autodiff = Val{true}(), standardtag = Val{true}(), concrete_jac = nothing,
        diff_type = Val{:forward},
        linsolve = nothing, precs = DEFAULT_PRECS, nlsolve = NLNewton(), κ = nothing,
        tol = nothing,
        extrapolant = :linear, controller = :Standard) where {MO}
    DFBDF{MO, _unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
        typeof(nlsolve), typeof(precs), diff_type, _unwrap_val(standardtag),
        _unwrap_val(concrete_jac),
        typeof(κ), typeof(tol)}(max_order, linsolve, nlsolve, precs, κ, tol, extrapolant,
        controller)
end

TruncatedStacktraces.@truncate_stacktrace DFBDF
