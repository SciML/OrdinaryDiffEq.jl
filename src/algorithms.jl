abstract type OrdinaryDiffEqAlgorithm <: DiffEqBase.AbstractODEAlgorithm end
abstract type OrdinaryDiffEqAdaptiveAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqCompositeAlgorithm <: OrdinaryDiffEqAlgorithm end

abstract type OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD} <: OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD} <: OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD} end
abstract type OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD} <: OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD} end

abstract type OrdinaryDiffEqImplicitAlgorithm{CS,AD} <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqNewtonAlgorithm{CS,AD} <:  OrdinaryDiffEqImplicitAlgorithm{CS,AD} end
abstract type OrdinaryDiffEqRosenbrockAlgorithm{CS,AD} <:  OrdinaryDiffEqImplicitAlgorithm{CS,AD} end
const NewtonAlgorithm = Union{OrdinaryDiffEqNewtonAlgorithm,OrdinaryDiffEqNewtonAdaptiveAlgorithm}
const RosenbrockAlgorithm = Union{OrdinaryDiffEqRosenbrockAlgorithm,OrdinaryDiffEqRosenbrockAdaptiveAlgorithm}

abstract type OrdinaryDiffEqExponentialAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqAdaptiveExponentialAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqLinearExponentialAlgorithm <: OrdinaryDiffEqExponentialAlgorithm end
const ExponentialAlgorithm = Union{OrdinaryDiffEqExponentialAlgorithm,OrdinaryDiffEqAdaptiveExponentialAlgorithm}

abstract type OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
abstract type OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS,AD} <: OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD} end

struct FunctionMap{scale_by_time} <: OrdinaryDiffEqAlgorithm end
FunctionMap(;scale_by_time=false) = FunctionMap{scale_by_time}()

function DiffEqBase.remake(thing::OrdinaryDiffEqAlgorithm, kwargs...)
  T = DiffEqBase.remaker_of(thing)
  T(; DiffEqBase.struct_as_namedtuple(thing)...,kwargs...)
end

function DiffEqBase.remake(thing::OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS,AD}, kwargs...) where {CS, AD}
  T = DiffEqBase.remaker_of(thing)
  T(; chunk_size=CS,autodiff=AD,DiffEqBase.struct_as_namedtuple(thing)...,kwargs...)
end

function DiffEqBase.remake(thing::OrdinaryDiffEqImplicitAlgorithm{CS,AD}, kwargs...) where {CS, AD}
  T = DiffEqBase.remaker_of(thing)
  T(; chunk_size=CS,autodiff=AD,DiffEqBase.struct_as_namedtuple(thing)...,kwargs...)
end
###############################################################################

# RK methods

struct ExplicitRK{TabType} <: OrdinaryDiffEqAdaptiveAlgorithm
  tableau::TabType
end
ExplicitRK(;tableau=ODE_DEFAULT_TABLEAU) = ExplicitRK(tableau)

@inline trivial_limiter!(u, f, p, t) = nothing
"""
Euler - The canonical forward Euler method. Fixed timestep only.
"""
struct Euler <: OrdinaryDiffEqAlgorithm end
"""
KuttaPRK2p5: Parallel Explicit Runge-Kutta Method
  A 5 parallel, 2 processor explicit Runge-Kutta method of 5th order.

  These methods utilize multithreading on the f calls to parallelize the problem.
  This requires that simultaneous calls to f are thread-safe.
"""
struct KuttaPRK2p5{TO} <: OrdinaryDiffEqAlgorithm
  threading::TO
end
KuttaPRK2p5(;threading=true) = KuttaPRK2p5(threading)
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
AitkenNeville(;max_order=10,min_order=1,init_order=5,threading=true) = AitkenNeville(max_order,min_order,init_order,threading)
"""
ImplicitEulerExtrapolation: Parallelized Implicit Extrapolation Method
   Extrapolation of implicit Euler method with Romberg sequence.
   Similar to Hairer's SEULEX.
"""
struct ImplicitEulerExtrapolation{CS,AD,F,FDT,TO} <: OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS,AD}
  linsolve::F
  n_max::Int
  n_min::Int
  n_init::Int
  threading::TO
  diff_type::FDT
  sequence::Symbol # Name of the subdividing sequence
end

function ImplicitEulerExtrapolation(;chunk_size=0,autodiff=true,
    diff_type=Val{:forward},linsolve=DEFAULT_LINSOLVE,
    max_order=12,min_order=3,init_order=5,threading=true,sequence = :bulirsch)

    n_min = max(3,min_order)
    n_init = max(n_min + 1,init_order)
    n_max = max(n_init + 1, max_order)
    if threading
      @warn "Threading in `ImplicitEulerExtrapolation` is currently disabled. Thus `threading` has been changed from `true` to `false`."
      threading = false
    end

    # Warn user if sequence has been changed:
    if sequence != :harmonic && sequence != :romberg && sequence != :bulirsch
      @warn "The `sequence` given to the `ImplicitEulerExtrapolation` algorithm
          is not valid: it must match `:harmonic`, `:romberg` or `:bulirsch`.
          Thus it has been changed
        :$(sequence) --> :bulirsch"
      sequence = :bulirsch
    end
    ImplicitEulerExtrapolation{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),typeof(threading)}(
      linsolve,n_max,n_min,n_init,threading,diff_type,sequence)
end
"""
ExtrapolationMidpointDeuflhard: Parallelized Explicit Extrapolation Method
   Midpoint extrapolation using Barycentric coordinates
"""
struct ExtrapolationMidpointDeuflhard{TO} <: OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm
  n_min::Int # Minimal extrapolation order
  n_init::Int # Initial extrapolation order
  n_max::Int # Maximal extrapolation order
  sequence::Symbol # Name of the subdividing sequence
  threading::TO
  sequence_factor::Int # An even factor by which sequence is scaled for midpoint extrapolation
end
function ExtrapolationMidpointDeuflhard(;min_order=1,init_order=5, max_order=10, sequence = :harmonic, threading = true, sequence_factor = 2)
  # Enforce 1 <=  min_order <= init_order <= max_order:
  n_min = max(1,min_order)
  n_init = max(n_min,init_order)
  n_max = max(n_init,max_order)

  # Warn user if orders have been changed
  if (min_order, init_order, max_order) != (n_min,n_init,n_max)
    @warn "The range of extrapolation orders and/or the initial order given to the
      `ExtrapolationMidpointDeuflhard` algorithm are not valid and have been changed:
      Minimal order: " * lpad(min_order,2," ") * " --> "  * lpad(n_min,2," ") * "
      Maximal order: " * lpad(max_order,2," ") * " --> "  * lpad(n_max,2," ") * "
      Initial order: " * lpad(init_order,2," ") * " --> "  * lpad(n_init,2," ")
  end

  # Warn user if sequence_factor is not even
  if sequence_factor%2 != 0
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
  ExtrapolationMidpointDeuflhard(n_min,n_init,n_max,sequence,threading,sequence_factor)
end
"""
ImplicitDeuflhardExtrapolation: Parallelized Implicit Extrapolation Method
   Midpoint extrapolation using Barycentric coordinates
"""
struct ImplicitDeuflhardExtrapolation{CS,AD,F,FDT,TO} <: OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS,AD}
  linsolve::F
  n_min::Int # Minimal extrapolation order
  n_init::Int # Initial extrapolation order
  n_max::Int # Maximal extrapolation order
  sequence::Symbol # Name of the subdividing sequence
  diff_type::FDT
  threading::TO
end
function ImplicitDeuflhardExtrapolation(;chunk_size=0,autodiff=true,
  linsolve=DEFAULT_LINSOLVE,diff_type=Val{:forward},
  min_order=1,init_order=5,max_order=10,sequence = :harmonic,threading=false)
  # Enforce 1 <=  min_order <= init_order <= max_order:
  n_min = max(1,min_order)
  n_init = max(n_min,init_order)
  n_max = max(n_init,max_order)

  # Warn user if orders have been changed
  if (min_order, init_order, max_order) != (n_min,n_init,n_max)
    @warn "The range of extrapolation orders and/or the initial order given to the
      `ImplicitDeuflhardExtrapolation` algorithm are not valid and have been changed:
      Minimal order: " * lpad(min_order,2," ") * " --> "  * lpad(n_min,2," ") * "
      Maximal order: " * lpad(max_order,2," ") * " --> "  * lpad(n_max,2," ") * "
      Initial order: " * lpad(init_order,2," ") * " --> "  * lpad(n_init,2," ")
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
  ImplicitDeuflhardExtrapolation{chunk_size, autodiff,
      typeof(linsolve), typeof(diff_type), typeof(threading)}(linsolve,n_min,n_init,n_max,sequence,diff_type,threading)
end
"""
ExtrapolationMidpointHairerWanner: Parallelized Explicit Extrapolation Method
  Midpoint extrapolation using Barycentric coordinates, following Hairer's ODEX in the adaptivity behavior.
"""
struct ExtrapolationMidpointHairerWanner{TO} <: OrdinaryDiffEqExtrapolationVarOrderVarStepAlgorithm
  n_min::Int # Minimal extrapolation order
  n_init::Int # Initial extrapolation order
  n_max::Int # Maximal extrapolation order
  sequence::Symbol # Name of the subdividing sequence
  threading::TO
  sequence_factor::Int # An even factor by which sequence is scaled for midpoint extrapolation
end
function ExtrapolationMidpointHairerWanner(;min_order=2,init_order=5, max_order=10, sequence = :harmonic, threading = true, sequence_factor = 2)
  # Enforce 2 <=  min_order
  # and min_order + 1 <= init_order <= max_order - 1:
  n_min = max(2, min_order)
  n_init = max(n_min + 1, init_order)
  n_max = max(n_init + 1, max_order)

  # Warn user if orders have been changed
  if (min_order, init_order, max_order) != (n_min,n_init,n_max)
    @warn "The range of extrapolation orders and/or the initial order given to the
      `ExtrapolationMidpointHairerWanner` algorithm are not valid and have been changed:
      Minimal order: " * lpad(min_order,2," ") * " --> "  * lpad(n_min,2," ") * "
      Maximal order: " * lpad(max_order,2," ") * " --> "  * lpad(n_max,2," ") * "
      Initial order: " * lpad(init_order,2," ") * " --> "  * lpad(n_init,2," ")
  end

  # Warn user if sequence_factor is not even
  if sequence_factor%2 != 0
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
  ExtrapolationMidpointHairerWanner(n_min,n_init,n_max,sequence,threading,sequence_factor)
end
"""
ImplicitHairerWannerExtrapolation: Parallelized Implicit Extrapolation Method
  Midpoint extrapolation using Barycentric coordinates, following Hairer's SODEX in the adaptivity behavior.
"""
struct ImplicitHairerWannerExtrapolation{CS,AD,F,FDT,TO} <: OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS,AD}
  linsolve::F
  n_min::Int # Minimal extrapolation order
  n_init::Int # Initial extrapolation order
  n_max::Int # Maximal extrapolation order
  sequence::Symbol # Name of the subdividing sequence
  diff_type::FDT
  threading::TO
end
function ImplicitHairerWannerExtrapolation(;chunk_size=0,autodiff=true,
  linsolve=DEFAULT_LINSOLVE,diff_type=Val{:forward},
  min_order=2,init_order=5,max_order=10,sequence = :harmonic,threading=false)
  # Enforce 2 <=  min_order
  # and min_order + 1 <= init_order <= max_order - 1:
  n_min = max(2, min_order)
  n_init = max(n_min + 1, init_order)
  n_max = max(n_init + 1, max_order)

  # Warn user if orders have been changed
  if (min_order, init_order, max_order) != (n_min,n_init,n_max)
    @warn "The range of extrapolation orders and/or the initial order given to the
      `ImplicitHairerWannerExtrapolation` algorithm are not valid and have been changed:
      Minimal order: " * lpad(min_order,2," ") * " --> "  * lpad(n_min,2," ") * "
      Maximal order: " * lpad(max_order,2," ") * " --> "  * lpad(n_max,2," ") * "
      Initial order: " * lpad(init_order,2," ") * " --> "  * lpad(n_init,2," ")
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
  ImplicitHairerWannerExtrapolation{chunk_size, autodiff,
      typeof(linsolve), typeof(diff_type), typeof(threading)}(linsolve,n_min,n_init,n_max,
      sequence,diff_type,threading)
end

struct ImplicitEulerBarycentricExtrapolation{CS,AD,F,FDT,TO} <: OrdinaryDiffEqImplicitExtrapolationAlgorithm{CS,AD}
  linsolve::F
  n_min::Int # Minimal extrapolation order
  n_init::Int # Initial extrapolation order
  n_max::Int # Maximal extrapolation order
  sequence::Symbol # Name of the subdividing sequence
  diff_type::FDT
  threading::TO
  sequence_factor::Int
end
function ImplicitEulerBarycentricExtrapolation(;chunk_size=0,autodiff=true,
  linsolve=DEFAULT_LINSOLVE,diff_type=Val{:forward},
  min_order=3,init_order=5,max_order=12,sequence = :harmonic,threading=false,sequence_factor = 2)
  # Enforce 2 <=  min_order
  # and min_order + 1 <= init_order <= max_order - 1:
  n_min = max(3, min_order)
  n_init = max(n_min + 1, init_order)
  n_max = max(n_init + 1, max_order)

  # Warn user if orders have been changed
  if (min_order, init_order, max_order) != (n_min,n_init,n_max)
    @warn "The range of extrapolation orders and/or the initial order given to the
      `ImplicitEulerBarycentricExtrapolation` algorithm are not valid and have been changed:
      Minimal order: " * lpad(min_order,2," ") * " --> "  * lpad(n_min,2," ") * "
      Maximal order: " * lpad(max_order,2," ") * " --> "  * lpad(n_max,2," ") * "
      Initial order: " * lpad(init_order,2," ") * " --> "  * lpad(n_init,2," ")
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
  ImplicitEulerBarycentricExtrapolation{chunk_size, autodiff,
      typeof(linsolve), typeof(diff_type), typeof(threading)}(linsolve,n_min,n_init,n_max,
      sequence,diff_type,threading,sequence_factor)
end

"""
Julien Berland, Christophe Bogey, Christophe Bailly. Low-Dissipation and Low-Dispersion
Fourth-Order Runge-Kutta Algorithm. Computers & Fluids, 35(10), pp 1459-1463, 2006.
doi: https://doi.org/10.1016/j.compfluid.2005.04.003
"""
"""
RK46NL: 6-stage, fourth order low-stage, low-dissipation, low-dispersion scheme.
        Fixed timestep only.
"""
struct RK46NL <: OrdinaryDiffEqAlgorithm end
"""
Heun: 2nd order Heun's method
"""
struct Heun <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Ralston <: OrdinaryDiffEqAdaptiveAlgorithm end
struct Midpoint <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
@article{shampine2005solving,
  title={Solving ODEs and DDEs with residual control},
  author={Shampine, LF},
  journal={Applied Numerical Mathematics},
  volume={52},
  number={1},
  pages={113--127},
  year={2005},
  publisher={Elsevier}
}
"""
struct RK4 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct RKM <: OrdinaryDiffEqAlgorithm end
struct Anas5{T} <: OrdinaryDiffEqAlgorithm
  w::T
end
Anas5(; w=1) = Anas5(w)

"""
Matteo Bernardini, Sergio Pirozzoli. A General Strategy for the Optimization of
Runge-Kutta Schemes for Wave Propagation Phenomena. Journal of Computational Physics,
228(11), pp 4182-4199, 2009. doi: https://doi.org/10.1016/j.jcp.2009.02.032
"""
struct ORK256{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  williamson_condition::Bool
  ORK256(stage_limiter! =trivial_limiter!, step_limiter! =trivial_limiter!; williamson_condition=true) = new{typeof(stage_limiter!), typeof(step_limiter!)}(stage_limiter!, step_limiter!, williamson_condition)
end

"""
@article{carpenter1994fourth,
  title={Fourth-order 2N-storage Runge-Kutta schemes},
  author={Carpenter, Mark H and Kennedy, Christopher A},
  year={1994}
}
"""
struct CarpenterKennedy2N54{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  williamson_condition::Bool
  CarpenterKennedy2N54(stage_limiter! =trivial_limiter!, step_limiter! =trivial_limiter!; williamson_condition=true) = new{typeof(stage_limiter!), typeof(step_limiter!)}(stage_limiter!, step_limiter!, williamson_condition)
end


"""
D. Stanescu, W. G. Habashi. 2N-Storage Low Dissipation and Dispersion Runge-Kutta Schemes for
Computational Acoustics. Journal of Computational Physics, 143(2), pp 674-681, 1998. doi:
https://doi.org/10.1006/jcph.1998.5986
"""
struct SHLDDRK64{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  williamson_condition::Bool
  SHLDDRK64(stage_limiter! =trivial_limiter!, step_limiter! =trivial_limiter!; williamson_condition=true) = new{typeof(stage_limiter!), typeof(step_limiter!)}(stage_limiter!, step_limiter!, williamson_condition)
end
struct SHLDDRK52 <: OrdinaryDiffEqAlgorithm end
struct SHLDDRK_2N <: OrdinaryDiffEqAlgorithm end
"""
Deprecated SHLDDRK64 scheme from 'D. Stanescu, W. G. Habashi. 2N-Storage Low Dissipation and Dispersion Runge-Kutta Schemes for
Computational Acoustics'
"""
struct HSLDDRK64{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  williamson_condition::Bool
  function HSLDDRK64(stage_limiter! =trivial_limiter!, step_limiter! =trivial_limiter!; williamson_condition=true)
    Base.depwarn("HSLDDRK64 is deprecated, use SHLDDRK64 instead.", :HSLDDRK64)
    SHLDDRK64(stage_limiter!, step_limiter!; williamson_condition=williamson_condition)
  end
end

"""
T. Toulorge, W. Desmet. Optimal Runge–Kutta Schemes for Discontinuous Galerkin Space
Discretizations Applied to Wave Propagation Problems. Journal of Computational Physics, 231(4),
pp 2067-2091, 2012. doi: https://doi.org/10.1016/j.jcp.2011.11.024
"""
struct DGLDDRK73_C{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  williamson_condition::Bool
  DGLDDRK73_C(stage_limiter! =trivial_limiter!, step_limiter! =trivial_limiter!; williamson_condition=true) = new{typeof(stage_limiter!), typeof(step_limiter!)}(stage_limiter!, step_limiter!, williamson_condition)
end

"""
T. Toulorge, W. Desmet. Optimal Runge–Kutta Schemes for Discontinuous Galerkin Space
Discretizations Applied to Wave Propagation Problems. Journal of Computational Physics, 231(4),
pp 2067-2091, 2012. doi: https://doi.org/10.1016/j.jcp.2011.11.024
"""
struct DGLDDRK84_C{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  williamson_condition::Bool
  DGLDDRK84_C(stage_limiter! =trivial_limiter!, step_limiter! =trivial_limiter!; williamson_condition=true) = new{typeof(stage_limiter!), typeof(step_limiter!)}(stage_limiter!, step_limiter!, williamson_condition)
end

"""
T. Toulorge, W. Desmet. Optimal Runge–Kutta Schemes for Discontinuous Galerkin Space
Discretizations Applied to Wave Propagation Problems. Journal of Computational Physics, 231(4),
pp 2067-2091, 2012. doi: https://doi.org/10.1016/j.jcp.2011.11.024
"""
struct DGLDDRK84_F{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  williamson_condition::Bool
  DGLDDRK84_F(stage_limiter! =trivial_limiter!, step_limiter! =trivial_limiter!; williamson_condition=true) = new{typeof(stage_limiter!), typeof(step_limiter!)}(stage_limiter!, step_limiter!, williamson_condition)
end

"""
Jens Niegemann, Richard Diehl, Kurt Busch. Efficient Low-Storage Runge–Kutta Schemes with
Optimized Stability Regions. Journal of Computational Physics, 231, pp 364-372, 2012.
doi: https://doi.org/10.1016/j.jcp.2011.09.003
"""
struct NDBLSRK124{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  williamson_condition::Bool
  NDBLSRK124(stage_limiter! =trivial_limiter!, step_limiter! =trivial_limiter!; williamson_condition=true) = new{typeof(stage_limiter!), typeof(step_limiter!)}(stage_limiter!, step_limiter!, williamson_condition)
end

"""
Jens Niegemann, Richard Diehl, Kurt Busch. Efficient Low-Storage Runge–Kutta Schemes with
Optimized Stability Regions. Journal of Computational Physics, 231, pp 364-372, 2012.
doi: https://doi.org/10.1016/j.jcp.2011.09.003
"""
struct NDBLSRK134{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  williamson_condition::Bool
  NDBLSRK134(stage_limiter! =trivial_limiter!, step_limiter! =trivial_limiter!; williamson_condition=true) = new{typeof(stage_limiter!), typeof(step_limiter!)}(stage_limiter!, step_limiter!, williamson_condition)
end

"""
Jens Niegemann, Richard Diehl, Kurt Busch. Efficient Low-Storage Runge–Kutta Schemes with
Optimized Stability Regions. Journal of Computational Physics, 231, pp 364-372, 2012.
doi: https://doi.org/10.1016/j.jcp.2011.09.003
"""
struct NDBLSRK144{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
  williamson_condition::Bool
  NDBLSRK144(stage_limiter! =trivial_limiter!, step_limiter! =trivial_limiter!; williamson_condition=true) = new{typeof(stage_limiter!), typeof(step_limiter!)}(stage_limiter!, step_limiter!, williamson_condition)
end

"""
M. Calvo, J. M. Franco, L. Randez. A New Minimum Storage Runge–Kutta Scheme
for Computational Acoustics. Journal of Computational Physics, 201, pp 1-12, 2004.
doi: https://doi.org/10.1016/j.jcp.2004.05.012
"""
struct CFRLDDRK64 <: OrdinaryDiffEqAlgorithm end

"""
Kostas Tselios, T. E. Simos. Optimized Runge–Kutta Methods with Minimal Dispersion and Dissipation
for Problems arising from Computational Ccoustics. Physics Letters A, 393(1-2), pp 38-47, 2007.
doi: https://doi.org/10.1016/j.physleta.2006.10.072
"""
struct TSLDDRK74 <: OrdinaryDiffEqAlgorithm end
struct CKLLSRK43_2 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3C <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK95_4S <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK95_4C <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK95_4M <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3C_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3M_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3N_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK85_4C_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK85_4M_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK85_4P_3R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3N_4R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK54_3M_4R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK65_4M_4R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK85_4FM_4R <: OrdinaryDiffEqAdaptiveAlgorithm end
struct CKLLSRK75_4M_5R <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
"Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems."
SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
doi: https://doi.org/10.1137/120885899
"""
struct ParsaniKetchesonDeconinck3S32 <: OrdinaryDiffEqAlgorithm end

"""
Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
"Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems."
SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
doi: https://doi.org/10.1137/120885899
"""
struct ParsaniKetchesonDeconinck3S82 <: OrdinaryDiffEqAlgorithm end

"""
Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
"Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems."
SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
doi: https://doi.org/10.1137/120885899
"""
struct ParsaniKetchesonDeconinck3S53 <: OrdinaryDiffEqAlgorithm end

"""
Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
"Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems."
SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
doi: https://doi.org/10.1137/120885899
"""
struct ParsaniKetchesonDeconinck3S173 <: OrdinaryDiffEqAlgorithm end

"""
Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
"Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems."
SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
doi: https://doi.org/10.1137/120885899
"""
struct ParsaniKetchesonDeconinck3S94 <: OrdinaryDiffEqAlgorithm end

"""
Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
"Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems."
SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
doi: https://doi.org/10.1137/120885899
"""
struct ParsaniKetchesonDeconinck3S184 <: OrdinaryDiffEqAlgorithm end

"""
Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
"Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems."
SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
doi: https://doi.org/10.1137/120885899
"""
struct ParsaniKetchesonDeconinck3S105 <: OrdinaryDiffEqAlgorithm end

"""
Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
"Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems."
SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
doi: https://doi.org/10.1137/120885899
"""
struct ParsaniKetchesonDeconinck3S205 <: OrdinaryDiffEqAlgorithm end

"""
    RDPK3Sp35()

A third-order, five-stage explicit Runge-Kutta method with embedded error estimator
designed for spectral element discretizations of compressible fluid mechanics.

## References
- Ranocha, Dalcin, Parsani, Ketcheson (2021)
  Optimized Runge-Kutta Methods with Automatic Step Size Control for
  Compressible Computational Fluid Dynamics
  [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)
"""
struct RDPK3Sp35 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
    RDPK3SpFSAL35()

A third-order, five-stage explicit Runge-Kutta method with embedded error estimator
using the FSAL property designed for spectral element discretizations of
compressible fluid mechanics.

## References
- Ranocha, Dalcin, Parsani, Ketcheson (2021)
  Optimized Runge-Kutta Methods with Automatic Step Size Control for
  Compressible Computational Fluid Dynamics
  [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)
"""
struct RDPK3SpFSAL35 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
    RDPK3Sp49()

A fourth-order, nine-stage explicit Runge-Kutta method with embedded error estimator
designed for spectral element discretizations of compressible fluid mechanics.

## References
- Ranocha, Dalcin, Parsani, Ketcheson (2021)
  Optimized Runge-Kutta Methods with Automatic Step Size Control for
  Compressible Computational Fluid Dynamics
  [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)
"""
struct RDPK3Sp49 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
    RDPK3SpFSAL49()

A fourth-order, nine-stage explicit Runge-Kutta method with embedded error estimator
using the FSAL property designed for spectral element discretizations of
compressible fluid mechanics.

## References
- Ranocha, Dalcin, Parsani, Ketcheson (2021)
  Optimized Runge-Kutta Methods with Automatic Step Size Control for
  Compressible Computational Fluid Dynamics
  [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)
"""
struct RDPK3SpFSAL49 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
    RDPK3Sp510()

A fifth-order, ten-stage explicit Runge-Kutta method with embedded error estimator
designed for spectral element discretizations of compressible fluid mechanics.

## References
- Ranocha, Dalcin, Parsani, Ketcheson (2021)
  Optimized Runge-Kutta Methods with Automatic Step Size Control for
  Compressible Computational Fluid Dynamics
  [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)
"""
struct RDPK3Sp510 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
    RDPK3SpFSAL510()

A fifth-order, ten-stage explicit Runge-Kutta method with embedded error estimator
using the FSAL property designed for spectral element discretizations of
compressible fluid mechanics.

## References
- Ranocha, Dalcin, Parsani, Ketcheson (2021)
  Optimized Runge-Kutta Methods with Automatic Step Size Control for
  Compressible Computational Fluid Dynamics
  [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)
"""
struct RDPK3SpFSAL510 <: OrdinaryDiffEqAdaptiveAlgorithm end

struct KYK2014DGSSPRK_3S2 <: OrdinaryDiffEqAlgorithm end

"""
Tsitouras, Ch. "Explicit Runge–Kutta methods for starting integration of
Lane–Emden problem." Applied Mathematics and Computation 354 (2019): 353-364.
doi: https://doi.org/10.1016/j.amc.2019.02.047
"""
struct RKO65 <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

"""
Shu, Chi-Wang, and Stanley Osher. "Efficient implementation of essentially
non-oscillatory shock-capturing schemes." Journal of Computational Physics
77.2 (1988): 439-471.
"""
struct SSPRK22{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK22(stage_limiter! = trivial_limiter!) = SSPRK22(stage_limiter!, trivial_limiter!)

"""
Shu, Chi-Wang, and Stanley Osher. "Efficient implementation of essentially
non-oscillatory shock-capturing schemes." Journal of Computational Physics
77.2 (1988): 439-471.
"""
struct SSPRK33{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK33(stage_limiter! = trivial_limiter!) = SSPRK33(stage_limiter!, trivial_limiter!)

"""
Ruuth, Steven. "Global optimization of explicit strong-stability-preserving
Runge-Kutta methods." Mathematics of Computation 75.253 (2006): 183-207.
"""
struct SSPRK53{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
struct KYKSSPRK42 <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end
SSPRK53(stage_limiter! = trivial_limiter!) = SSPRK53(stage_limiter!, trivial_limiter!)

"""
Higueras and T. Roldán. "New third order low-storage SSP explicit Runge–Kutta methods". arXiv:1809.04807v1.
"""
struct SSPRK53_2N1{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK53_2N1(stage_limiter! = trivial_limiter!) = SSPRK53_2N1(stage_limiter!, trivial_limiter!)

"""
Higueras and T. Roldán. "New third order low-storage SSP explicit Runge–Kutta methods". arXiv:1809.04807v1.
"""
struct SSPRK53_2N2{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK53_2N2(stage_limiter! = trivial_limiter!) = SSPRK53_2N2(stage_limiter!, trivial_limiter!)

"""
Higueras and T. Roldán. "New third order low-storage SSP explicit Runge–Kutta methods". arXiv:1809.04807v1.
"""
struct SSPRK53_H{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK53_H(stage_limiter! = trivial_limiter!) = SSPRK53_H(stage_limiter!, trivial_limiter!)

"""
Ruuth, Steven. "Global optimization of explicit strong-stability-preserving
Runge-Kutta methods." Mathematics of Computation 75.253 (2006): 183-207.
"""
struct SSPRK63{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK63(stage_limiter! = trivial_limiter!) = SSPRK63(stage_limiter!, trivial_limiter!)

"""
Ruuth, Steven. "Global optimization of explicit strong-stability-preserving
Runge-Kutta methods." Mathematics of Computation 75.253 (2006): 183-207.
"""
struct SSPRK73{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK73(stage_limiter! = trivial_limiter!) = SSPRK73(stage_limiter!, trivial_limiter!)

"""
Ruuth, Steven. "Global optimization of explicit strong-stability-preserving
Runge-Kutta methods." Mathematics of Computation 75.253 (2006): 183-207.
"""
struct SSPRK83{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK83(stage_limiter! = trivial_limiter!) = SSPRK83(stage_limiter!, trivial_limiter!)

"""
Optimal third-order explicit SSP method with four stages discovered by
- J. F. B. M. Kraaijevanger.
  "Contractivity of Runge-Kutta methods."
  In: BIT Numerical Mathematics 31.3 (1991), pp. 482–528.
  [DOI: 10.1007/BF01933264](https://doi.org/10.1007/BF01933264).

Embedded method constructed by
- Sidafa Conde, Imre Fekete, John N. Shadid.
  "Embedded error estimation and adaptive step-size control for
  optimal explicit strong stability preserving Runge–Kutta methods."
  [arXiv: 1806.08693](https://arXiv.org/abs/1806.08693)

Efficient implementation (and optimized controller) described in
- Ranocha, Dalcin, Parsani, Ketcheson (2021)
  Optimized Runge-Kutta Methods with Automatic Step Size Control for
  Compressible Computational Fluid Dynamics
  [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)
"""
struct SSPRK43{StageLimiter,StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK43(stage_limiter! = trivial_limiter!) = SSPRK43(stage_limiter!, trivial_limiter!)

"""
Gottlieb, Sigal, David I. Ketcheson, and Chi-Wang Shu.
Strong stability preserving Runge-Kutta and multistep time discretizations.
World Scientific, 2011.
Example 6.1.

Consider using `SSPRK43` instead, which uses the same main method and an improved embedded method.
"""
struct SSPRK432{StageLimiter,StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK432(stage_limiter! = trivial_limiter!) = SSPRK432(stage_limiter!, trivial_limiter!)

struct SSPRKMSVS43{StageLimiter,StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRKMSVS43(stage_limiter! = trivial_limiter!) = SSPRKMSVS43(stage_limiter!, trivial_limiter!)
struct SSPRKMSVS32{StageLimiter,StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRKMSVS32(stage_limiter! = trivial_limiter!) = SSPRKMSVS32(stage_limiter!, trivial_limiter!)

"""
Gottlieb, Sigal, David I. Ketcheson, and Chi-Wang Shu. Strong stability
preserving Runge-Kutta and multistep time discretizations. World Scientific,
2011.
"""
struct SSPRK932{StageLimiter,StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK932(stage_limiter! = trivial_limiter!) = SSPRK932(stage_limiter!, trivial_limiter!)

"""
Ruuth, Steven. "Global optimization of explicit strong-stability-preserving
Runge-Kutta methods." Mathematics of Computation 75.253 (2006): 183-207.
"""
struct SSPRK54{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK54(stage_limiter! = trivial_limiter!) = SSPRK54(stage_limiter!, trivial_limiter!)

"""
Ketcheson, David I. "Highly efficient strong stability-preserving Runge–Kutta
methods with low-storage implementations." SIAM Journal on Scientific
Computing 30.4 (2008): 2113-2136.
"""
struct SSPRK104{StageLimiter,StepLimiter} <: OrdinaryDiffEqAlgorithm
  stage_limiter!::StageLimiter
  step_limiter!::StepLimiter
end
SSPRK104(stage_limiter! = trivial_limiter!) = SSPRK104(stage_limiter!, trivial_limiter!)

"""
@article{owren1992derivation,
  title={Derivation of efficient, continuous, explicit Runge--Kutta methods},
  author={Owren, Brynjulf and Zennaro, Marino},
  journal={SIAM journal on scientific and statistical computing},
  volume={13},
  number={6},
  pages={1488--1501},
  year={1992},
  publisher={SIAM}
}
"""
struct OwrenZen3 <: OrdinaryDiffEqAdaptiveAlgorithm end
"""
@article{owren1992derivation,
  title={Derivation of efficient, continuous, explicit Runge--Kutta methods},
  author={Owren, Brynjulf and Zennaro, Marino},
  journal={SIAM journal on scientific and statistical computing},
  volume={13},
  number={6},
  pages={1488--1501},
  year={1992},
  publisher={SIAM}
}
"""
struct OwrenZen4 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
@article{owren1992derivation,
  title={Derivation of efficient, continuous, explicit Runge--Kutta methods},
  author={Owren, Brynjulf and Zennaro, Marino},
  journal={SIAM journal on scientific and statistical computing},
  volume={13},
  number={6},
  pages={1488--1501},
  year={1992},
  publisher={SIAM}
}
"""
struct OwrenZen5 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
@article{bogacki19893,
  title={A 3 (2) pair of Runge-Kutta formulas},
  author={Bogacki, Przemyslaw and Shampine, Lawrence F},
  journal={Applied Mathematics Letters},
  volume={2},
  number={4},
  pages={321--325},
  year={1989},
  publisher={Elsevier}
}
"""
struct BS3 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
@article{dormand1980family,
  title={A family of embedded Runge-Kutta formulae},
  author={Dormand, John R and Prince, Peter J},
  journal={Journal of computational and applied mathematics},
  volume={6},
  number={1},
  pages={19--26},
  year={1980},
  publisher={Elsevier}
}
"""
struct DP5 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
@article{tsitouras2011runge,
  title={Runge--Kutta pairs of order 5 (4) satisfying only the first column simplifying assumption},
  author={Tsitouras, Ch},
  journal={Computers \\& Mathematics with Applications},
  volume={62},
  number={2},
  pages={770--775},
  year={2011},
  publisher={Elsevier}
}
"""
struct Tsit5 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
Springer-Verlag.
"""
struct DP8 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
Tanaka M., Muramatsu S., Yamashita S., (1992), On the Optimization of Some Nine-Stage
Seventh-order Runge-Kutta Method, Information Processing Society of Japan,
33 (12), pp. 1512-1526.
"""
struct TanYam7 <: OrdinaryDiffEqAdaptiveAlgorithm end
struct TsitPap8 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
@article{feagin2012high,
  title={High-order explicit Runge-Kutta methods using m-symmetry},
  author={Feagin, Terry},
  year={2012},
  publisher={Neural, Parallel \\& Scientific Computations}
}
"""
struct Feagin10 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
@article{feagin2012high,
  title={High-order explicit Runge-Kutta methods using m-symmetry},
  author={Feagin, Terry},
  year={2012},
  publisher={Neural, Parallel \\& Scientific Computations}
}
"""
struct Feagin12 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
Feagin, T., “An Explicit Runge-Kutta Method of Order Fourteen,” Numerical
Algorithms, 2009
"""
struct Feagin14 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
@article{bogacki1996efficient,
  title={An efficient runge-kutta (4, 5) pair},
  author={Bogacki, P and Shampine, Lawrence F},
  journal={Computers \\& Mathematics with Applications},
  volume={32},
  number={6},
  pages={15--28},
  year={1996},
  publisher={Elsevier}
}
"""
struct BS5 <: OrdinaryDiffEqAdaptiveAlgorithm
  lazy::Bool
  BS5(;lazy=true) = new(lazy)
end

"""
@article{verner2010numerically,
  title={Numerically optimal Runge--Kutta pairs with interpolants},
  author={Verner, James H},
  journal={Numerical Algorithms},
  volume={53},
  number={2-3},
  pages={383--396},
  year={2010},
  publisher={Springer}
}
"""
struct Vern6 <: OrdinaryDiffEqAdaptiveAlgorithm
  lazy::Bool
  Vern6(;lazy=true) = new(lazy)
end

"""
@article{verner2010numerically,
  title={Numerically optimal Runge--Kutta pairs with interpolants},
  author={Verner, James H},
  journal={Numerical Algorithms},
  volume={53},
  number={2-3},
  pages={383--396},
  year={2010},
  publisher={Springer}
}
"""
struct Vern7 <: OrdinaryDiffEqAdaptiveAlgorithm
  lazy::Bool
  Vern7(;lazy=true) = new(lazy)
end

"""
@article{verner2010numerically,
  title={Numerically optimal Runge--Kutta pairs with interpolants},
  author={Verner, James H},
  journal={Numerical Algorithms},
  volume={53},
  number={2-3},
  pages={383--396},
  year={2010},
  publisher={Springer}
}
"""
struct Vern8 <: OrdinaryDiffEqAdaptiveAlgorithm
  lazy::Bool
  Vern8(;lazy=true) = new(lazy)
end

"""
@article{verner2010numerically,
  title={Numerically optimal Runge--Kutta pairs with interpolants},
  author={Verner, James H},
  journal={Numerical Algorithms},
  volume={53},
  number={2-3},
  pages={383--396},
  year={2010},
  publisher={Springer}
}
"""
struct Vern9 <: OrdinaryDiffEqAdaptiveAlgorithm
  lazy::Bool
  Vern9(;lazy=true) = new(lazy)
end
struct FRK65{T} <: OrdinaryDiffEqAdaptiveAlgorithm
  omega::T
  FRK65(omega=0.0) = new{typeof(omega)}(omega)
end

struct PFRK87{T} <: OrdinaryDiffEqAdaptiveAlgorithm
  omega::T
  PFRK87(omega=0.0) = new{typeof(omega)}(omega)
end


################################################################################

# Symplectic methods

struct SymplecticEuler <: OrdinaryDiffEqAlgorithm end

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
struct VelocityVerlet <: OrdinaryDiffEqAlgorithm end

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
struct VerletLeapfrog <: OrdinaryDiffEqAlgorithm end

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
struct PseudoVerletLeapfrog <: OrdinaryDiffEqAlgorithm end

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
struct McAte2 <: OrdinaryDiffEqAlgorithm end

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
struct Ruth3 <: OrdinaryDiffEqAlgorithm end

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
struct McAte3 <: OrdinaryDiffEqAlgorithm end

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
struct CandyRoz4 <: OrdinaryDiffEqAlgorithm end
struct McAte4 <: OrdinaryDiffEqAlgorithm end

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
struct CalvoSanz4 <: OrdinaryDiffEqAlgorithm end

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
struct McAte42 <: OrdinaryDiffEqAlgorithm end

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
struct McAte5 <: OrdinaryDiffEqAlgorithm end

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
struct Yoshida6 <: OrdinaryDiffEqAlgorithm end

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
struct KahanLi6 <: OrdinaryDiffEqAlgorithm end

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
struct McAte8 <: OrdinaryDiffEqAlgorithm end

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
struct KahanLi8 <: OrdinaryDiffEqAlgorithm end

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
struct SofSpa10 <: OrdinaryDiffEqAlgorithm end

# Nyström methods

"""
@article{rabiei2012numerical,
  title={Numerical Solution of Second-Order Ordinary Differential Equations by Improved Runge-Kutta Nystrom Method},
  author={Rabiei, Faranak and Ismail, Fudziah and Norazak, S and Emadi, Saeid},
  publisher={Citeseer}
}
"""
struct IRKN3 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
  Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
  Springer-Verlag.
"""
struct Nystrom4 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
  Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
  Springer-Verlag.
"""
struct Nystrom4VelocityIndependent <: OrdinaryDiffEqAlgorithm end

"""
@article{rabiei2012numerical,
  title={Numerical Solution of Second-Order Ordinary Differential Equations by Improved Runge-Kutta Nystrom Method},
  author={Rabiei, Faranak and Ismail, Fudziah and Norazak, S and Emadi, Saeid},
  publisher={Citeseer}
}
"""
struct IRKN4 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
  Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
  Springer-Verlag.
"""
struct Nystrom5VelocityIndependent <: OrdinaryDiffEqAlgorithm end

"""
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
struct DPRKN6 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
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
struct DPRKN8 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
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
struct DPRKN12 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
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
struct ERKN4 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
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
struct ERKN5 <: OrdinaryDiffEqAdaptiveAlgorithm end

################################################################################

# Adams Bashforth and Adams moulton methods

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct AB3 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct AB4 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct AB5 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct ABM32 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct ABM43 <: OrdinaryDiffEqAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct ABM54 <: OrdinaryDiffEqAlgorithm end


# Variable Step Size Adams methods

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct VCAB3 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct VCAB4 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct VCAB5 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct VCABM3 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct VCABM4 <: OrdinaryDiffEqAdaptiveAlgorithm end

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct VCABM5 <: OrdinaryDiffEqAdaptiveAlgorithm end

# Variable Order and Variable Step Size Adams methods

"""
E. Hairer, S. P. Norsett, G. Wanner, Solving Ordinary Differential Equations I, Nonstiff
Problems. Computational Mathematics (2nd revised ed.), Springer (1996) doi:
https://doi.org/10.1007/978-3-540-78862-1
"""
struct VCABM <: OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm end

# IMEX Multistep methods

struct CNAB2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
CNAB2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      CNAB2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct CNLF2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
CNLF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      CNLF2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct QNDF1{CS,AD,F,F2,FDT,κType} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  kappa::κType
  controller::Symbol
end
QNDF1(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                 linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                                  extrapolant=:linear,kappa = -0.1850,
                 controller = :Standard) =
                 QNDF1{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),
                 typeof(kappa)}(
                 linsolve,nlsolve,diff_type,extrapolant,kappa,controller)

QBDF1(;kwargs...) = QNDF1(;kappa=0,kwargs...)

struct QNDF2{CS,AD,F,F2,FDT,κType} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  kappa::κType
  controller::Symbol
end
QNDF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                 linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                                  extrapolant=:linear,kappa = -1//9,
                 controller = :Standard) =
                 QNDF2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),
                 typeof(kappa)}(
                 linsolve,nlsolve,diff_type,extrapolant,kappa,controller)

QBDF2(;kwargs...) = QNDF2(;kappa=0,kwargs...)


struct QNDF{MO,CS,AD,F,F2,FDT,K,T,κType} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  max_order::Val{MO}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  kappa::κType
  controller::Symbol
end
Base.@pure QNDF(;max_order::Val{MO}=Val(5),chunk_size=0,autodiff=true,diff_type=Val{:forward},
                linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),κ=nothing,tol=nothing,
                extrapolant=:linear,kappa=promote(-0.1850,-1//9,-0.0823,-0.0415,0),
                controller = :Standard) where {MO} =
                QNDF{MO,chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),
                typeof(κ),typeof(tol),typeof(kappa)}(
                max_order,linsolve,nlsolve,diff_type,κ,tol,extrapolant,kappa,controller)

Base.@pure QBDF(;kwargs...) = QNDF(;kappa=tuple(0//1,0//1,0//1,0//1,0//1),kwargs...)

"""
Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton. Implicit-Explicit Methods for Time-
Dependent Partial Differential Equations. 1995 Society for Industrial and Applied Mathematics
Journal on Numerical Analysis, 32(3), pp 797-823, 1995. doi: https://doi.org/10.1137/0732037
"""
struct SBDF{CS,AD,F,F2,FDT,K,T} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  order::Int
end
SBDF(order;chunk_size=0,autodiff=true,diff_type=Val{:forward},
     linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),κ=nothing,tol=nothing,
     extrapolant=:linear) =
     SBDF{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),
     typeof(κ),typeof(tol)}(
     linsolve,nlsolve,diff_type,κ,tol,extrapolant,order)

 """
 Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton. Implicit-Explicit Methods for Time-
 Dependent Partial Differential Equations. 1995 Society for Industrial and Applied Mathematics
 Journal on Numerical Analysis, 32(3), pp 797-823, 1995. doi: https://doi.org/10.1137/0732037
 """
IMEXEuler(;kwargs...) = SBDF(1;kwargs...)

"""
Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton. Implicit-Explicit Methods for Time-
Dependent Partial Differential Equations. 1995 Society for Industrial and Applied Mathematics
Journal on Numerical Analysis, 32(3), pp 797-823, 1995. doi: https://doi.org/10.1137/0732037
"""
SBDF2(;kwargs...) = SBDF(2;kwargs...)

"""
Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton. Implicit-Explicit Methods for Time-
Dependent Partial Differential Equations. 1995 Society for Industrial and Applied Mathematics
Journal on Numerical Analysis, 32(3), pp 797-823, 1995. doi: https://doi.org/10.1137/0732037
"""
SBDF3(;kwargs...) = SBDF(3;kwargs...)

"""
Uri M. Ascher, Steven J. Ruuth, Brian T. R. Wetton. Implicit-Explicit Methods for Time-
Dependent Partial Differential Equations. 1995 Society for Industrial and Applied Mathematics
Journal on Numerical Analysis, 32(3), pp 797-823, 1995. doi: https://doi.org/10.1137/0732037
"""
SBDF4(;kwargs...) = SBDF(4;kwargs...)

# Adams/BDF methods in Nordsieck forms
struct AN5   <: OrdinaryDiffEqAdaptiveAlgorithm end
struct JVODE{bType,aType} <: OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm
  algorithm::Symbol
  bias1::bType
  bias2::bType
  bias3::bType
  addon::aType
end

JVODE(algorithm=:Adams;bias1=6, bias2=6,bias3=10,
                 addon=1//10^6) = JVODE(algorithm,bias1,bias2,bias3,addon)
JVODE_Adams(;kwargs...) = JVODE(:Adams;kwargs...)
JVODE_BDF(;kwargs...) = JVODE(:BDF;kwargs...)

# ROCK methods

"""
Assyr Abdulle, Alexei A. Medovikov. Second Order Chebyshev Methods based on Orthogonal Polynomials.
Numerische Mathematik, 90 (1), pp 1-18, 2001. doi: https://dx.doi.org/10.1007/s002110100292
"""
struct ROCK2{E} <: OrdinaryDiffEqAdaptiveAlgorithm
  min_stages::Int
  max_stages::Int
  eigen_est::E
end
ROCK2(;min_stages=0,max_stages=200,eigen_est=nothing) = ROCK2(min_stages,max_stages,eigen_est)

"""
Assyr Abdulle. Fourth Order Chebyshev Methods With Recurrence Relation. 2002 Society for
Industrial and Applied Mathematics Journal on Scientific Computing, 23(6), pp 2041-2054, 2001.
doi: https://doi.org/10.1137/S1064827500379549
"""
struct ROCK4{E} <: OrdinaryDiffEqAdaptiveAlgorithm
  min_stages::Int
  max_stages::Int
  eigen_est::E
end
ROCK4(;min_stages=0,max_stages=152,eigen_est=nothing) = ROCK4(min_stages,max_stages,eigen_est)

# SERK methods

#=
RKC

B. P. Sommeijer, L. F. Shampine, J. G. Verwer. RKC: An Explicit Solver for Parabolic PDEs,
  Journal of Computational and Applied Mathematics, 88(2), pp 315-326, 1998. doi:
  https://doi.org/10.1016/S0377-0427(97)00219-7
=#

for Alg in [:ESERK4, :ESERK5, :RKC]
  @eval begin
    struct $Alg{E} <: OrdinaryDiffEqAdaptiveAlgorithm
      eigen_est::E
    end
    $Alg(;eigen_est=nothing) = $Alg(eigen_est)
  end
end
struct SERK2{E} <: OrdinaryDiffEqAdaptiveAlgorithm
  controller::Symbol
  eigen_est::E
end
SERK2(;controller=:PI,eigen_est=nothing) = SERK2(controller,eigen_est)

struct IRKC{CS,AD,F,F2,FDT,K,T,E} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  κ::K
  tol::T
  extrapolant::Symbol
  controller::Symbol
  eigen_est::E
end
IRKC(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                 linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),κ=nothing,tol=nothing,
                 extrapolant=:linear,controller = :Standard,eigen_est=nothing) =
  IRKC{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),typeof(κ),typeof(tol),typeof(eigen_est)}(
                 linsolve,nlsolve,diff_type,κ,tol,extrapolant,controller,eigen_est)

################################################################################

# Linear Methods

for Alg in [:MagnusMidpoint,:MagnusLeapfrog,:LieEuler,:MagnusGauss4,:MagnusNC6,:MagnusGL6,:MagnusGL8,:MagnusNC8,:MagnusGL4,:RKMK2,:RKMK4,:LieRK4,:CG2,:CG3]
  @eval struct $Alg <: OrdinaryDiffEqLinearExponentialAlgorithm
    krylov::Bool
    m::Int
    iop::Int
  end
  @eval $Alg(;krylov=false, m=30, iop=0) = $Alg(krylov, m, iop)
end

struct MagnusAdapt4 <: OrdinaryDiffEqAdaptiveAlgorithm end

struct LinearExponential <: OrdinaryDiffEqExponentialAlgorithm
  krylov::Symbol
  m::Int
  iop::Int
end
LinearExponential(;krylov=:off, m=10, iop=0) = LinearExponential(krylov, m, iop)

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
"""
struct RadauIIA3{CS,AD,F,FDT,Tol,C1,C2} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
  extrapolant::Symbol
  κ::Tol
  maxiters::Int
  fast_convergence_cutoff::C1
  new_W_γdt_cutoff::C2
  controller::Symbol
end

RadauIIA3(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                                linsolve=DEFAULT_LINSOLVE,
                                extrapolant=:dense,fast_convergence_cutoff=1//5,new_W_γdt_cutoff=1//5,
                                controller=:Predictive,κ=nothing,maxiters=10) =
                                RadauIIA3{chunk_size,autodiff,typeof(linsolve),
                                typeof(diff_type),
                                typeof(κ),typeof(fast_convergence_cutoff),typeof(new_W_γdt_cutoff)}(
                                  linsolve,diff_type,extrapolant,κ,maxiters,fast_convergence_cutoff,new_W_γdt_cutoff,controller)

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
"""
struct RadauIIA5{CS,AD,F,FDT,Tol,C1,C2} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  κ::Tol
  maxiters::Int
  fast_convergence_cutoff::C1
  new_W_γdt_cutoff::C2
  controller::Symbol
end
RadauIIA5(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                          linsolve=DEFAULT_LINSOLVE,
                          extrapolant=:dense,fast_convergence_cutoff=1//5,new_W_γdt_cutoff=1//5,
                          controller=:Predictive,κ=nothing,maxiters=10,smooth_est=true) =
                          RadauIIA5{chunk_size,autodiff,typeof(linsolve),
                          typeof(diff_type),
                          typeof(κ),typeof(fast_convergence_cutoff),typeof(new_W_γdt_cutoff)}(
                            linsolve,diff_type,smooth_est,extrapolant,κ,maxiters,fast_convergence_cutoff,new_W_γdt_cutoff,controller)

################################################################################

# SDIRK Methods

struct ImplicitEuler{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
ImplicitEuler(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                          linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          controller=:PI) =
                          ImplicitEuler{chunk_size,autodiff,typeof(linsolve),
                          typeof(nlsolve),typeof(diff_type)}(linsolve,
                          nlsolve,diff_type,extrapolant,controller)

struct ImplicitMidpoint{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
ImplicitMidpoint(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      ImplicitMidpoint{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

"""
Andre Vladimirescu. 1994. The Spice Book. John Wiley & Sons, Inc., New York,
  NY, USA.
"""
struct Trapezoid{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
Trapezoid(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                                            extrapolant=:linear,
                      controller = :PI) =
                      Trapezoid{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant,controller)

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
"""
struct TRBDF2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
TRBDF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                 linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                 smooth_est=true,extrapolant=:linear,
                 controller = :PI) =
TRBDF2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
      linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

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
"""
struct SDIRK2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
SDIRK2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 SDIRK2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct SDIRK22{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
SDIRK22(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                                            extrapolant=:linear,
                      controller = :PI) =
                      Trapezoid{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant,controller)


struct SSPSDIRK2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD} # Not adaptive
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
SSPSDIRK2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:constant,
                   controller = :PI) =
 SSPSDIRK2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

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
"""
struct Kvaerno3{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
Kvaerno3(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 Kvaerno3{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

"""
@book{kennedy2001additive,
  title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
  author={Kennedy, Christopher Alan},
  year={2001},
  publisher={National Aeronautics and Space Administration, Langley Research Center}
}
"""
struct KenCarp3{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
KenCarp3(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 KenCarp3{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct CFNLIRK3{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
CFNLIRK3(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      CFNLIRK3{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

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
"""
struct Cash4{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  embedding::Int
  controller::Symbol
end
Cash4(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI,embedding=3) =
 Cash4{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,embedding,controller)


struct SFSDIRK4{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
SFSDIRK4(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      SFSDIRK4{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

                      struct SFSDIRK5{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
                        linsolve::F
                        nlsolve::F2
                        diff_type::FDT
                        extrapolant::Symbol
                      end
SFSDIRK5(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      SFSDIRK5{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct SFSDIRK6{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end

SFSDIRK6(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      SFSDIRK6{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct SFSDIRK7{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end

SFSDIRK7(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:linear) =
                      SFSDIRK7{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

struct SFSDIRK8{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
linsolve::F
nlsolve::F2
diff_type::FDT
extrapolant::Symbol
end

SFSDIRK8(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                    linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                    extrapolant=:linear) =
                    SFSDIRK8{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                    linsolve,nlsolve,diff_type,extrapolant)

"""
E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.),
  Springer (1996)
"""
struct Hairer4{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
Hairer4(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 Hairer4{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

"""
E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.),
  Springer (1996)
"""
struct Hairer42{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
Hairer42(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 Hairer42{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

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
"""
struct Kvaerno4{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
Kvaerno4(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 Kvaerno4{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

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
"""
struct Kvaerno5{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
Kvaerno5(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 Kvaerno5{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

"""
@book{kennedy2001additive,
  title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
  author={Kennedy, Christopher Alan},
  year={2001},
  publisher={National Aeronautics and Space Administration, Langley Research Center}
}
"""
struct KenCarp4{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
KenCarp4(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 KenCarp4{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct KenCarp47{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
KenCarp47(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 KenCarp47{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

"""
@book{kennedy2001additive,
  title={Additive Runge-Kutta schemes for convection-diffusion-reaction equations},
  author={Kennedy, Christopher Alan},
  year={2001},
  publisher={National Aeronautics and Space Administration, Langley Research Center}
}
"""
struct KenCarp5{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
KenCarp5(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 KenCarp5{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

struct KenCarp58{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
KenCarp58(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:linear,
                   controller = :PI) =
 KenCarp58{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
        linsolve,nlsolve,diff_type,smooth_est,extrapolant,controller)

# `smooth_est` is not necessary, as the embedded method is also L-stable
struct ESDIRK54I8L2SA{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
ESDIRK54I8L2SA(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   extrapolant=:linear,controller = :PI) =
 ESDIRK54I8L2SA{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(linsolve,nlsolve,diff_type,extrapolant,controller)

################################################################################

# Rosenbrock Methods

#=
#### Rosenbrock23, Rosenbrock32, ode23s, ModifiedRosenbrockIntegrator

- Shampine L.F. and Reichelt M., (1997) The MATLAB ODE Suite, SIAM Journal of
Scientific Computing, 18 (1), pp. 1-22.

#### ROS3P

- Lang, J. & Verwer, ROS3P—An Accurate Third-Order Rosenbrock Solver Designed for
  Parabolic Problems J. BIT Numerical Mathematics (2001) 41: 731. doi:10.1023/A:1021900219772

#### Rodas3, Ros4LStab, Rodas4, Rodas42

- E. Hairer, G. Wanner, Solving ordinary differential equations II, stiff and
  differential-algebraic problems. Computational mathematics (2nd revised ed.), Springer (1996)

#### RosShamp4

- L. F. Shampine, Implementation of Rosenbrock Methods, ACM Transactions on
  Mathematical Software (TOMS), 8: 2, 93-113. doi:10.1145/355993.355994

#### Veldd4, Velds4

- van Veldhuizen, D-stability and Kaps-Rentrop-methods, M. Computing (1984) 32: 229.
  doi:10.1007/BF02243574

#### GRK4T, GRK4A

- Kaps, P. & Rentrop, Generalized Runge-Kutta methods of order four with stepsize control
  for stiff ordinary differential equations. P. Numer. Math. (1979) 33: 55. doi:10.1007/BF01396495

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
=#

for Alg in [:Rosenbrock23, :Rosenbrock32, :ROS3P, :Rodas3, :ROS34PW1a, :ROS34PW1b, :ROS34PW2, :ROS34PW3, :RosShamp4, :Veldd4, :Velds4, :GRK4T, :GRK4A, :Ros4LStab, :Rodas4, :Rodas42, :Rodas4P, :Rodas4P2, :Rodas5]
  @eval begin
    struct $Alg{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
      linsolve::F
      diff_type::FDT
    end
    $Alg(;chunk_size=0,autodiff=true,diff_type=Val{:forward},linsolve=DEFAULT_LINSOLVE) = $Alg{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)
  end
end

struct GeneralRosenbrock{CS,AD,F,TabType} <: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS,AD}
  tableau::TabType
  factorization::F
end

GeneralRosenbrock(;chunk_size=0,autodiff=true,
                    factorization=lu!,tableau=ROSENBROCK_DEFAULT_TABLEAU) =
                    GeneralRosenbrock{chunk_size,autodiff,typeof(factorization),typeof(tableau)}(tableau,factorization)

struct RosenbrockW6S4OS{CS,AD,F,FDT} <: OrdinaryDiffEqRosenbrockAlgorithm{CS,AD}
  linsolve::F
  diff_type::FDT
end
RosenbrockW6S4OS(;chunk_size=0,autodiff=true,diff_type=Val{:central},linsolve=DEFAULT_LINSOLVE) = RosenbrockW6S4OS{chunk_size,autodiff,typeof(linsolve),typeof(diff_type)}(linsolve,diff_type)
######################################

for Alg in [:LawsonEuler, :NorsettEuler, :ETDRK2, :ETDRK3, :ETDRK4, :HochOst4]

  """
  Hochbruck, Marlis, and Alexander Ostermann. “Exponential Integrators.” Acta
    Numerica 19 (2010): 209–86. doi:10.1017/S0962492910000048.
  """
  @eval struct $Alg{FDT} <: OrdinaryDiffEqExponentialAlgorithm
    krylov::Bool
    m::Int
    iop::Int
    autodiff::Bool
    chunksize::Int
    diff_type::FDT
  end
  @eval $Alg(;krylov=false, m=30, iop=0, autodiff=true, chunksize=0,
            diff_type = Val{:forward}) = $Alg(krylov, m, iop, autodiff,
            chunksize, diff_type)
end
const ETD1 = NorsettEuler # alias
for Alg in [:Exprb32, :Exprb43]
  @eval struct $Alg{FDT} <: OrdinaryDiffEqAdaptiveExponentialAlgorithm
    m::Int
    iop::Int
    autodiff::Bool
    chunksize::Int
    diff_type::FDT
  end
  @eval $Alg(;m=30, iop=0, autodiff=true, chunksize=0,
            diff_type = Val{:forward}) = $Alg(m, iop, autodiff, chunksize, diff_type)
end
for Alg in [:Exp4, :EPIRK4s3A, :EPIRK4s3B, :EPIRK5s3, :EXPRB53s3, :EPIRK5P1, :EPIRK5P2]
  @eval struct $Alg{FDT} <: OrdinaryDiffEqExponentialAlgorithm
    adaptive_krylov::Bool
    m::Int
    iop::Int
    autodiff::Bool
    chunksize::Int
    diff_type::FDT
  end
  @eval $Alg(;adaptive_krylov=true, m=30, iop=0, autodiff=true,
              chunksize=0, diff_type = Val{:forward}) =
              $Alg(adaptive_krylov, m, iop, autodiff, chunksize, diff_type)
end
struct SplitEuler <: OrdinaryDiffEqExponentialAlgorithm end
struct ETD2 <: OrdinaryDiffEqExponentialAlgorithm end

#########################################

"""
E. Alberdi Celayaa, J. J. Anza Aguirrezabalab, P. Chatzipantelidisc. Implementation of
an Adaptive BDF2 Formula and Comparison with The MATLAB Ode15s. Procedia Computer Science,
29, pp 1014-1026, 2014. doi: https://doi.org/10.1016/j.procs.2014.05.091
"""
struct ABDF2{CS,AD,F,F2,FDT,K,T} <: OrdinaryDiffEqNewtonAdaptiveAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  κ::K
  tol::T
  smooth_est::Bool
  extrapolant::Symbol
  controller::Symbol
end
ABDF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
      κ=nothing,tol=nothing,linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
      smooth_est=true,extrapolant=:linear,
      controller=:Standard) =
ABDF2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),
      typeof(κ),typeof(tol)}(
      linsolve,nlsolve,diff_type,κ,tol,smooth_est,extrapolant,controller)

#########################################

struct CompositeAlgorithm{T,F} <: OrdinaryDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end

################################################################################

struct MEBDF2{CS,AD,F,F2,FDT} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end
MEBDF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:constant) =
                      MEBDF2{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
                      linsolve,nlsolve,diff_type,extrapolant)

#################################################

struct PDIRK44{CS,AD,F,F2,FDT,TO} <: OrdinaryDiffEqNewtonAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  threading::TO
end
PDIRK44(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:constant,threading=true) =
                      PDIRK44{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type),typeof(threading)}(
                      linsolve,nlsolve,diff_type,extrapolant,threading)
### Algorithm Groups

const MultistepAlgorithms = Union{IRKN3,IRKN4,
                                  ABDF2,
                                  AB3,AB4,AB5,ABM32,ABM43,ABM54}

const SplitAlgorithms = Union{CNAB2,CNLF2,IRKC,SBDF,
                              KenCarp3,KenCarp4,KenCarp47,KenCarp5,KenCarp58,CFNLIRK3}


# DAE Specific Algorithms
abstract type DAEAlgorithm{CS,AD} <: DiffEqBase.AbstractDAEAlgorithm end

#=
struct DBDF{CS,AD,F,F2,FDT} <: DAEAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
end

DBDF(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
     linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),extrapolant=:linear) =
     DBDF{chunk_size,autodiff,typeof(linsolve),typeof(nlsolve),typeof(diff_type)}(
     linsolve,nlsolve,diff_type,extrapolant)
=#

struct DImplicitEuler{CS,AD,F,F2,FDT} <: DAEAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
DImplicitEuler(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                          linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          controller=:Standard) =
                          DImplicitEuler{chunk_size,autodiff,typeof(linsolve),
                          typeof(nlsolve),typeof(diff_type)}(linsolve,
                          nlsolve,diff_type,extrapolant,controller)


struct DABDF2{CS,AD,F,F2,FDT} <: DAEAlgorithm{CS,AD}
  linsolve::F
  nlsolve::F2
  diff_type::FDT
  extrapolant::Symbol
  controller::Symbol
end
DABDF2(;chunk_size=0,autodiff=true,diff_type=Val{:forward},
                          linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          controller=:Standard) =
                          DABDF2{chunk_size,autodiff,typeof(linsolve),
                          typeof(nlsolve),typeof(diff_type)}(linsolve,
                          nlsolve,diff_type,extrapolant,controller)
