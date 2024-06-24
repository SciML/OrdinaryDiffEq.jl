using OrdinaryDiffEq: explicit_rk_docstring
using Static: False

@inline trivial_limiter!(u, integrator, p, t) = nothing

@doc explicit_rk_docstring(
    "A second-order, five-stage explicit Runge-Kutta method for wave propagation
equations. Fixed timestep only.", "ORK256",
    references = "Matteo Bernardini, Sergio Pirozzoli.
    A General Strategy for the Optimization of Runge-Kutta Schemes for Wave
    Propagation Phenomena.
    Journal of Computational Physics, 228(11), pp 4182-4199, 2009.
    doi: https://doi.org/10.1016/j.jcp.2009.02.032",
    extra_keyword_description = """- `williamson_condition`: allows for an optimization that allows fusing broadcast expressions with the function call `f`. However, it only works for `Array` types.
                    """,
    extra_keyword_default = "williamson_condition = true")
Base.@kwdef struct ORK256{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    williamson_condition::Bool = true
end
# for backwards compatibility
function ORK256(stage_limiter!,
        step_limiter! = trivial_limiter!;
        williamson_condition = true)
    ORK256(stage_limiter!, step_limiter!, False(), williamson_condition)
end

@doc explicit_rk_docstring(
    "7-stage, third order low-storage low-dissipation, low-dispersion scheme for
discontinuous Galerkin space discretizations applied to wave propagation problems.
Optimized for PDE discretizations when maximum spatial step is small due to
geometric features of computational domain. Fixed timestep only.",
    "DGLDDRK73_C",
    references = "T. Toulorge, W. Desmet.
    Optimal Runge–Kutta Schemes for Discontinuous Galerkin Space Discretizations
    Applied to Wave Propagation Problems.
    Journal of Computational Physics, 231(4), pp 2067-2091, 2012.
    doi: https://doi.org/10.1016/j.jcp.2011.11.024",
    extra_keyword_description = """- `williamson_condition`: allows for an optimization that allows fusing broadcast expressions with the function call `f`. However, it only works for `Array` types.
                    """,
    extra_keyword_default = "williamson_condition = true")
Base.@kwdef struct DGLDDRK73_C{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    williamson_condition::Bool = true
end
# for backwards compatibility
function DGLDDRK73_C(stage_limiter!, step_limiter! = trivial_limiter!;
        williamson_condition = true)
    DGLDDRK73_C(stage_limiter!,
        step_limiter!,
        False(),
        williamson_condition)
end

@doc explicit_rk_docstring(
    "A fourth-order, five-stage explicit low-storage method of Carpenter and Kennedy
(free 3rd order Hermite interpolant). Fixed timestep only. Designed for
hyperbolic PDEs (stability properties).",
    "CarpenterKennedy2N54",
    references = "@article{carpenter1994fourth,
    title={Fourth-order 2N-storage Runge-Kutta schemes},
    author={Carpenter, Mark H and Kennedy, Christopher A},
    year={1994}
    }",
    extra_keyword_description = """- `williamson_condition`: allows for an optimization that allows fusing broadcast expressions with the function call `f`. However, it only works for `Array` types.
                    """,
    extra_keyword_default = "williamson_condition = true")
    Base.@kwdef struct CarpenterKennedy2N54{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    williamson_condition::Bool = true
end
# for backwards compatibility
function CarpenterKennedy2N54(stage_limiter!,
        step_limiter! = trivial_limiter!;
        williamson_condition = true)
    CarpenterKennedy2N54(stage_limiter!, step_limiter!, False(), williamson_condition)
end

@doc explicit_rk_docstring(
    "12-stage, fourth order low-storage method with optimized stability regions for
advection-dominated problems. Fixed timestep only.",
    "NDBLSRK124",
    references = "Jens Niegemann, Richard Diehl, Kurt Busch.
    Efficient Low-Storage Runge–Kutta Schemes with Optimized Stability Regions.
    Journal of Computational Physics, 231, pp 364-372, 2012.
    doi: https://doi.org/10.1016/j.jcp.2011.09.003",
    extra_keyword_description = """- `williamson_condition`: allows for an optimization that allows fusing broadcast expressions with the function call `f`. However, it only works for `Array` types.
                    """,
    extra_keyword_default = "williamson_condition = true")
Base.@kwdef struct NDBLSRK124{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    williamson_condition::Bool = true
end
# for backwards compatibility
function NDBLSRK124(stage_limiter!, step_limiter! = trivial_limiter!;
        williamson_condition = true)
    NDBLSRK124(stage_limiter!,
        step_limiter!, False(),
        williamson_condition)
end

@doc explicit_rk_docstring(
    "14-stage, fourth order low-storage method with optimized stability regions for
advection-dominated problems. Fixed timestep only.",
    "NDBLSRK144",
    references = "Jens Niegemann, Richard Diehl, Kurt Busch.
    Efficient Low-Storage Runge–Kutta Schemes with Optimized Stability Regions.
    Journal of Computational Physics, 231, pp 364-372, 2012.
    doi: https://doi.org/10.1016/j.jcp.2011.09.003",
    extra_keyword_description = """- `williamson_condition`: allows for an optimization that allows fusing broadcast expressions with the function call `f`. However, it only works for `Array` types.
                    """,
    extra_keyword_default = "williamson_condition = true")
Base.@kwdef struct NDBLSRK144{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    williamson_condition::Bool = true
end
# for backwards compatibility
function NDBLSRK144(stage_limiter!, step_limiter! = trivial_limiter!;
        williamson_condition = true)
    NDBLSRK144{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!, False(),
        williamson_condition)
end

@doc explicit_rk_docstring("Low-Storage Method
6-stage, fourth order low-storage, low-dissipation, low-dispersion scheme.
Fixed timestep only.", "CFRLDDRK64",
    references = "M. Calvo, J. M. Franco, L. Randez. A New Minimum Storage Runge–Kutta Scheme
    for Computational Acoustics. Journal of Computational Physics, 201, pp 1-12, 2004.
    doi: https://doi.org/10.1016/j.jcp.2004.05.012")
Base.@kwdef struct CFRLDDRK64{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CFRLDDRK64(stage_limiter!, step_limiter! = trivial_limiter!)
    CFRLDDRK64(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
7-stage, fourth order low-storage low-dissipation, low-dispersion scheme with maximal accuracy and stability limit along the imaginary axes.
Fixed timestep only.",
    "TSLDDRK74",
    references = "Kostas Tselios, T. E. Simos. Optimized Runge–Kutta Methods with Minimal Dispersion and Dissipation
    for Problems arising from Computational Acoustics. Physics Letters A, 393(1-2), pp 38-47, 2007.
    doi: https://doi.org/10.1016/j.physleta.2006.10.072")
Base.@kwdef struct TSLDDRK74{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function TSLDDRK74(stage_limiter!, step_limiter! = trivial_limiter!)
    TSLDDRK74(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "8-stage, fourth order low-storage low-dissipation, low-dispersion scheme for
discontinuous Galerkin space discretizations applied to wave propagation problems.
Optimized for PDE discretizations when maximum spatial step is small due to
geometric features of computational domain. Fixed timestep only.",
    "DGLDDRK84_C",
    references = "T. Toulorge, W. Desmet.
    Optimal Runge–Kutta Schemes for Discontinuous Galerkin Space Discretizations
    Applied to Wave Propagation Problems.
    Journal of Computational Physics, 231(4), pp 2067-2091, 2012.
    doi: https://doi.org/10.1016/j.jcp.2011.11.024",
    extra_keyword_description = """- `williamson_condition`: allows for an optimization that allows fusing broadcast expressions with the function call `f`. However, it only works for `Array` types.
                    """,
    extra_keyword_default = "williamson_condition = true")
Base.@kwdef struct DGLDDRK84_C{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    williamson_condition::Bool = true
end
# for backwards compatibility
function DGLDDRK84_C(stage_limiter!, step_limiter! = trivial_limiter!;
        williamson_condition = true)
    DGLDDRK84_C(stage_limiter!,
        step_limiter!,
        False(),
        williamson_condition)
end

@doc explicit_rk_docstring(
    "8-stage, fourth order low-storage low-dissipation, low-dispersion scheme for
discontinuous Galerkin space discretizations applied to wave propagation problems.
Optimized for PDE discretizations when the maximum spatial step size is not
constrained. Fixed timestep only.",
    "DGLDDRK84_F",
    references = "T. Toulorge, W. Desmet.
    Optimal Runge–Kutta Schemes for Discontinuous Galerkin Space Discretizations
    Applied to Wave Propagation Problems.
    Journal of Computational Physics, 231(4), pp 2067-2091, 2012.
    doi: https://doi.org/10.1016/j.jcp.2011.11.024",
    extra_keyword_description = """- `williamson_condition`: allows for an optimization that allows fusing broadcast expressions with the function call `f`. However, it only works for `Array` types.
                    """,
    extra_keyword_default = "williamson_condition = true")
Base.@kwdef struct DGLDDRK84_F{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    williamson_condition::Bool = true
end
# for backwards compatibility
function DGLDDRK84_F(stage_limiter!, step_limiter! = trivial_limiter!;
        williamson_condition = true)
    DGLDDRK84_F(stage_limiter!,
        step_limiter!,
        False(),
        williamson_condition)
end

@doc explicit_rk_docstring(
    "A fourth-order, six-stage explicit low-storage method. Fixed timestep only.",
    "SHLDDRK64",
    references = "D. Stanescu, W. G. Habashi.
    2N-Storage Low Dissipation and Dispersion Runge-Kutta Schemes for Computational
    Acoustics.
    Journal of Computational Physics, 143(2), pp 674-681, 1998.
    doi: https://doi.org/10.1006/jcph.1998.5986
    }",
    extra_keyword_description = """- `williamson_condition`: allows for an optimization that allows fusing broadcast expressions with the function call `f`. However, it only works for `Array` types.
                    """,
    extra_keyword_default = "williamson_condition = true")
Base.@kwdef struct SHLDDRK64{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    williamson_condition::Bool = true
end
# for backwards compatibility
function SHLDDRK64(stage_limiter!,
        step_limiter! = trivial_limiter!;
        williamson_condition = true)
    SHLDDRK64(stage_limiter!, step_limiter!, False(), williamson_condition)
end

@doc explicit_rk_docstring(
    "6-stage, fourth order low-stage, low-dissipation, low-dispersion scheme.
Fixed timestep only.", "RK46NL",
    references = "Julien Berland, Christophe Bogey, Christophe Bailly. Low-Dissipation and Low-Dispersion Fourth-Order Runge-Kutta Algorithm. Computers & Fluids, 35(10), pp 1459-1463, 2006. doi: https://doi.org/10.1016/j.compfluid.2005.04.003")
Base.@kwdef struct RK46NL{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function RK46NL(stage_limiter!, step_limiter! = trivial_limiter!)
    RK46NL(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
3-stage, second order (3S) low-storage scheme, optimized  the spectral difference method applied to wave propagation problems.",
    "ParsaniKetchesonDeconinck3S32",
    references = "Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
    Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems.
    SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
    doi: https://doi.org/10.1137/120885899")
Base.@kwdef struct ParsaniKetchesonDeconinck3S32{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function ParsaniKetchesonDeconinck3S32(stage_limiter!, step_limiter! = trivial_limiter!)
    ParsaniKetchesonDeconinck3S32{typeof(stage_limiter!), typeof(step_limiter!), False}(
        stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
8-stage, second order (3S) low-storage scheme, optimized for the spectral difference method applied to wave propagation problems.",
    "ParsaniKetchesonDeconinck3S82",
    references = "Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
    Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems.
    SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
    doi: https://doi.org/10.1137/120885899")
Base.@kwdef struct ParsaniKetchesonDeconinck3S82{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function ParsaniKetchesonDeconinck3S82(stage_limiter!, step_limiter! = trivial_limiter!)
    ParsaniKetchesonDeconinck3S82{typeof(stage_limiter!), typeof(step_limiter!), False}(
        stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
5-stage, third order (3S) low-storage scheme, optimized for the spectral difference method applied to wave propagation problems.",
    "ParsaniKetchesonDeconinck3S53",
    references = "Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
    Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems.
    SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
    doi: https://doi.org/10.1137/120885899")
Base.@kwdef struct ParsaniKetchesonDeconinck3S53{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function ParsaniKetchesonDeconinck3S53(stage_limiter!, step_limiter! = trivial_limiter!)
    ParsaniKetchesonDeconinck3S53{typeof(stage_limiter!), typeof(step_limiter!), False}(
        stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
17-stage, third order (3S) low-storage scheme, optimized for the spectral difference method applied to wave propagation problems.",
    "ParsaniKetchesonDeconinck3S173",
    references = "Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
    Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems.
    SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
    doi: https://doi.org/10.1137/120885899")
Base.@kwdef struct ParsaniKetchesonDeconinck3S173{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function ParsaniKetchesonDeconinck3S173(stage_limiter!, step_limiter! = trivial_limiter!)
    ParsaniKetchesonDeconinck3S173{typeof(stage_limiter!), typeof(step_limiter!), False}(
        stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
9-stage, fourth order (3S) low-storage scheme, optimized for the spectral difference method applied to wave propagation problems.",
    "ParsaniKetchesonDeconinck3S94",
    references = "Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
    Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems.
    SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
    doi: https://doi.org/10.1137/120885899")
Base.@kwdef struct ParsaniKetchesonDeconinck3S94{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function ParsaniKetchesonDeconinck3S94(stage_limiter!, step_limiter! = trivial_limiter!)
    ParsaniKetchesonDeconinck3S94{typeof(stage_limiter!), typeof(step_limiter!), False}(
        stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
18-stage, fourth order (3S) low-storage scheme, optimized for the spectral difference method applied to wave propagation problems.",
    "ParsaniKetchesonDeconinck3S184",
    references = "Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
    Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems.
    SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
    doi: https://doi.org/10.1137/120885899")
Base.@kwdef struct ParsaniKetchesonDeconinck3S184{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function ParsaniKetchesonDeconinck3S184(stage_limiter!, step_limiter! = trivial_limiter!)
    ParsaniKetchesonDeconinck3S184{typeof(stage_limiter!), typeof(step_limiter!), False}(
        stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
10-stage, fifth order (3S) low-storage scheme, optimized for the spectral difference method applied to wave propagation problems.",
    "ParsaniKetchesonDeconinck3S105",
    references = "Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
    Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems.
    SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
    doi: https://doi.org/10.1137/120885899")
Base.@kwdef struct ParsaniKetchesonDeconinck3S105{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function ParsaniKetchesonDeconinck3S105(stage_limiter!, step_limiter! = trivial_limiter!)
    ParsaniKetchesonDeconinck3S105{typeof(stage_limiter!), typeof(step_limiter!), False}(
        stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
20-stage, fifth order (3S) low-storage scheme, optimized for the spectral difference method applied to wave propagation problems.",
    "ParsaniKetchesonDeconinck3S205",
    references = "Parsani, Matteo, David I. Ketcheson, and W. Deconinck.
    Optimized explicit Runge--Kutta schemes for the spectral difference method applied to wave propagation problems.
    SIAM Journal on Scientific Computing 35.2 (2013): A957-A986.
    doi: https://doi.org/10.1137/120885899")
Base.@kwdef struct ParsaniKetchesonDeconinck3S205{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function ParsaniKetchesonDeconinck3S205(stage_limiter!, step_limiter! = trivial_limiter!)
    ParsaniKetchesonDeconinck3S205{typeof(stage_limiter!), typeof(step_limiter!), False}(
        stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
4-stage, third order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK43_2")
Base.@kwdef struct CKLLSRK43_2{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK43_2(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK43_2{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
5-stage, fourth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK54_3C")
Base.@kwdef struct CKLLSRK54_3C{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK54_3C(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK54_3C{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
9-stage, fifth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK95_4S")
Base.@kwdef struct CKLLSRK95_4S{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK95_4S(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK95_4S{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
9-stage, fifth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK95_4C")
Base.@kwdef struct CKLLSRK95_4C{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK95_4C(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK95_4C{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
9-stage, fifth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK95_4M")
Base.@kwdef struct CKLLSRK95_4M{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK95_4M(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK95_4M{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
5-stage, fourth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK54_3C_3R")
Base.@kwdef struct CKLLSRK54_3C_3R{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK54_3C_3R(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK54_3C_3R{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
5-stage, fourth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK54_3M_3R")
Base.@kwdef struct CKLLSRK54_3M_3R{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK54_3M_3R(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK54_3M_3R{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
5-stage, fourth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK54_3N_3R")
Base.@kwdef struct CKLLSRK54_3N_3R{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK54_3N_3R(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK54_3N_3R{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
8-stage, fifth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK85_4C_3R")
Base.@kwdef struct CKLLSRK85_4C_3R{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK85_4C_3R(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK85_4C_3R{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
8-stage, fifth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK85_4M_3R")
Base.@kwdef struct CKLLSRK85_4M_3R{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK85_4M_3R(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK85_4M_3R{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
8-stage, fifth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK85_4P_3R")
Base.@kwdef struct CKLLSRK85_4P_3R{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK85_4P_3R(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK85_4P_3R{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
5-stage, fourth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK54_3N_4R")
Base.@kwdef struct CKLLSRK54_3N_4R{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK54_3N_4R(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK54_3N_4R{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
5-stage, fourth order low-storage scheme, optimized for compressible Navier–Stokes equations.
", "CKLLSRK54_3M_4R")
Base.@kwdef struct CKLLSRK54_3M_4R{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK54_3M_4R(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK54_3M_4R{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "6-stage, fifth order low-storage scheme, optimized for compressible Navier–Stokes equations.",
    "CKLLSRK65_4M_4R")
Base.@kwdef struct CKLLSRK65_4M_4R{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK65_4M_4R(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK65_4M_4R(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "Low-Storage Method
8-stage, fifth order low-storage scheme, optimized for compressible Navier–Stokes equations.",
    "CKLLSRK85_4FM_4R")
Base.@kwdef struct CKLLSRK85_4FM_4R{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK85_4FM_4R(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK85_4FM_4R(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "CKLLSRK75_4M_5R: Low-Storage Method
7-stage, fifth order low-storage scheme, optimized for compressible Navier–Stokes equations.",
    "CKLLSRK75_4M_5R")
Base.@kwdef struct CKLLSRK75_4M_5R{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function CKLLSRK75_4M_5R(stage_limiter!, step_limiter! = trivial_limiter!)
    CKLLSRK75_4M_5R{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "A third-order, five-stage explicit Runge-Kutta method with embedded error estimator
designed for spectral element discretizations of compressible fluid mechanics.",
    "RDPK3Sp35",
    references = "Ranocha, Dalcin, Parsani, Ketcheson (2021)
    Optimized Runge-Kutta Methods with Automatic Step Size Control for
    Compressible Computational Fluid Dynamics
    [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)")
Base.@kwdef struct RDPK3Sp35{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function RDPK3Sp35(stage_limiter!, step_limiter! = trivial_limiter!)
    RDPK3Sp35{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A third-order, five-stage explicit Runge-Kutta method with embedded error estimator
using the FSAL property designed for spectral element discretizations of
compressible fluid mechanics.",
    "RDPK3SpFSAL35",
    references = "Ranocha, Dalcin, Parsani, Ketcheson (2021)
    Optimized Runge-Kutta Methods with Automatic Step Size Control for
    Compressible Computational Fluid Dynamics
    [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)")
Base.@kwdef struct RDPK3SpFSAL35{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function RDPK3SpFSAL35(stage_limiter!, step_limiter! = trivial_limiter!)
    RDPK3SpFSAL35{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "A fourth-order, nine-stage explicit Runge-Kutta method with embedded error estimator
designed for spectral element discretizations of compressible fluid mechanics.",
    "RDPK3Sp49",
    references = "Ranocha, Dalcin, Parsani, Ketcheson (2021)
    Optimized Runge-Kutta Methods with Automatic Step Size Control for
    Compressible Computational Fluid Dynamics
    [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)")
Base.@kwdef struct RDPK3Sp49{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function RDPK3Sp49(stage_limiter!, step_limiter! = trivial_limiter!)
    RDPK3Sp49{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A fourth-order, nine-stage explicit Runge-Kutta method with embedded error estimator
using the FSAL property designed for spectral element discretizations of
compressible fluid mechanics.",
    "RDPK3SpFSAL49",
    references = "Ranocha, Dalcin, Parsani, Ketcheson (2021)
    Optimized Runge-Kutta Methods with Automatic Step Size Control for
    Compressible Computational Fluid Dynamics
    [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)")
Base.@kwdef struct RDPK3SpFSAL49{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function RDPK3SpFSAL49(stage_limiter!, step_limiter! = trivial_limiter!)
    RDPK3SpFSAL49{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "A fifth-order, ten-stage explicit Runge-Kutta method with embedded error estimator
designed for spectral element discretizations of compressible fluid mechanics.",
    "RDPK3Sp510",
    references = "Ranocha, Dalcin, Parsani, Ketcheson (2021)
    Optimized Runge-Kutta Methods with Automatic Step Size Control for
    Compressible Computational Fluid Dynamics
    [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)")
Base.@kwdef struct RDPK3Sp510{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function RDPK3Sp510(stage_limiter!, step_limiter! = trivial_limiter!)
    RDPK3Sp510{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A fifth-order, ten-stage explicit Runge-Kutta method with embedded error estimator
using the FSAL property designed for spectral element discretizations of
compressible fluid mechanics.",
    "RDPK3SpFSAL510",
    references = "Ranocha, Dalcin, Parsani, Ketcheson (2021)
    Optimized Runge-Kutta Methods with Automatic Step Size Control for
    Compressible Computational Fluid Dynamics
    [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)")
Base.@kwdef struct RDPK3SpFSAL510{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function RDPK3SpFSAL510(stage_limiter!, step_limiter! = trivial_limiter!)
    RDPK3SpFSAL510{typeof(stage_limiter!), typeof(step_limiter!), False}(stage_limiter!,
        step_limiter!,
        False())
end

#Low Storage Explicit Runge-Kutta Methods

@doc explicit_rk_docstring("Low-Storage Method
6-stage, fourth order low-stage, low-dissipation, low-dispersion scheme.
Fixed timestep only.", "HSLDDRK64",
    references = "D. Stanescu, W. G. Habashi.
    2N-Storage Low Dissipation and Dispersion Runge-Kutta Schemes for Computational
    Acoustics.
    Journal of Computational Physics, 143(2), pp 674-681, 1998.
    doi: https://doi.org/10.1006/jcph.1998.5986
    }",
    extra_keyword_description = """- `williamson_condition`: allows for an optimization that allows fusing broadcast expressions with the function call `f`. However, it only works for `Array` types.
                    """,
    extra_keyword_default = "williamson_condition = true")
struct HSLDDRK64{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
    williamson_condition::Bool
    function HSLDDRK64(stage_limiter! = trivial_limiter!, step_limiter! = trivial_limiter!;
            williamson_condition = true)
        Base.depwarn("HSLDDRK64 is deprecated, use SHLDDRK64 instead.", :HSLDDRK64)
        SHLDDRK64(stage_limiter!, step_limiter!, thread;
            williamson_condition = williamson_condition)
    end
end

@doc explicit_rk_docstring(
    "13-stage, fourth order low-storage method with optimized stability regions for
advection-dominated problems. Fixed timestep only.",
    "NDBLSRK134",
    references = "Jens Niegemann, Richard Diehl, Kurt Busch.
    Efficient Low-Storage Runge–Kutta Schemes with Optimized Stability Regions.
    Journal of Computational Physics, 231, pp 364-372, 2012.
    doi: https://doi.org/10.1016/j.jcp.2011.09.003",
    extra_keyword_description = """- `williamson_condition`: allows for an optimization that allows fusing broadcast expressions with the function call `f`. However, it only works for `Array` types.
                    """,
    extra_keyword_default = "williamson_condition = true")
Base.@kwdef struct NDBLSRK134{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    williamson_condition::Bool = true
end
# for backwards compatibility
function NDBLSRK134(stage_limiter!, step_limiter! = trivial_limiter!;
        williamson_condition = true)
    NDBLSRK134(stage_limiter!,
        step_limiter!, False(),
        williamson_condition)
end

#SSP Optimized Runge-Kutta Methods

@doc explicit_rk_docstring("TBD",
    "KYK2014DGSSPRK_3S2")
Base.@kwdef struct KYK2014DGSSPRK_3S2{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function KYK2014DGSSPRK_3S2(stage_limiter!, step_limiter! = trivial_limiter!)
    KYK2014DGSSPRK_3S2(stage_limiter!,
        step_limiter!,
        False())
end