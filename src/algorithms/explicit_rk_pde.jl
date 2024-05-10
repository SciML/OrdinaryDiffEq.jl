#Low Storage Explicit Runge-Kutta Methods

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

@doc explicit_rk_docstring("TBD", "SHLDDRK52")
Base.@kwdef struct SHLDDRK52{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SHLDDRK52(stage_limiter!, step_limiter! = trivial_limiter!)
    SHLDDRK52(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring("TBD", "SHLDDRK_2N")
Base.@kwdef struct SHLDDRK_2N{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SHLDDRK_2N(stage_limiter!, step_limiter! = trivial_limiter!)
    SHLDDRK_2N(stage_limiter!,
        step_limiter!,
        False())
end

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

@doc explicit_rk_docstring(
    "A second-order, two-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK22",
    references = "Shu, Chi-Wang, and Stanley Osher.
    Efficient implementation of essentially non-oscillatory shock-capturing schemes.
    Journal of Computational Physics 77.2 (1988): 439-471.
    https://doi.org/10.1016/0021-9991(88)90177-5")
Base.@kwdef struct SSPRK22{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK22(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK22(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A third-order, three-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK33",
    references = "Shu, Chi-Wang, and Stanley Osher.
    Efficient implementation of essentially non-oscillatory shock-capturing schemes.
    Journal of Computational Physics 77.2 (1988): 439-471.
    https://doi.org/10.1016/0021-9991(88)90177-5")
Base.@kwdef struct SSPRK33{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK33(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK33(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A third-order, five-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK53",
    references = "Ruuth, Steven.
    Global optimization of explicit strong-stability-preserving Runge-Kutta methods.
    Mathematics of Computation 75.253 (2006): 183-207")
Base.@kwdef struct SSPRK53{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK53(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK53(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring("TBD",
    "KYKSSPRK42")
Base.@kwdef struct KYKSSPRK42{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function KYKSSPRK42(stage_limiter!, step_limiter! = trivial_limiter!)
    KYKSSPRK42(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "A third-order, five-stage explicit strong stability preserving (SSP) low-storage method.
Fixed timestep only.",
    "SSPRK53_2N1",
    references = "Higueras and T. Roldán.
    New third order low-storage SSP explicit Runge–Kutta methods
    arXiv:1809.04807v1.")
Base.@kwdef struct SSPRK53_2N1{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK53_2N1(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK53_2N1(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "A third-order, five-stage explicit strong stability preserving (SSP) low-storage method.
Fixed timestep only.",
    "SSPRK53_2N2",
    references = "Higueras and T. Roldán.
    New third order low-storage SSP explicit Runge–Kutta methods
    arXiv:1809.04807v1.")
Base.@kwdef struct SSPRK53_2N2{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK53_2N2(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK53_2N2(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "A third-order, five-stage explicit strong stability preserving (SSP) low-storage method.
Fixed timestep only.",
    "SSPRK53_H",
    references = "Higueras and T. Roldán.
    New third order low-storage SSP explicit Runge–Kutta methods
    arXiv:1809.04807v1.")
Base.@kwdef struct SSPRK53_H{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK53_H(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK53_H(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A third-order, six-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK63",
    references = "Ruuth, Steven.
    Global optimization of explicit strong-stability-preserving Runge-Kutta methods.
    Mathematics of Computation 75.253 (2006): 183-207")
Base.@kwdef struct SSPRK63{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK63(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK63(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A third-order, seven-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK73",
    references = "Ruuth, Steven.
    Global optimization of explicit strong-stability-preserving Runge-Kutta methods.
    Mathematics of Computation 75.253 (2006): 183-207")
Base.@kwdef struct SSPRK73{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK73(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK73(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A third-order, eight-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK83",
    references = "Ruuth, Steven.
    Global optimization of explicit strong-stability-preserving Runge-Kutta methods.
    Mathematics of Computation 75.253 (2006): 183-207")
Base.@kwdef struct SSPRK83{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK83(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK83(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A third-order, four-stage explicit strong stability preserving (SSP) method.",
    "SSPRK43",
    references = """Optimal third-order explicit SSP method with four stages discovered by

    - J. F. B. M. Kraaijevanger.
      "Contractivity of Runge-Kutta methods."
      In: BIT Numerical Mathematics 31.3 (1991), pp. 482–528.
      [DOI: 10.1007/BF01933264](https://doi.org/10.1007/BF01933264).

  Embedded method constructed by

    - Sidafa Conde, Imre Fekete, John N. Shadid.
      "Embedded error estimation and adaptive step-size control for
      optimal explicit strong stability preserving Runge–Kutta methods."
      [arXiv: 1806.08693](https://arXiv.org/abs/1806.08693)

  Efficient implementation (and optimized controller) developed by

    - Hendrik Ranocha, Lisandro Dalcin, Matteo Parsani, David I. Ketcheson (2021)
      Optimized Runge-Kutta Methods with Automatic Step Size Control for
      Compressible Computational Fluid Dynamics
      [arXiv:2104.06836](https://arxiv.org/abs/2104.06836)""")
Base.@kwdef struct SSPRK43{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK43(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK43(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A third-order, four-stage explicit strong stability preserving (SSP) method.",
    "SSPRK432",
    references = "Gottlieb, Sigal, David I. Ketcheson, and Chi-Wang Shu.
    Strong stability preserving Runge-Kutta and multistep time discretizations.
    World Scientific, 2011.
    Example 6.1")
Base.@kwdef struct SSPRK432{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK432(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK432(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A third-order, four-step explicit strong stability preserving (SSP) linear multistep method.
This method does not come with an error estimator and requires a fixed time step
size.",
    "SSPRKMSVS43",
    references = "Shu, Chi-Wang.
    Total-variation-diminishing time discretizations.
    SIAM Journal on Scientific and Statistical Computing 9, no. 6 (1988): 1073-1084.
    [DOI: 10.1137/0909073](https://doi.org/10.1137/0909073)")
Base.@kwdef struct SSPRKMSVS43{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRKMSVS43(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRKMSVS43(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "A second-order, three-step explicit strong stability preserving (SSP) linear multistep method.
This method does not come with an error estimator and requires a fixed time step
size.",
    "SSPRKMSVS32",
    references = "Shu, Chi-Wang.
    Total-variation-diminishing time discretizations.
    SIAM Journal on Scientific and Statistical Computing 9, no. 6 (1988): 1073-1084.
    [DOI: 10.1137/0909073](https://doi.org/10.1137/0909073)")
Base.@kwdef struct SSPRKMSVS32{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRKMSVS32(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRKMSVS32(stage_limiter!,
        step_limiter!,
        False())
end

@doc explicit_rk_docstring(
    "A third-order, nine-stage explicit strong stability preserving (SSP) method.

Consider using `SSPRK43` instead, which uses the same main method and an
improved embedded method.",
    "SSPRK932",
    references = "Gottlieb, Sigal, David I. Ketcheson, and Chi-Wang Shu.
    Strong stability preserving Runge-Kutta and multistep time discretizations.
    World Scientific, 2011.")
Base.@kwdef struct SSPRK932{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK932(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK932(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A fourth-order, five-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK54",
    references = "Ruuth, Steven.
    Global optimization of explicit strong-stability-preserving Runge-Kutta methods.
    Mathematics of Computation 75.253 (2006): 183-207.")
Base.@kwdef struct SSPRK54{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK54(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK54(stage_limiter!,
        step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A fourth-order, ten-stage explicit strong stability preserving (SSP) method.
Fixed timestep only.",
    "SSPRK104",
    references = "Ketcheson, David I.
    Highly efficient strong stability-preserving Runge–Kutta methods with
    low-storage implementations.
    SIAM Journal on Scientific Computing 30.4 (2008): 2113-2136.")
Base.@kwdef struct SSPRK104{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function SSPRK104(stage_limiter!, step_limiter! = trivial_limiter!)
    SSPRK104(stage_limiter!,
        step_limiter!, False())
end
