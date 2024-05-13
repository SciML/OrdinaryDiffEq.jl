function Base.show(io::IO, alg::OrdinaryDiffEqAlgorithm)
    print(io, String(typeof(alg).name.name), "(;")
    for fieldname in fieldnames(typeof(alg))
        print(io, " ", fieldname, " = ", getfield(alg, fieldname), ",")
    end
    print(io, ")")
end
function explicit_rk_docstring(description::String,
        name::String;
        references::String = "",
        extra_keyword_description = "",
        extra_keyword_default = "")
    if !isempty(extra_keyword_default)
        extra_keyword_default = "\n" * repeat(" ", 8) * extra_keyword_default
    end
    start_docstring = """
       ```julia
       $name(; stage_limiter! = OrdinaryDiffEq.trivial_limiter!,
               step_limiter! = OrdinaryDiffEq.trivial_limiter!,
               thread = OrdinaryDiffEq.False(),$extra_keyword_default)
       ```

       Explicit Runge-Kutta Method.
       """
    keyword_docstring = """

        ### Keyword Arguments

         - `stage_limiter!`: function of the form `limiter!(u, integrator, p, t)`
         - `step_limiter!`: function of the form `limiter!(u, integrator, p, t)`
         - `thread`: determines whether internal broadcasting on
            appropriate CPU arrays should be serial (`thread = OrdinaryDiffEq.False()`,
            default) or use multiple threads (`thread = OrdinaryDiffEq.True()`) when
            Julia is started with multiple threads.
        """
    start_docstring * description * keyword_docstring * extra_keyword_description *
    "## References\n" * references
end

@doc explicit_rk_docstring(
    "The second order Heun's method. Uses embedded Euler method for adaptivity.",
    "Heun")
Base.@kwdef struct Heun{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function Heun(stage_limiter!, step_limiter! = trivial_limiter!)
    Heun(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "The optimized second order midpoint method. Uses embedded Euler method for adaptivity.",
    "Ralston")
Base.@kwdef struct Ralston{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function Ralston(stage_limiter!, step_limiter! = trivial_limiter!)
    Ralston(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "The second order midpoint method. Uses embedded Euler method for adaptivity.",
    "Midpoint")
Base.@kwdef struct Midpoint{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function Midpoint(stage_limiter!, step_limiter! = trivial_limiter!)
    Midpoint(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring("The canonical Runge-Kutta Order 4 method.
Uses a defect control for adaptive stepping using maximum error over the whole interval.",
    "RK4",
    references = "@article{shampine2005solving,
      title={Solving ODEs and DDEs with residual control},
      author={Shampine, LF},
      journal={Applied Numerical Mathematics},
      volume={52},
      number={1},
      pages={113--127},
      year={2005},
      publisher={Elsevier}
      }")
Base.@kwdef struct RK4{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function RK4(stage_limiter!, step_limiter! = trivial_limiter!)
    RK4(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring("TBD", "RKM")
Base.@kwdef struct RKM{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function RKM(stage_limiter!, step_limiter! = trivial_limiter!)
    RKM(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring("4-stage Pseudo-Symplectic Explicit RK method.", "3p5q(4)",
    references = "@article{Aubry1998,
    author = {A. Aubry and P. Chartier},
    journal = {BIT Numer. Math.},
    title =  {Pseudo-symplectic {R}unge-{K}utta methods},
    year = {1998},
    },
    @article{Capuano2017,
    title = {Explicit {R}unge–{K}utta schemes for incompressible flow with improved energy-conservation properties},
    journal = {J. Comput. Phys.},
    year = {2017},
    author = {F. Capuano and G. Coppola and L. Rández and L. {de Luca}},}")
Base.@kwdef struct PSRK3p5q4{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end

@doc explicit_rk_docstring("5-stage Pseudo-Symplectic Explicit RK method.", "3p6q(5)",
    references = "@article{Aubry1998,
    author = {A. Aubry and P. Chartier},
    journal = {BIT Numer. Math.},
    title =  {Pseudo-symplectic {R}unge-{K}utta methods},
    year = {1998},
    },
    @article{Capuano2017,
    title = {Explicit {R}unge–{K}utta schemes for incompressible flow with improved energy-conservation properties},
    journal = {J. Comput. Phys.},
    year = {2017},
    author = {F. Capuano and G. Coppola and L. Rández and L. {de Luca}},}")
Base.@kwdef struct PSRK3p6q5{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end

@doc explicit_rk_docstring("6-stage Pseudo-Symplectic Explicit RK method.", "4p7q(6)",
    references = "@article{Aubry1998,
    author = {A. Aubry and P. Chartier},
    journal = {BIT Numer. Math.},
    title =  {Pseudo-symplectic {R}unge-{K}utta methods},
    volume = {38},
    PAGES = {439-461},
    year = {1998},
    },
    @article{Capuano2017,
    title = {Explicit {R}unge–{K}utta schemes for incompressible flow with improved energy-conservation properties},
    journal = {J. Comput. Phys.},
    volume = {328},
    pages = {86-94},
    year = {2017},
    issn = {0021-9991},
    doi = {https://doi.org/10.1016/j.jcp.2016.10.040},
    author = {F. Capuano and G. Coppola and L. Rández and L. {de Luca}},}")
Base.@kwdef struct PSRK4p7q6{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end

@doc explicit_rk_docstring("5th order Explicit RK method.", "MSRK5",
    references = "Misha Stepanov - https://arxiv.org/pdf/2202.08443.pdf : Figure 3.")
Base.@kwdef struct MSRK5{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function MSRK5(stage_limiter!, step_limiter! = trivial_limiter!)
    MSRK5(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring("6th order Explicit RK method.", "MSRK6",
    references = "Misha Stepanov - https://arxiv.org/pdf/2202.08443.pdf : Table4")
Base.@kwdef struct MSRK6{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function MSRK6(stage_limiter!, step_limiter! = trivial_limiter!)
    MSRK6(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring("5th order Explicit RK method.",
    "Stepanov5",
    references = "@article{Stepanov2021Embedded5,
    title={Embedded (4, 5) pairs of explicit 7-stage Runge–Kutta methods with FSAL property},
    author={Misha Stepanov},
    journal={Calcolo},
    year={2021},
    volume={59}
    }")
Base.@kwdef struct Stepanov5{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function Stepanov5(stage_limiter!, step_limiter! = trivial_limiter!)
    Stepanov5(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring("4th order Runge-Kutta method designed for periodic problems.",
    "Anas5",
    extra_keyword_description = """- `w`: a periodicity estimate, which when accurate the method becomes 5th order
                              (and is otherwise 4th order with less error for better estimates).
                    """,
    extra_keyword_default = "w = 1")
Base.@kwdef struct Anas5{StageLimiter, StepLimiter, Thread, T} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    w::T = 1
end
# for backwards compatibility
function Anas5(stage_limiter!, step_limiter! = trivial_limiter!; w = 1)
    Anas5(stage_limiter!, step_limiter!, False(), w)
end

@doc explicit_rk_docstring("5th order Explicit RK method.", "RKO5",
    references = "Tsitouras, Ch. \"Explicit Runge–Kutta methods for starting integration of
    Lane–Emden problem.\" Applied Mathematics and Computation 354 (2019): 353-364.
    doi: https://doi.org/10.1016/j.amc.2019.02.047")
Base.@kwdef struct RKO65{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function RKO65(stage_limiter!, step_limiter! = trivial_limiter!)
    RKO65(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "Owren-Zennaro optimized interpolation 3/2 method (free 3rd order interpolant).",
    "OwrenZen3",
    references = "@article{owren1992derivation,
    title={Derivation of efficient, continuous, explicit Runge--Kutta methods},
    author={Owren, Brynjulf and Zennaro, Marino},
    journal={SIAM journal on scientific and statistical computing},
    volume={13},
    number={6},
    pages={1488--1501},
    year={1992},
    publisher={SIAM}
    }")
Base.@kwdef struct OwrenZen3{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function OwrenZen3(stage_limiter!, step_limiter! = trivial_limiter!)
    OwrenZen3(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "Owren-Zennaro optimized interpolation 4/3 method (free 4th order interpolant).",
    "OwrenZen4",
    references = "@article{owren1992derivation,
    title={Derivation of efficient, continuous, explicit Runge--Kutta methods},
    author={Owren, Brynjulf and Zennaro, Marino},
    journal={SIAM journal on scientific and statistical computing},
    volume={13},
    number={6},
    pages={1488--1501},
    year={1992},
    publisher={SIAM}
    }")
Base.@kwdef struct OwrenZen4{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function OwrenZen4(stage_limiter!, step_limiter! = trivial_limiter!)
    OwrenZen4(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "Owren-Zennaro optimized interpolation 5/4 method (free 5th order interpolant).",
    "OwrenZen5",
    references = "@article{owren1992derivation,
    title={Derivation of efficient, continuous, explicit Runge--Kutta methods},
    author={Owren, Brynjulf and Zennaro, Marino},
    journal={SIAM journal on scientific and statistical computing},
    volume={13},
    number={6},
    pages={1488--1501},
    year={1992},
    publisher={SIAM}
    }")
Base.@kwdef struct OwrenZen5{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function OwrenZen5(stage_limiter!, step_limiter! = trivial_limiter!)
    OwrenZen5(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A third-order, four-stage explicit FSAL Runge-Kutta method with embedded error
estimator of Bogacki and Shampine.",
    "BS3",
    references = "@article{bogacki19893,
    title={A 3 (2) pair of Runge-Kutta formulas},
    author={Bogacki, Przemyslaw and Shampine, Lawrence F},
    journal={Applied Mathematics Letters},
    volume={2},
    number={4},
    pages={321--325},
    year={1989},
    publisher={Elsevier}
    }")
Base.@kwdef struct BS3{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function BS3(stage_limiter!, step_limiter! = trivial_limiter!)
    BS3(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "Dormand-Prince's 5/4 Runge-Kutta method. (free 4th order interpolant).",
    "DP5",
    references = "@article{dormand1980family,
    title={A family of embedded Runge-Kutta formulae},
    author={Dormand, John R and Prince, Peter J},
    journal={Journal of computational and applied mathematics},
    volume={6},
    number={1},
    pages={19--26},
    year={1980},
    publisher={Elsevier}
    }")
Base.@kwdef struct DP5{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function DP5(stage_limiter!, step_limiter! = trivial_limiter!)
    DP5(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "A fifth-order explicit Runge-Kutta method with embedded error
estimator of Tsitouras. Free 4th order interpolant.", "Tsit5",
    references = "@article{tsitouras2011runge,
    title={Runge--Kutta pairs of order 5 (4) satisfying only the first column simplifying assumption},
    author={Tsitouras, Ch},
    journal={Computers \\& Mathematics with Applications},
    volume={62},
    number={2},
    pages={770--775},
    year={2011},
    publisher={Elsevier}
    }")
Base.@kwdef struct Tsit5{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
TruncatedStacktraces.@truncate_stacktrace Tsit5 3
# for backwards compatibility
function Tsit5(stage_limiter!, step_limiter! = trivial_limiter!)
    Tsit5(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "Hairer's 8/5/3 adaption of the Dormand-Prince Runge-Kutta method. (7th order interpolant).",
    "DP8",
    references = "E. Hairer, S.P. Norsett, G. Wanner, (1993) Solving Ordinary Differential Equations I.
    Nonstiff Problems. 2nd Edition. Springer Series in Computational Mathematics,
    Springer-Verlag.")
Base.@kwdef struct DP8{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function DP8(stage_limiter!, step_limiter! = trivial_limiter!)
    DP8(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring(
    "Tanaka-Yamashita 7 Runge-Kutta method. (7th order interpolant).",
    "TanYam7",
    references = "Tanaka M., Muramatsu S., Yamashita S., (1992), On the Optimization of Some Nine-Stage
    Seventh-order Runge-Kutta Method, Information Processing Society of Japan,
    33 (12), pp. 1512-1526.")
Base.@kwdef struct TanYam7{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function TanYam7(stage_limiter!, step_limiter! = trivial_limiter!)
    TanYam7(stage_limiter!, step_limiter!, False())
end

@doc explicit_rk_docstring("Tsitouras-Papakostas 8/7 Runge-Kutta method.", "TsitPap8")
Base.@kwdef struct TsitPap8{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function TsitPap8(stage_limiter!, step_limiter! = trivial_limiter!)
    TsitPap8(stage_limiter!, step_limiter!, False())
end

"""
@article{feagin2012high,
title={High-order explicit Runge-Kutta methods using m-symmetry},
author={Feagin, Terry},
year={2012},
publisher={Neural, Parallel \\& Scientific Computations}
}

Feagin10: Explicit Runge-Kutta Method
Feagin's 10th-order Runge-Kutta method.
"""
Base.@kwdef struct Feagin10{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm 
    step_limiter!::StepLimiter = trivial_limiter!
end

"""
@article{feagin2012high,
title={High-order explicit Runge-Kutta methods using m-symmetry},
author={Feagin, Terry},
year={2012},
publisher={Neural, Parallel \\& Scientific Computations}
}

Feagin12: Explicit Runge-Kutta Method
Feagin's 12th-order Runge-Kutta method.
"""
Base.@kwdef struct Feagin12{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm 
    step_limiter!::StepLimiter = trivial_limiter!
end

"""
Feagin, T., “An Explicit Runge-Kutta Method of Order Fourteen,” Numerical
Algorithms, 2009

Feagin14: Explicit Runge-Kutta Method
Feagin's 14th-order Runge-Kutta method.
"""
Base.@kwdef struct Feagin14{StepLimiter} <: OrdinaryDiffEqAdaptiveAlgorithm 
    step_limiter!::StepLimiter = trivial_limiter!
end

@doc explicit_rk_docstring("Zero Dissipation Runge-Kutta of 6th order.", "FRK65",
    extra_keyword_description = """- `omega`: a periodicity phase estimate,
                                   when accurate this method results in zero numerical dissipation.
                    """,
    extra_keyword_default = "omega = 0.0")
Base.@kwdef struct FRK65{StageLimiter, StepLimiter, Thread, T} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    omega::T = 0.0
end

# for backwards compatibility
function FRK65(stage_limiter!, step_limiter! = trivial_limiter!; omega = 0.0)
    FRK65(stage_limiter!, step_limiter!, False(), omega)
end

@doc explicit_rk_docstring("Phase-fitted Runge-Kutta of 8th order.", "PFRK87",
    extra_keyword_description = """- `omega`: a periodicity phase estimate,
                                   when accurate this method results in zero numerical dissipation.
                    """,
    extra_keyword_default = "omega = 0.0")
Base.@kwdef struct PFRK87{StageLimiter, StepLimiter, Thread, T} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    omega::T = 0.0
end
# for backwards compatibility
function PFRK87(stage_limiter!, step_limiter! = trivial_limiter!; omega = 0.0)
    PFRK87(stage_limiter!, step_limiter!, False(), omega)
end

@doc explicit_rk_docstring(
    "Bogacki-Shampine 5/4 Runge-Kutta method. (lazy 5th order interpolant).",
    "BS5",
    references = "@article{bogacki1996efficient,
    title={An efficient runge-kutta (4, 5) pair},
    author={Bogacki, P and Shampine, Lawrence F},
    journal={Computers \\& Mathematics with Applications},
    volume={32},
    number={6},
    pages={15--28},
    year={1996},
    publisher={Elsevier}
    }",
    extra_keyword_description = """- `lazy`: determines if the lazy interpolant is used.
                    """,
    extra_keyword_default = "lazy = true")
Base.@kwdef struct BS5{StageLimiter, StepLimiter, Thread} <: OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::Bool = true
end
# for backwards compatibility
function BS5(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    BS5(stage_limiter!, step_limiter!, False(), lazy)
end

@doc explicit_rk_docstring(
    "Verner's “Most Efficient” 6/5 Runge-Kutta method. (lazy 6th order interpolant).",
    "Vern6",
    references = "@article{verner2010numerically,
    title={Numerically optimal Runge--Kutta pairs with interpolants},
    author={Verner, James H},
    journal={Numerical Algorithms},
    volume={53},
    number={2-3},
    pages={383--396},
    year={2010},
    publisher={Springer}
    }",
    extra_keyword_description = """- `lazy`: determines if the lazy interpolant is used.
                    """,
    extra_keyword_default = "lazy = true")
Base.@kwdef struct Vern6{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::Bool = true
end
TruncatedStacktraces.@truncate_stacktrace Vern6 3
# for backwards compatibility
function Vern6(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    Vern6(stage_limiter!, step_limiter!, False(), lazy)
end

@doc explicit_rk_docstring(
    "Verner's “Most Efficient” 7/6 Runge-Kutta method. (lazy 7th order interpolant).",
    "Vern7",
    references = "@article{verner2010numerically,
    title={Numerically optimal Runge--Kutta pairs with interpolants},
    author={Verner, James H},
    journal={Numerical Algorithms},
    volume={53},
    number={2-3},
    pages={383--396},
    year={2010},
    publisher={Springer}
    }",
    extra_keyword_description = """- `lazy`: determines if the lazy interpolant is used.
                    """,
    extra_keyword_default = "lazy = true")
Base.@kwdef struct Vern7{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::Bool = true
end
TruncatedStacktraces.@truncate_stacktrace Vern7 3
# for backwards compatibility
function Vern7(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    Vern7(stage_limiter!, step_limiter!, False(), lazy)
end

@doc explicit_rk_docstring(
    "Verner's “Most Efficient” 8/7 Runge-Kutta method. (lazy 8th order interpolant).",
    "Vern8",
    references = "@article{verner2010numerically,
    title={Numerically optimal Runge--Kutta pairs with interpolants},
    author={Verner, James H},
    journal={Numerical Algorithms},
    volume={53},
    number={2-3},
    pages={383--396},
    year={2010},
    publisher={Springer}
    }",
    extra_keyword_description = """- `lazy`: determines if the lazy interpolant is used.
                    """,
    extra_keyword_default = "lazy = true")
Base.@kwdef struct Vern8{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::Bool = true
end
TruncatedStacktraces.@truncate_stacktrace Vern8 3
# for backwards compatibility
function Vern8(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    Vern8(stage_limiter!, step_limiter!, False(), lazy)
end

@doc explicit_rk_docstring(
    "Verner's “Most Efficient” 9/8 Runge-Kutta method. (lazy9th order interpolant).",
    "Vern9",
    references = "@article{verner2010numerically,
    title={Numerically optimal Runge--Kutta pairs with interpolants},
    author={Verner, James H},
    journal={Numerical Algorithms},
    volume={53},
    number={2-3},
    pages={383--396},
    year={2010},
    publisher={Springer}
    }",
    extra_keyword_description = """- `lazy`: determines if the lazy interpolant is used.
                    """, extra_keyword_default = "lazy = true")
Base.@kwdef struct Vern9{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
    lazy::Bool = true
end
TruncatedStacktraces.@truncate_stacktrace Vern9 3
# for backwards compatibility
function Vern9(stage_limiter!, step_limiter! = trivial_limiter!; lazy = true)
    Vern9(stage_limiter!, step_limiter!, False(), lazy)
end

"""
Euler - The canonical forward Euler method. Fixed timestep only.
"""
struct Euler <: OrdinaryDiffEqAlgorithm end

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

"""
KuttaPRK2p5: Parallel Explicit Runge-Kutta Method
A 5 parallel, 2 processor explicit Runge-Kutta method of 5th order.

These methods utilize multithreading on the f calls to parallelize the problem.
This requires that simultaneous calls to f are thread-safe.
"""
Base.@kwdef struct KuttaPRK2p5{TO} <: OrdinaryDiffEqAlgorithm
    threading::TO = true
end

@doc explicit_rk_docstring(
    "Runge–Kutta pairs of orders 9(8) for use in quadruple precision computations", "QPRK98",
    references = "Kovalnogov VN, Fedorov RV, Karpukhina TV, Simos TE, Tsitouras C. Runge–Kutta pairs 
    of orders 9 (8) for use in quadruple precision computations. Numerical Algorithms, 2023. 
    doi: https://doi.org/10.1007/s11075-023-01632-8")
Base.@kwdef struct QPRK98{StageLimiter, StepLimiter, Thread} <:
                   OrdinaryDiffEqAdaptiveAlgorithm
    stage_limiter!::StageLimiter = trivial_limiter!
    step_limiter!::StepLimiter = trivial_limiter!
    thread::Thread = False()
end
# for backwards compatibility
function QPRK98(stage_limiter!, step_limiter! = trivial_limiter!)
    QPRK98(stage_limiter!, step_limiter!, False())
end
