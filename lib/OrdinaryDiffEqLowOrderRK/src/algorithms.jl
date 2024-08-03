function Base.show(io::IO, alg::OrdinaryDiffEqAlgorithm)
    print(io, String(typeof(alg).name.name), "(;")
    for fieldname in fieldnames(typeof(alg))
        print(io, " ", fieldname, " = ", getfield(alg, fieldname), ",")
    end
    print(io, ")")
end

"""
Euler - The canonical forward Euler method. Fixed timestep only.
"""
struct Euler <: OrdinaryDiffEqAlgorithm end

struct SplitEuler <:
    OrdinaryDiffEqExponentialAlgorithm{0, false, Val{:forward}, Val{true}, nothing} end

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

AutoDP5(alg; kwargs...) = AutoAlgSwitch(DP5(), alg; kwargs...)

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
