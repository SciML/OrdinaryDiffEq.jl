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

#### ROK4a
- Tranquilli, Paul and Sandu, Adrian (2014): 
  Rosenbrock--Krylov Methods for Large Systems of Differential Equations
  https://doi.org/10.1137/130923336

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
    :ROK4a,
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
        struct $Alg{CS, AD, F, P, FDT, ST, CJ, StepLimiter, StageLimiter} <:
               OrdinaryDiffEqRosenbrockAdaptiveAlgorithm{CS, AD, FDT, ST, CJ}
            linsolve::F
            precs::P
            step_limiter!::StepLimiter
            stage_limiter!::StageLimiter
        end
        function $Alg(; chunk_size = Val{0}(), autodiff = Val{true}(),
                standardtag = Val{true}(), concrete_jac = nothing,
                diff_type = Val{:forward}, linsolve = nothing,
                precs = DEFAULT_PRECS, step_limiter! = trivial_limiter!,
                stage_limiter! = trivial_limiter!)
            $Alg{_unwrap_val(chunk_size), _unwrap_val(autodiff), typeof(linsolve),
                typeof(precs), diff_type, _unwrap_val(standardtag),
                _unwrap_val(concrete_jac), typeof(step_limiter!),
                typeof(stage_limiter!)}(linsolve, precs, step_limiter!,
                stage_limiter!)
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

@doc rosenbrock_wanner_docstring(
    """
    A 4th order L-stable Rosenbrock-W method (fixed step only).
    """,
    "RosenbrockW6S4OS",
    references = """
    https://doi.org/10.1016/j.cam.2009.09.017
    """)
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

#########################################

struct CompositeAlgorithm{CS, T, F} <: OrdinaryDiffEqCompositeAlgorithm
    algs::T
    choice_function::F
    function CompositeAlgorithm(algs::T, choice_function::F) where {T, F}
        CS = mapreduce(alg -> has_chunksize(alg) ? get_chunksize_int(alg) : 0, max, algs)
        new{CS, T, F}(algs, choice_function)
    end
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
    function AutoSwitchCache(count::Int,
            successive_switches::Int,
            nonstiffalg::nAlg,
            stiffalg::sAlg,
            is_stiffalg::Bool,
            maxstiffstep::Int,
            maxnonstiffstep::Int,
            nonstifftol::tolType,
            stifftol::tolType,
            dtfac::T,
            stiffalgfirst::Bool,
            switch_max::Int,
            current::Int = 0) where {nAlg, sAlg, tolType, T}
        new{nAlg, sAlg, tolType, T}(count,
            successive_switches,
            nonstiffalg,
            stiffalg,
            is_stiffalg,
            maxstiffstep,
            maxnonstiffstep,
            nonstifftol,
            stifftol,
            dtfac,
            stiffalgfirst,
            switch_max,
            current)
    end
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

# ABDF2
const MultistepAlgorithms = Union{
    AB3, AB4, AB5, ABM32, ABM43, ABM54}

const SplitAlgorithms = Union{CNAB2, CNLF2}

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
