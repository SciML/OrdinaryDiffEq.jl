# Initial dt selection algorithms
abstract type InitDtAlg end
struct DefaultInitDt <: InitDtAlg end
struct StiffInitDt <: InitDtAlg end

abstract type OrdinaryDiffEqAlgorithm <: SciMLBase.AbstractODEAlgorithm end
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
const NewtonAlgorithm = Union{
    OrdinaryDiffEqNewtonAlgorithm,
    OrdinaryDiffEqNewtonAdaptiveAlgorithm,
}
const RosenbrockAlgorithm = Union{
    OrdinaryDiffEqRosenbrockAlgorithm,
    OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
}

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
    nothing,
} end
const ExponentialAlgorithm = Union{
    OrdinaryDiffEqExponentialAlgorithm,
    OrdinaryDiffEqAdaptiveExponentialAlgorithm,
}

abstract type OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end

# DAE Specific Algorithms
abstract type DAEAlgorithm{CS, AD, FDT, ST, CJ} <: SciMLBase.AbstractDAEAlgorithm end

# Partitioned ODE Specific Algorithms
abstract type OrdinaryDiffEqPartitionedAlgorithm <: OrdinaryDiffEqAlgorithm end
abstract type OrdinaryDiffEqAdaptivePartitionedAlgorithm <: OrdinaryDiffEqAdaptiveAlgorithm end
const PartitionedAlgorithm = Union{
    OrdinaryDiffEqPartitionedAlgorithm,
    OrdinaryDiffEqAdaptivePartitionedAlgorithm,
}

# Second order ODE Specific Algorithms
abstract type OrdinaryDiffEqImplicitSecondOrderAlgorithm{CS, AD, FDT, ST, CJ} <:
OrdinaryDiffEqImplicitAlgorithm{CS, AD, FDT, ST, CJ} end
abstract type OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm{CS, AD, FDT, ST, CJ} <:
OrdinaryDiffEqAdaptiveImplicitAlgorithm{CS, AD, FDT, ST, CJ} end
const ImplicitSecondOrderAlgorithm = Union{
    OrdinaryDiffEqImplicitSecondOrderAlgorithm,
    OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm,
}

function SciMLBase.remake(thing::OrdinaryDiffEqAlgorithm; kwargs...)
    T = SciMLBase.remaker_of(thing)
    return T(; SciMLBase.struct_as_namedtuple(thing)..., kwargs...)
end

function SciMLBase.remake(
        thing::Union{
            OrdinaryDiffEqAdaptiveImplicitAlgorithm{
                CS, AD, FDT,
                ST, CJ,
            },
            OrdinaryDiffEqImplicitAlgorithm{
                CS, AD, FDT, ST, CJ,
            },
            DAEAlgorithm{CS, AD, FDT, ST, CJ},
        };
        kwargs...
    ) where {CS, AD, FDT, ST, CJ}
    if haskey(kwargs, :autodiff) && kwargs[:autodiff] isa AutoForwardDiff
        chunk_size = _get_fwd_chunksize(kwargs[:autodiff])
    else
        chunk_size = Val{CS}()
    end

    T = SciMLBase.remaker_of(thing)
    return T(;
        SciMLBase.struct_as_namedtuple(thing)...,
        chunk_size = chunk_size, autodiff = thing.autodiff, standardtag = Val{ST}(),
        concrete_jac = CJ === nothing ? CJ : Val{CJ}(),
        kwargs...
    )
end

###############################################################################

################################################################################

@inline trivial_limiter!(u, integrator, p, t) = nothing

################################################################################

################################################################################

################################################################################

######################################

#########################################

#########################################

"""
    CompositeAlgorithm(algs, choice_function)

A composite algorithm that chooses between multiple ODE solvers based on a user-defined choice function.
This allows for adaptive algorithm switching based on problem characteristics or performance metrics.

# Arguments

  - `algs`: Tuple or array of ODE algorithms to choose from
  - `choice_function`: Function that determines which algorithm to use at each step

The choice function receives the integrator and should return an index indicating which algorithm to use.
This enables sophisticated algorithm switching strategies based on solution behavior, step size, or other criteria.
"""
struct CompositeAlgorithm{CS, T, F} <: OrdinaryDiffEqCompositeAlgorithm
    algs::T
    choice_function::F
    function CompositeAlgorithm(algs::T, choice_function::F) where {T, F}
        CS = mapreduce(alg -> 0, max, algs)
        return new{CS, T, F}(algs, choice_function)
    end
end

@truncate_stacktrace CompositeAlgorithm 1

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
    function AutoSwitchCache(
            count::Int,
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
            current::Int = 0
        ) where {nAlg, sAlg, tolType, T}
        return new{nAlg, sAlg, tolType, T}(
            count,
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
            current
        )
    end
end

"""
    AutoSwitch(nonstiffalg, stiffalg; kwargs...)

An automatic algorithm switching method that dynamically chooses between a nonstiff and stiff solver
based on the problem's stiffness detection. This provides robust performance across a wide range of problems
without requiring the user to know the problem's stiffness characteristics a priori.

# Arguments

  - `nonstiffalg`: Algorithm to use for nonstiff regions (default: Tsit5())
  - `stiffalg`: Algorithm to use for stiff regions (default: Rodas5P())

# Keywords

  - `maxstiffstep`: Maximum number of consecutive steps before switching from nonstiff to stiff (default: 10)
  - `maxnonstiffstep`: Maximum number of consecutive steps before switching from stiff to nonstiff (default: 3)
  - `nonstifftol`: Tolerance for detecting nonstiff behavior (default: 3//4)
  - `stifftol`: Tolerance for detecting stiff behavior (default: 9//10)
  - `dtfac`: Factor for step size adjustment during switches (default: 2.0)
  - `stiffalgfirst`: Whether to start with the stiff algorithm (default: false)
  - `switch_max`: Maximum number of algorithm switches allowed (default: 10)

The switching decision is based on step size rejections and stability estimates.
"""
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
