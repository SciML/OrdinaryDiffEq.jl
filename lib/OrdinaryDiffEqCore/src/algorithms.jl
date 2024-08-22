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

################################################################################

@inline trivial_limiter!(u, integrator, p, t) = nothing

################################################################################

################################################################################

################################################################################

######################################

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
