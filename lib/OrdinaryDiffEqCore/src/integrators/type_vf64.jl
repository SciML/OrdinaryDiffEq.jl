# Type aliases for the AutoSpecialize FunctionWrappersWrapper with ForwardDiff AD.
# These encode the exact types used by SciMLBase's function wrapping for in-place
# Float64 ODE problems, parameterized only by the parameter type P.
const _ODEDualTag = ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}
const _ODEDualType = ForwardDiff.Dual{_ODEDualTag, Float64, 1}

const _ODEFunctionWrapperType{P} = FunctionWrappersWrappers.FunctionWrappersWrapper{
    Tuple{
        FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{Float64}, Vector{Float64}, P, Float64}},
        FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{_ODEDualType}, Vector{_ODEDualType}, P, Float64}},
        FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{_ODEDualType}, Vector{Float64}, P, _ODEDualType}},
        FunctionWrappers.FunctionWrapper{Nothing, Tuple{Vector{_ODEDualType}, Vector{_ODEDualType}, P, _ODEDualType}},
    },
    false,
}

# Full ODEFunction type for the VF64 path: in-place, AutoSpecialize, default fields.
const _ODEFunctionVF64Type{P} = SciMLBase.ODEFunction{
    true, SciMLBase.AutoSpecialize,
    _ODEFunctionWrapperType{P},
    UniformScaling{Bool},
    Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing,
    Nothing, Nothing,
    typeof(SciMLBase.DEFAULT_OBSERVED),
    Nothing, Nothing, Nothing, Nothing,
}

# ODEProblem type for the VF64 path: default kwargs, StandardODEProblem.
const _ODEProblemVF64Type{P} = SciMLBase.ODEProblem{
    Vector{Float64}, Tuple{Float64, Float64}, true, P,
    _ODEFunctionVF64Type{P},
    Base.Pairs{Symbol, Union{}, Nothing, @NamedTuple{}},
    SciMLBase.StandardODEProblem,
}

# Concrete InterpolationData for VF64: 2 type parameters instead of 7.
struct InterpolationDataVF64{pType, CacheType} <: OrdinaryDiffEqInterpolation{CacheType}
    f::_ODEFunctionVF64Type{pType}
    timeseries::Vector{Vector{Float64}}
    ts::Vector{Float64}
    ks::Vector{Vector{Vector{Float64}}}
    alg_choice::Nothing
    dense::Bool
    cache::CacheType
    differential_vars::Nothing
    sensitivitymode::Bool
end

# InterpolationDataVF64Type union is defined in interp_func.jl (after InterpolationData).

# Type aliases for the wrapper types used in Rosenbrock/NLSolver VF64 caches.
# These compute the concrete wrapper type from the parameter type P.
const _TimeGradientWrapperVF64Type{P} = DiffEqBase.TimeGradientWrapper{
    true, _ODEFunctionVF64Type{P}, Vector{Float64}, P,
}
const _UJacobianWrapperVF64Type{P} = DiffEqBase.UJacobianWrapper{
    true, _ODEFunctionVF64Type{P}, Float64, P,
}

# Concrete ODESolution for VF64: 3 type parameters instead of 16.
mutable struct ODESolutionVF64{algType, pType, CacheType} <:
    SciMLBase.AbstractODESolution{Float64, 2, Vector{Vector{Float64}}}
    u::Vector{Vector{Float64}}
    u_analytic::Nothing
    errors::Nothing
    t::Vector{Float64}
    k::Vector{Vector{Vector{Float64}}}
    discretes::Nothing
    prob::_ODEProblemVF64Type{pType}
    alg::algType
    interp::InterpolationDataVF64{pType, CacheType}
    dense::Bool
    tslocation::Int
    stats::SciMLBase.DEStats
    alg_choice::Nothing
    retcode::SciMLBase.ReturnCode.T
    resid::Nothing
    original::Nothing
    saved_subsystem::Nothing
end

const ODESolutionVF64Type = Union{SciMLBase.ODESolution, ODESolutionVF64}

# SciMLBase interface methods for ODESolutionVF64.
function SciMLBase.solution_new_retcode(sol::ODESolutionVF64, retcode)
    return @reset sol.retcode = retcode
end
function SciMLBase.solution_new_tslocation(sol::ODESolutionVF64, tslocation)
    return @reset sol.tslocation = tslocation
end

# Parametric type alias for DEOptions in the common Float64 case.
# CallbackType and SaveCacheType are left as parameters so that using saveat
# or callbacks does not force a fallback to the full 21-parameter ODEIntegrator.
const DefaultDEVerbosity = typeof(DEFAULT_VERBOSE)
const DEOptionsVF64{ControllerType, CallbackType, SaveCacheType} = DEOptions{
    Float64, Float64, Float64, Float64,
    ControllerType,
    typeof(DiffEqBase.ODE_DEFAULT_NORM),
    typeof(opnorm),
    Nothing,
    CallbackType,
    typeof(DiffEqBase.ODE_DEFAULT_ISOUTOFDOMAIN),
    typeof(DiffEqBase.ODE_DEFAULT_PROG_MESSAGE),
    typeof(DiffEqBase.ODE_DEFAULT_UNSTABLE_CHECK),
    DataStructures.BinaryHeap{Float64, DataStructures.FasterForward},
    DataStructures.BinaryHeap{Float64, DataStructures.FasterForward},
    Nothing, Nothing, Int64, Tuple{}, SaveCacheType, Tuple{},
    DefaultDEVerbosity,
}

"""
    ODEIntegratorVF64

Specialized version of [`ODEIntegrator`](@ref) for the common case of in-place
problems with `u0::Vector{Float64}` and `tspan::Tuple{Float64, Float64}` using
default options and AutoSpecialize wrapping. This type has 8 type parameters
instead of 21, resulting in significantly shorter stack traces.

The type parameters are:
- `algType`: the algorithm type
- `pType`: the parameter type (e.g. `NullParameters`, `Vector{Float64}`)
- `CacheType`: the algorithm cache type
- `FSALType`: the FSAL type (`Vector{Float64}` or `Nothing`)
- `ControllerType`: the step size controller type (e.g. `PIController{Float64}`, `DummyController`)
- `CallbackType`: the callback type in DEOptions
- `CallbackCacheType`: the callback cache type
- `SaveCacheType`: the saveat cache type in DEOptions

The `sol`, `f`, and `p` field types are all computed from the type parameters.
"""
mutable struct ODEIntegratorVF64{
        algType <: Union{OrdinaryDiffEqAlgorithm, DAEAlgorithm},
        pType, CacheType, FSALType,
        ControllerType, CallbackType, CallbackCacheType, SaveCacheType,
    } <:
    SciMLBase.AbstractODEIntegrator{algType, true, Vector{Float64}, Float64}
    sol::ODESolutionVF64{algType, pType, CacheType}
    u::Vector{Float64}
    du::Nothing
    k::Vector{Vector{Float64}}
    t::Float64
    dt::Float64
    f::_ODEFunctionVF64Type{pType}
    p::pType
    uprev::Vector{Float64}
    uprev2::Vector{Float64}
    duprev::Nothing
    tprev::Float64
    alg::algType
    dtcache::Float64
    dtchangeable::Bool
    dtpropose::Float64
    tdir::Float64
    eigen_est::Float64
    EEst::Float64
    qold::Float64
    q11::Float64
    erracc::Float64
    dtacc::Float64
    controller_cache::ControllerType
    success_iter::Int
    iter::Int
    saveiter::Int
    saveiter_dense::Int
    cache::CacheType
    callback_cache::CallbackCacheType
    kshortsize::Int
    force_stepfail::Bool
    last_stepfail::Bool
    just_hit_tstop::Bool
    do_error_check::Bool
    event_last_time::Int
    vector_event_last_time::Int
    last_event_error::Float64
    accept_step::Bool
    isout::Bool
    reeval_fsal::Bool
    u_modified::Bool
    reinitialize::Bool
    isdae::Bool
    opts::DEOptionsVF64{ControllerType, CallbackType, SaveCacheType}
    stats::SciMLBase.DEStats
    initializealg::DiffEqBase.DefaultInit
    differential_vars::Nothing
    fsalfirst::FSALType
    fsallast::FSALType
end

"""Union type for dispatching on both ODEIntegrator variants."""
const ODEIntegratorType = Union{ODEIntegrator, ODEIntegratorVF64}

# Forward controller reinit! methods for ODEIntegratorVF64.
# The originals in controllers.jl dispatch on ::ODEIntegrator; we duplicate them here
# to avoid ambiguities with SDE/DDE integrators.
SciMLBase.reinit!(integrator::ODEIntegratorVF64, controller::AbstractController) = nothing
SciMLBase.reinit!(integrator::ODEIntegratorVF64, controller::AbstractControllerCache) = nothing
SciMLBase.reinit!(integrator::ODEIntegratorVF64, cache::IControllerCache{T}) where {T} = cache.q = one(T)
function SciMLBase.reinit!(integrator::ODEIntegratorVF64, cache::PIControllerCache{T}) where {T}
    cache.q = one(T)
    cache.q11 = one(T)
    return cache.errold = T(cache.controller.qoldinit)
end
function SciMLBase.reinit!(integrator::ODEIntegratorVF64, cache::PIDControllerCache{T}) where {T}
    cache.err = MVector{3, T}(true, true, true)
    return cache.dt_factor = T(1 // 10^4)
end
function SciMLBase.reinit!(integrator::ODEIntegratorVF64, cache::PredictiveControllerCache{T}) where {T}
    cache.dtacc = one(T)
    cache.erracc = one(T)
    cache.qold = one(T)
    return cache.q = one(T)
end
