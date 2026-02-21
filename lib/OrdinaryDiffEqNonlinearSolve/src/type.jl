# algorithms
struct NLFunctional{K, C} <: AbstractNLSolverAlgorithm
    κ::K
    fast_convergence_cutoff::C
    max_iter::Int
end

function NLFunctional(; κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5)
    return NLFunctional(κ, fast_convergence_cutoff, max_iter)
end

struct NLAnderson{K, D, C} <: AbstractNLSolverAlgorithm
    κ::K
    fast_convergence_cutoff::C
    max_iter::Int
    max_history::Int
    aa_start::Int
    droptol::D
end

function NLAnderson(;
        κ = 1 // 100, max_iter = 10, max_history::Int = 5, aa_start::Int = 1,
        droptol = nothing, fast_convergence_cutoff = 1 // 5
    )
    return NLAnderson(κ, fast_convergence_cutoff, max_iter, max_history, aa_start, droptol)
end

struct NLNewton{K, C1, C2, R} <: AbstractNLSolverAlgorithm
    κ::K
    max_iter::Int
    fast_convergence_cutoff::C1
    new_W_dt_cutoff::C2
    always_new::Bool
    check_div::Bool
    relax::R
end

function NLNewton(;
        κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5,
        new_W_dt_cutoff = 1 // 5, always_new = false, check_div = true,
        relax = nothing
    )
    if relax isa Number && !(0 <= relax < 1)
        throw(ArgumentError("The relaxation parameter must be in [0, 1), got `relax = $relax`"))
    end

    return NLNewton(
        κ, max_iter, fast_convergence_cutoff, new_W_dt_cutoff, always_new, check_div,
        relax
    )
end

struct NonlinearSolveAlg{K, C1, C2, A} <: AbstractNLSolverAlgorithm
    κ::K
    max_iter::Int
    fast_convergence_cutoff::C1
    new_W_dt_cutoff::C2
    always_new::Bool
    check_div::Bool
    alg::A
end

function NonlinearSolveAlg(
        alg = NewtonRaphson(autodiff = AutoFiniteDiff());
        κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5,
        new_W_dt_cutoff = 1 // 5, always_new = false, check_div = true
    )
    return NonlinearSolveAlg(
        κ, max_iter, fast_convergence_cutoff, new_W_dt_cutoff, always_new, check_div,
        alg
    )
end

# solver

mutable struct NLSolver{
        algType, iip, uType, gamType, tmpType, tType,
        C <: AbstractNLSolverCache, E,
    } <: AbstractNLSolver{algType, iip}
    z::uType
    tmp::uType # DIRK and multistep methods only use tmp
    tmp2::tmpType # for GLM if necessary
    ztmp::uType
    γ::gamType
    c::tType
    α::tType
    alg::algType
    κ::tType
    fast_convergence_cutoff::tType
    ηold::tType
    iter::Int
    maxiters::Int
    status::NLStatus
    cache::C
    method::MethodType
    nfails::Int
    prev_θ::E
end

# default to DIRK
function NLSolver{iip, tType}(
        z, tmp, ztmp, γ, c, α, alg, κ, fast_convergence_cutoff, ηold,
        iter, maxiters, status, cache, method = DIRK, tmp2 = nothing,
        nfails::Int = 0
    ) where {iip, tType}
    RT = real(eltype(z))
    return NLSolver{typeof(alg), iip, typeof(z), typeof(γ), typeof(tmp2), tType, typeof(cache), RT}(
        z,
        tmp,
        tmp2,
        ztmp,
        float(γ),
        convert(tType, c),
        convert(tType, α),
        alg,
        convert(tType, κ),
        convert(tType, fast_convergence_cutoff),
        convert(tType, ηold),
        iter,
        maxiters,
        status,
        cache,
        method,
        nfails,
        one(RT)
    )
end

# caches

mutable struct NLNewtonCache{
        uType,
        tType,
        tType2,
        rateType,
        J,
        W,
        ufType,
        jcType,
        lsType,
    } <: AbstractNLSolverCache
    ustep::uType
    tstep::tType
    k::rateType
    atmp::uType
    dz::uType
    J::J
    W::W
    new_W::Bool
    firststage::Bool
    firstcall::Bool
    W_γdt::tType
    du1::uType
    uf::ufType
    jac_config::jcType
    linsolve::lsType
    weight::uType
    invγdt::tType2
    new_W_γdt_cutoff::tType
    J_t::tType
end

mutable struct NLNewtonConstantCache{tType, tType2, J, W, ufType} <: AbstractNLSolverCache
    tstep::tType
    J::J
    W::W
    new_W::Bool
    firststage::Bool
    firstcall::Bool
    W_γdt::tType
    uf::ufType
    invγdt::tType2
    new_W_γdt_cutoff::tType
    J_t::tType
end

mutable struct NLFunctionalCache{uType, tType, rateType} <: AbstractNLSolverCache
    ustep::uType
    tstep::tType
    k::rateType
    atmp::uType
    dz::uType
end

mutable struct NLFunctionalConstantCache{tType} <: AbstractNLSolverCache
    tstep::tType
end

mutable struct NLAndersonCache{uType, tType, rateType, uEltypeNoUnits} <:
    AbstractNLSolverCache
    ustep::uType
    tstep::tType
    k::rateType
    atmp::uType
    dz::uType
    """residuals `g(zprev) - zprev` of previous fixed-point iteration"""
    dzold::uType
    """value `g(zprev)` of previous fixed-point iteration"""
    z₊old::uType
    Δz₊s::Vector{uType}
    Q::Matrix{uEltypeNoUnits}
    R::Matrix{uEltypeNoUnits}
    γs::Vector{uEltypeNoUnits}
    history::Int
    aa_start::Int
    droptol::Union{Nothing, tType}
end

mutable struct NLAndersonConstantCache{uType, tType, uEltypeNoUnits} <:
    AbstractNLSolverCache
    tstep::tType
    dz::uType
    """residuals `g(zprev) - zprev` of previous fixed-point iteration"""
    dzold::uType
    """value `g(zprev)` of previous fixed-point iteration"""
    z₊old::uType
    Δz₊s::Vector{uType}
    Q::Matrix{uEltypeNoUnits}
    R::Matrix{uEltypeNoUnits}
    γs::Vector{uEltypeNoUnits}
    history::Int
    aa_start::Int
    droptol::Union{Nothing, tType}
end

mutable struct NonlinearSolveCache{uType, tType, rateType, tType2, P, C} <:
    AbstractNLSolverCache
    ustep::uType
    tstep::tType
    k::rateType
    atmp::uType
    invγdt::tType2
    prob::P
    cache::C
end

# VF64 specialized types: fewer type parameters for shorter stack traces.

"""
    NLNewtonCacheVF64{pType, LSType, JCType}

Specialized NLNewtonCache for Vector{Float64} in-place problems with AutoSpecialize.
3 type parameters instead of 10. jac_config type varies with autodiff choice
(ForwardDiff vs FiniteDiff). linsolve varies with verbose settings.
"""
mutable struct NLNewtonCacheVF64{pType, LSType, JCType} <: AbstractNLSolverCache
    ustep::Vector{Float64}
    tstep::Float64
    k::Vector{Float64}
    atmp::Vector{Float64}
    dz::Vector{Float64}
    J::Matrix{Float64}
    W::Matrix{Float64}
    new_W::Bool
    firststage::Bool
    firstcall::Bool
    W_γdt::Float64
    du1::Vector{Float64}
    uf::OrdinaryDiffEqCore._UJacobianWrapperVF64Type{pType}
    jac_config::JCType
    linsolve::LSType
    weight::Vector{Float64}
    invγdt::Float64
    new_W_γdt_cutoff::Float64
    J_t::Float64
end

"""
    NLNewtonCacheVF64FiniteDiff{pType, LSType}

Specialized NLNewtonCache for Vector{Float64} in-place problems using AutoFiniteDiff.
2 type parameters instead of 10. The FiniteDiff jacobian prep type is function-independent,
allowing jac_config to be hardcoded. Used by DefaultODEAlgorithm which selects AutoFiniteDiff.
"""
mutable struct NLNewtonCacheVF64FiniteDiff{pType, LSType} <: AbstractNLSolverCache
    ustep::Vector{Float64}
    tstep::Float64
    k::Vector{Float64}
    atmp::Vector{Float64}
    dz::Vector{Float64}
    J::Matrix{Float64}
    W::Matrix{Float64}
    new_W::Bool
    firststage::Bool
    firstcall::Bool
    W_γdt::Float64
    du1::Vector{Float64}
    uf::OrdinaryDiffEqCore._UJacobianWrapperVF64Type{pType}
    jac_config::_JacConfigFiniteDiff
    linsolve::LSType
    weight::Vector{Float64}
    invγdt::Float64
    new_W_γdt_cutoff::Float64
    J_t::Float64
end

const NLNewtonCacheType = Union{NLNewtonCache, NLNewtonCacheVF64, NLNewtonCacheVF64FiniteDiff}

"""
    NLSolverVF64{algType, CacheType}

Specialized NLSolver for Vector{Float64} in-place problems.
2 type parameters instead of 8.
"""
mutable struct NLSolverVF64{
        algType, CacheType <: AbstractNLSolverCache,
    } <: AbstractNLSolver{algType, true}
    z::Vector{Float64}
    tmp::Vector{Float64}
    tmp2::Nothing
    ztmp::Vector{Float64}
    γ::Float64
    c::Float64
    α::Float64
    alg::algType
    κ::Float64
    fast_convergence_cutoff::Float64
    ηold::Float64
    iter::Int
    maxiters::Int
    status::NLStatus
    cache::CacheType
    method::MethodType
    nfails::Int
    prev_θ::Float64
end

const NLSolverType = Union{NLSolver, NLSolverVF64}
