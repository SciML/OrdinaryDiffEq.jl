abstract type AbstractNLSolverCache end

# Method type
@enum MethodType begin
    DIRK
    COEFFICIENT_MULTISTEP
    NORDSIECK_MULTISTEP
    GLM
end

@enum NLStatus::Int8 begin
    FastConvergence = 2
    Convergence = 1
    SlowConvergence = 0
    VerySlowConvergence = -1
    Divergence = -2
end

# solver

abstract type AbstractNLSolver{algType, iip} end

mutable struct NLSolver{algType, iip, uType, gamType, tmpType, tType,
                        C <: AbstractNLSolverCache} <: AbstractNLSolver{algType, iip}
    z::uType
    tmp::uType # DIRK and multistep methods only use tmp
    tmp2::tmpType # for GLM if neccssary
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
end

# default to DIRK
function NLSolver{iip, tType}(z, tmp, ztmp, γ, c, α, alg, κ, fast_convergence_cutoff, ηold,
                              iter, maxiters, status, cache, method = DIRK, tmp2 = nothing,
                              nfails::Int = 0) where {iip, tType}
    NLSolver{typeof(alg), iip, typeof(z), typeof(γ), typeof(tmp2), tType, typeof(cache)}(z,
                                                                                         tmp,
                                                                                         tmp2,
                                                                                         ztmp,
                                                                                         γ,
                                                                                         convert(tType,
                                                                                                 c),
                                                                                         convert(tType,
                                                                                                 α),
                                                                                         alg,
                                                                                         convert(tType,
                                                                                                 κ),
                                                                                         convert(tType,
                                                                                                 fast_convergence_cutoff),
                                                                                         convert(tType,
                                                                                                 ηold),
                                                                                         iter,
                                                                                         maxiters,
                                                                                         status,
                                                                                         cache,
                                                                                         method,
                                                                                         nfails)
end

# caches

mutable struct NLNewtonCache{uType, tType, tType2, rateType, J, W, ufType, jcType, lsType
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
