# algorithms
struct NLFunctional{K, C} <: AbstractNLSolverAlgorithm
    κ::K
    fast_convergence_cutoff::C
    max_iter::Int
end

function NLFunctional(; κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5)
    NLFunctional(κ, fast_convergence_cutoff, max_iter)
end

struct NLAnderson{K, D, C} <: AbstractNLSolverAlgorithm
    κ::K
    fast_convergence_cutoff::C
    max_iter::Int
    max_history::Int
    aa_start::Int
    droptol::D
end

function NLAnderson(; κ = 1 // 100, max_iter = 10, max_history::Int = 5, aa_start::Int = 1,
        droptol = nothing, fast_convergence_cutoff = 1 // 5)
    NLAnderson(κ, fast_convergence_cutoff, max_iter, max_history, aa_start, droptol)
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

function NLNewton(; κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5,
        new_W_dt_cutoff = 1 // 5, always_new = false, check_div = true,
        relax = nothing)
    if relax isa Number && !(0 <= relax < 1)
        throw(ArgumentError("The relaxation parameter must be in [0, 1), got `relax = $relax`"))
    end

    NLNewton(κ, max_iter, fast_convergence_cutoff, new_W_dt_cutoff, always_new, check_div,
        relax)
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

function NonlinearSolveAlg(alg = NewtonRaphson(autodiff = AutoFiniteDiff());
        κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5,
        new_W_dt_cutoff = 1 // 5, always_new = false, check_div = true)
    NonlinearSolveAlg(
        κ, max_iter, fast_convergence_cutoff, new_W_dt_cutoff, always_new, check_div,
        alg)
end

# solver

mutable struct NLSolver{algType, iip, uType, gamType, tmpType, tType,
    C <: AbstractNLSolverCache, E} <: AbstractNLSolver{algType, iip}
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
function NLSolver{iip, tType}(z, tmp, ztmp, γ, c, α, alg, κ, fast_convergence_cutoff, ηold,
        iter, maxiters, status, cache, method = DIRK, tmp2 = nothing,
        nfails::Int = 0) where {iip, tType}
    RT = real(eltype(z))
    NLSolver{typeof(alg), iip, typeof(z), typeof(γ), typeof(tmp2), tType, typeof(cache), RT}(
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
        one(RT))
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
    lsType
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
    new_W::Bool
end
