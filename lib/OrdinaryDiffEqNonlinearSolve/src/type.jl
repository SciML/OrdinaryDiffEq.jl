# algorithms
"""
    NLFunctional(; κ = 1//100, max_iter = 10, fast_convergence_cutoff = 1//5)

Functional (fixed-point) iteration solver for the implicit stage equations,
`z ← g(z)`. No Jacobian/`W` is formed, so it is cheap per iteration but only
converges for mildly stiff problems.

# Keyword arguments

  - `κ`: relative tolerance on the increment used in the convergence test.
  - `max_iter`: maximum number of fixed-point iterations per solve.
  - `fast_convergence_cutoff`: convergence-rate threshold for fast convergence.
"""
struct NLFunctional{K, C} <: AbstractNLSolverAlgorithm
    κ::K
    fast_convergence_cutoff::C
    max_iter::Int
end

function NLFunctional(; κ = 1 // 100, max_iter = 10, fast_convergence_cutoff = 1 // 5)
    return NLFunctional(κ, fast_convergence_cutoff, max_iter)
end

"""
    NLAnderson(; κ = 1//100, max_iter = 10, max_history = 5, aa_start = 1,
               droptol = nothing, fast_convergence_cutoff = 1//5)

Anderson-accelerated fixed-point iteration for the implicit stage equations. Like
[`NLFunctional`](@ref) but mixes in `max_history` previous residuals via a
least-squares update to accelerate convergence.

# Keyword arguments

  - `κ`, `max_iter`, `fast_convergence_cutoff`: as in [`NLFunctional`](@ref).
  - `max_history`: number of past iterates kept for the acceleration.
  - `aa_start`: iteration at which acceleration starts.
  - `droptol`: optional condition-number threshold for dropping history columns.
"""
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

"""
    NLNewton(; κ = 1//100, max_iter = 10, fast_convergence_cutoff = 1//5,
             new_W_dt_cutoff = 1//5, always_new = false, check_div = true, relax = nothing)

Quasi-Newton nonlinear solver for the implicit stage equations. Uses the
`W = M/(γΔt) - J` matrix (reused/refactorized across steps and stages when
possible) to solve `g(z) = 0`.

# Keyword arguments

  - `κ`: relative tolerance on the Newton increment used in the convergence test.
  - `max_iter`: maximum number of Newton iterations per solve.
  - `fast_convergence_cutoff`: convergence-rate threshold below which convergence
    is deemed fast (allowing `W` reuse).
  - `new_W_dt_cutoff`: relative change in `γΔt` above which `W` is refactorized.
  - `always_new`: force recomputation of `W` on every solve.
  - `check_div`: enable early divergence detection.
  - `relax`: optional relaxation parameter in `[0, 1)` damping the Newton update.
"""
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
        DC,
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
    dae_jacobians::DC
end

mutable struct NLNewtonConstantCache{tType, tType2, J, W, ufType, DC} <: AbstractNLSolverCache
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
    dae_jacobians::DC
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

mutable struct NonlinearSolveCache{uType, tType, rateType, tType2, P, C, JType, WType, ufType, jcType, du1Type, weightType} <:
    AbstractNLSolverCache
    ustep::uType
    tstep::tType
    k::rateType
    atmp::uType
    invγdt::tType2
    prob::P
    cache::C
    J::JType
    W::WType
    uf::ufType
    jac_config::jcType
    du1::du1Type
    weight::weightType
    W_γdt::tType
    new_W::Bool
end
