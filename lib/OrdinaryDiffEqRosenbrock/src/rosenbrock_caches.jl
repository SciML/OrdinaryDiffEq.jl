abstract type RosenbrockMutableCache <: OrdinaryDiffEqMutableCache end
abstract type RosenbrockConstantCache <: OrdinaryDiffEqConstantCache end

"""
    JacReuseState

Lightweight mutable state for tracking Jacobian reuse in Rosenbrock-W methods.
W-methods guarantee correctness with a stale Jacobian, so we can skip expensive
Jacobian recomputations when conditions allow it.

Fields:
- `last_dtgamma`: The dtgamma value from the last Jacobian computation
- `steps_since_jac`: Number of accepted steps since last Jacobian update
- `max_jac_age`: Maximum number of accepted steps between Jacobian updates (default 50)
- `cached_J`: Cached Jacobian for OOP reuse (type-erased for flexibility)
- `cached_dT`: Cached time derivative for OOP reuse
- `cached_W`: Cached factorized W matrix for OOP reuse
"""
mutable struct JacReuseState
    last_dtgamma::Float64
    pending_dtgamma::Float64
    last_naccept::Int
    max_jac_age::Int
    cached_J::Any
    cached_dT::Any
    cached_W::Any
end

function JacReuseState(dtgamma)
    return JacReuseState(Float64(dtgamma), Float64(dtgamma), 0, 50, nothing, nothing, nothing)
end

# Fake values since non-FSAL
get_fsalfirstlast(cache::RosenbrockMutableCache, u) = (nothing, nothing)

################################################################################

# Shampine's Low-order Rosenbrocks

mutable struct RosenbrockCache{
        uType, rateType, tabType, uNoUnitsType, JType, WType, TabType,
        TFType, UFType, F, JCType, GCType, RTolType, A, StepLimiter, StageLimiter,
    } <:
    RosenbrockMutableCache
    u::uType
    uprev::uType
    dense::Vector{rateType}
    du::rateType
    du1::rateType
    du2::rateType
    dtC::Matrix{tabType}
    dtd::Vector{tabType}
    ks::Vector{rateType}
    fsalfirst::rateType
    fsallast::rateType
    dT::rateType
    J::JType
    W::WType
    tmp::rateType
    atmp::uNoUnitsType
    weight::uNoUnitsType
    tab::TabType
    tf::TFType
    uf::UFType
    linsolve_tmp::rateType
    linsolve::F
    jac_config::JCType
    grad_config::GCType
    reltol::RTolType
    alg::A
    step_limiter!::StepLimiter
    stage_limiter!::StageLimiter
    interp_order::Int
    jac_reuse::JacReuseState
end
function full_cache(c::RosenbrockCache)
    return [
        c.u, c.uprev, c.dense..., c.du, c.du1, c.du2,
        c.ks..., c.fsalfirst, c.fsallast, c.dT, c.tmp, c.atmp, c.weight, c.linsolve_tmp,
    ]
end

struct RosenbrockCombinedConstantCache{TF, UF, Tab, JType, WType, F, AD} <:
    RosenbrockConstantCache
    tf::TF
    uf::UF
    tab::Tab
    J::JType
    W::WType
    linsolve::F
    autodiff::AD
    interp_order::Int
    jac_reuse::JacReuseState
end

@cache mutable struct Rosenbrock23Cache{
        uType, rateType, uNoUnitsType, JType, WType,
        TabType, TFType, UFType, F, JCType, GCType,
        RTolType, A, AV, StepLimiter, StageLimiter,
    } <: RosenbrockMutableCache
    u::uType
    uprev::uType
    k₁::rateType
    k₂::rateType
    k₃::rateType
    du1::rateType
    du2::rateType
    f₁::rateType
    fsalfirst::rateType
    fsallast::rateType
    dT::rateType
    J::JType
    W::WType
    tmp::rateType
    atmp::uNoUnitsType
    weight::uNoUnitsType
    tab::TabType
    tf::TFType
    uf::UFType
    linsolve_tmp::rateType
    linsolve::F
    jac_config::JCType
    grad_config::GCType
    reltol::RTolType
    alg::A
    algebraic_vars::AV
    step_limiter!::StepLimiter
    stage_limiter!::StageLimiter
    jac_reuse::JacReuseState
end

@cache mutable struct Rosenbrock32Cache{
        uType, rateType, uNoUnitsType, JType, WType,
        TabType, TFType, UFType, F, JCType, GCType,
        RTolType, A, AV, StepLimiter, StageLimiter,
    } <: RosenbrockMutableCache
    u::uType
    uprev::uType
    k₁::rateType
    k₂::rateType
    k₃::rateType
    du1::rateType
    du2::rateType
    f₁::rateType
    fsalfirst::rateType
    fsallast::rateType
    dT::rateType
    J::JType
    W::WType
    tmp::rateType
    atmp::uNoUnitsType
    weight::uNoUnitsType
    tab::TabType
    tf::TFType
    uf::UFType
    linsolve_tmp::rateType
    linsolve::F
    jac_config::JCType
    grad_config::GCType
    reltol::RTolType
    alg::A
    algebraic_vars::AV
    step_limiter!::StepLimiter
    stage_limiter!::StageLimiter
    jac_reuse::JacReuseState
end

function alg_cache(
        alg::Rosenbrock23, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k₁ = zero(rate_prototype)
    k₂ = zero(rate_prototype)
    k₃ = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    # f₀ = zero(u) fsalfirst
    f₁ = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rosenbrock23Tableau(constvalue(uBottomEltypeNoUnits))
    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)
    linsolve_tmp = zero(rate_prototype)

    grad_config = build_grad_config(alg, f, tf, du1, t)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, du2)

    J, W = build_J_W(alg, u, uprev, p, t, dt, f, jac_config, uEltypeNoUnits, Val(true))

    linprob = LinearProblem(W, _vec(linsolve_tmp); u0 = _vec(tmp))
    Pl,
        Pr = wrapprecs(
        alg.precs(
            W, nothing, u, p, t, nothing, nothing, nothing,
            nothing
        )..., weight, tmp
    )
    linsolve = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        Pl = Pl, Pr = Pr,
        assumptions = LinearSolve.OperatorAssumptions(true),
        verbose = verbose.linear_verbosity
    )

    algebraic_vars = f.mass_matrix === I ? nothing :
        [all(iszero, x) for x in eachcol(f.mass_matrix)]

    jac_reuse = JacReuseState(zero(dt))

    return Rosenbrock23Cache(
        u, uprev, k₁, k₂, k₃, du1, du2, f₁,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf,
        linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, algebraic_vars, alg.step_limiter!,
        alg.stage_limiter!, jac_reuse
    )
end

function alg_cache(
        alg::Rosenbrock32, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k₁ = zero(rate_prototype)
    k₂ = zero(rate_prototype)
    k₃ = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    # f₀ = zero(u) fsalfirst
    f₁ = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rosenbrock32Tableau(constvalue(uBottomEltypeNoUnits))

    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)
    linsolve_tmp = zero(rate_prototype)

    grad_config = build_grad_config(alg, f, tf, du1, t)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, du2)

    J, W = build_J_W(alg, u, uprev, p, t, dt, f, jac_config, uEltypeNoUnits, Val(true))

    linprob = LinearProblem(W, _vec(linsolve_tmp); u0 = _vec(tmp))

    Pl,
        Pr = wrapprecs(
        alg.precs(
            W, nothing, u, p, t, nothing, nothing, nothing,
            nothing
        )..., weight, tmp
    )
    linsolve = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        Pl = Pl, Pr = Pr,
        assumptions = LinearSolve.OperatorAssumptions(true),
        verbose = verbose.linear_verbosity
    )

    algebraic_vars = f.mass_matrix === I ? nothing :
        [all(iszero, x) for x in eachcol(f.mass_matrix)]

    jac_reuse = JacReuseState(zero(dt))

    return Rosenbrock32Cache(
        u, uprev, k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W,
        tmp, atmp, weight, tab, tf, uf, linsolve_tmp, linsolve, jac_config,
        grad_config, reltol, alg, algebraic_vars, alg.step_limiter!, alg.stage_limiter!,
        jac_reuse
    )
end

struct Rosenbrock23ConstantCache{T, TF, UF, JType, WType, F, AD} <:
    RosenbrockConstantCache
    c₃₂::T
    d::T
    tf::TF
    uf::UF
    J::JType
    W::WType
    linsolve::F
    autodiff::AD
    jac_reuse::JacReuseState
end

function Rosenbrock23ConstantCache(
        ::Type{T}, tf, uf, J, W, linsolve, autodiff, jac_reuse
    ) where {T}
    tab = Rosenbrock23Tableau(T)
    return Rosenbrock23ConstantCache(
        tab.c₃₂, tab.d, tf, uf, J, W, linsolve, autodiff, jac_reuse
    )
end

function alg_cache(
        alg::Rosenbrock23, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, nothing, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    jac_reuse = JacReuseState(zero(dt))
    return Rosenbrock23ConstantCache(
        constvalue(uBottomEltypeNoUnits), tf, uf, J, W, linsolve,
        alg_autodiff(alg), jac_reuse
    )
end

struct Rosenbrock32ConstantCache{T, TF, UF, JType, WType, F, AD} <:
    RosenbrockConstantCache
    c₃₂::T
    d::T
    tf::TF
    uf::UF
    J::JType
    W::WType
    linsolve::F
    autodiff::AD
    jac_reuse::JacReuseState
end

function Rosenbrock32ConstantCache(
        ::Type{T}, tf, uf, J, W, linsolve, autodiff, jac_reuse
    ) where {T}
    tab = Rosenbrock32Tableau(T)
    return Rosenbrock32ConstantCache(
        tab.c₃₂, tab.d, tf, uf, J, W, linsolve, autodiff, jac_reuse
    )
end

function alg_cache(
        alg::Rosenbrock32, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, nothing, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    jac_reuse = JacReuseState(zero(dt))
    return Rosenbrock32ConstantCache(
        constvalue(uBottomEltypeNoUnits), tf, uf, J, W, linsolve,
        alg_autodiff(alg), jac_reuse
    )
end

### Rodas4+ methods and consolidated Rosenbrock methods (using RodasTableau)

# Helper accessors for step_limiter!/stage_limiter! — algorithms that have these fields
# return them directly; algorithms without return trivial_limiter!.
_get_step_limiter(alg) = trivial_limiter!
_get_stage_limiter(alg) = trivial_limiter!
for Alg in (
        :Rosenbrock23, :Rosenbrock32, :ROS3P, :Rodas3, :Rodas23W, :Rodas3P,
        :Rodas4, :Rodas42, :Rodas4P, :Rodas4P2, :Rodas5, :Rodas5P,
        :Rodas5Pe, :Rodas5Pr, :Rodas6P,
    )
    @eval _get_step_limiter(alg::$Alg) = alg.step_limiter!
    @eval _get_stage_limiter(alg::$Alg) = alg.stage_limiter!
end

# Tableau type dispatch
tabtype(::Rodas4) = Rodas4Tableau
tabtype(::Rodas42) = Rodas42Tableau
tabtype(::Rodas4P) = Rodas4PTableau
tabtype(::Rodas4P2) = Rodas4P2Tableau
tabtype(::Rodas5) = Rodas5Tableau
tabtype(::Rodas5P) = Rodas5PTableau
tabtype(::Rodas5Pr) = Rodas5PTableau
tabtype(::Rodas5Pe) = Rodas5PTableau
tabtype(::Rodas6P) = Rodas6PTableau

# Consolidated methods: tableau type dispatch
tabtype(::ROS3P) = ROS3PRodasTableau
tabtype(::Rodas3) = Rodas3RodasTableau
tabtype(::Rodas3P) = Rodas3PRodasTableau
tabtype(::Rodas23W) = Rodas23WRodasTableau
tabtype(::ROS2) = ROS2RodasTableau
tabtype(::ROS2PR) = ROS2PRRodasTableau
tabtype(::ROS2S) = ROS2SRodasTableau
tabtype(::ROS3) = ROS3RodasTableau
tabtype(::ROS3PR) = ROS3PRRodasTableau
tabtype(::Scholz4_7) = Scholz4_7RodasTableau
tabtype(::ROS34PW1a) = ROS34PW1aRodasTableau
tabtype(::ROS34PW1b) = ROS34PW1bRodasTableau
tabtype(::ROS34PW2) = ROS34PW2RodasTableau
tabtype(::ROS34PW3) = ROS34PW3RodasTableau
tabtype(::ROS34PRw) = ROS34PRwRodasTableau
tabtype(::ROS3PRL) = ROS3PRLRodasTableau
tabtype(::ROS3PRL2) = ROS3PRL2RodasTableau
tabtype(::ROK4a) = ROK4aRodasTableau
tabtype(::RosShamp4) = RosShamp4RodasTableau
tabtype(::Veldd4) = Veldd4RodasTableau
tabtype(::Velds4) = Velds4RodasTableau
tabtype(::GRK4T) = GRK4TRodasTableau
tabtype(::GRK4A) = GRK4ARodasTableau
tabtype(::Ros4LStab) = Ros4LStabRodasTableau
tabtype(::RosenbrockW6S4OS) = RosenbrockW6S4OSRodasTableau

# Union of all algorithms using RodasTableau-based RosenbrockCache
const RodasTableauAlgorithms = Union{
    Rodas4, Rodas42, Rodas4P, Rodas4P2, Rodas5,
    Rodas5P, Rodas5Pe, Rodas5Pr, Rodas6P,
    ROS3P, Rodas3, Rodas3P, Rodas23W,
    ROS2, ROS2PR, ROS2S, ROS3, ROS3PR, Scholz4_7,
    ROS34PW1a, ROS34PW1b, ROS34PW2, ROS34PW3,
    ROS34PRw, ROS3PRL, ROS3PRL2, ROK4a,
    RosShamp4, Veldd4, Velds4, GRK4T, GRK4A, Ros4LStab,
    RosenbrockW6S4OS,
}

function alg_cache(
        alg::RodasTableauAlgorithms,
        u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, nothing, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    tab = tabtype(alg)(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    jac_reuse = JacReuseState(zero(dt))
    H_rows = size(tab.H, 1)
    # Rodas3P/Rodas23W: H has 3 rows but only 2 are for interpolation;
    # the 3rd row is for interpoldiff error estimation
    if alg isa Union{Rodas3P, Rodas23W}
        interp_order = 2
    else
        interp_order = H_rows > 0 ? H_rows : 2
    end
    return RosenbrockCombinedConstantCache(
        tf, uf,
        tab, J, W, linsolve,
        alg_autodiff(alg), interp_order, jac_reuse
    )
end

function alg_cache(
        alg::RodasTableauAlgorithms,
        u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = tabtype(alg)(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    # Initialize vectors: kshortsize depends on whether H has rows
    H_rows = size(tab.H, 1)
    kshortsize = H_rows > 0 ? H_rows : 2
    # Rodas3P/Rodas23W: H has 3 rows but only 2 are for interpolation;
    # the 3rd row is for interpoldiff error estimation
    if alg isa Union{Rodas3P, Rodas23W}
        interp_order = 2
    else
        interp_order = kshortsize
    end
    dense = [zero(rate_prototype) for _ in 1:kshortsize]
    ks = [zero(rate_prototype) for _ in 1:size(tab.A, 1)]
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)

    # Promote t-type for AD
    dtC = zero(tab.C) .* dt
    dtd = zero(tab.d) .* dt

    # Initialize other variables
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)

    # Temporary and helper variables
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)

    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)

    grad_config = build_grad_config(alg, f, tf, du1, t)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, du2)

    J, W = build_J_W(alg, u, uprev, p, t, dt, f, jac_config, uEltypeNoUnits, Val(true))

    Pl,
        Pr = wrapprecs(
        alg.precs(
            W, nothing, u, p, t, nothing, nothing, nothing,
            nothing
        )..., weight, tmp
    )

    linsolve_tmp = zero(rate_prototype)
    linprob = LinearProblem(W, _vec(linsolve_tmp); u0 = _vec(tmp))

    linsolve = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        Pl = Pl, Pr = Pr,
        assumptions = LinearSolve.OperatorAssumptions(true),
        verbose = verbose.linear_verbosity
    )

    jac_reuse = JacReuseState(zero(dt))

    # Return the cache struct with vectors
    return RosenbrockCache(
        u, uprev, dense, du, du1, du2, dtC, dtd, ks, fsalfirst, fsallast,
        dT, J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg,
        _get_step_limiter(alg), _get_stage_limiter(alg), interp_order, jac_reuse
    )
end

function get_fsalfirstlast(
        cache::Union{
            Rosenbrock23Cache, Rosenbrock32Cache,
            RosenbrockCache,
        },
        u
    )
    return (cache.fsalfirst, cache.fsallast)
end

################################################################################

### Tsit5DA - hybrid explicit/linear-implicit method for DAEs

struct HybridExplicitImplicitConstantCache{TF, UF, Tab, JType, WType, F, AD} <: RosenbrockConstantCache
    tf::TF
    uf::UF
    tab::Tab
    J::JType
    W::WType
    linsolve::F
    autodiff::AD
    interp_order::Int
end

mutable struct HybridExplicitImplicitCache{
        uType, rateType, uNoUnitsType, JType, WType, TabType,
        TFType, UFType, F, JCType, GCType, RTolType, A,
        StepLimiter, StageLimiter, DVType, AVType,
        GZType, GYType, WZType, FZ,
    } <: RosenbrockMutableCache
    u::uType
    uprev::uType
    dense::Vector{rateType}
    du::rateType
    du1::rateType
    du2::rateType
    ks::Vector{rateType}
    fsalfirst::rateType
    fsallast::rateType
    dT::rateType
    J::JType
    W::WType
    tmp::rateType
    atmp::uNoUnitsType
    weight::uNoUnitsType
    tab::TabType
    tf::TFType
    uf::UFType
    linsolve_tmp::rateType
    linsolve::F
    jac_config::JCType
    grad_config::GCType
    reltol::RTolType
    alg::A
    step_limiter!::StepLimiter
    stage_limiter!::StageLimiter
    interp_order::Int
    # DAE-specific fields
    diff_vars::DVType
    alg_vars::AVType
    g_z::GZType         # n_g x n_g algebraic Jacobian block
    g_y::GYType         # n_g x n_f coupling block
    W_z::WZType         # -gamma * g_z (used for linear solve)
    linsolve_tmp_z::FZ  # n_g-sized RHS for algebraic solve
end
function full_cache(c::HybridExplicitImplicitCache)
    return [
        c.u, c.uprev, c.dense..., c.du, c.du1, c.du2,
        c.ks..., c.fsalfirst, c.fsallast, c.dT, c.tmp, c.atmp, c.weight, c.linsolve_tmp,
    ]
end

function alg_cache(
        alg::HybridExplicitImplicitRK, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, nothing, uEltypeNoUnits, Val(false))
    tab = alg.tab(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return HybridExplicitImplicitConstantCache(
        tf, uf, tab, J, W, nothing, alg_autodiff(alg), size(tab.H, 1)
    )
end

function alg_cache(
        alg::HybridExplicitImplicitRK, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = alg.tab(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    num_stages = size(tab.A, 1)
    interp_order = size(tab.H, 1)

    # Initialize vectors
    dense = [zero(rate_prototype) for _ in 1:interp_order]
    ks = [zero(rate_prototype) for _ in 1:num_stages]
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    linsolve_tmp = zero(rate_prototype)

    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)

    grad_config = build_grad_config(alg, f, tf, du1, t)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, du2)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, jac_config, uEltypeNoUnits, Val(true))

    linprob = LinearProblem(W, _vec(linsolve_tmp); u0 = _vec(tmp))
    Pl, Pr = wrapprecs(
        alg.precs(
            W, nothing, u, p, t, nothing, nothing, nothing,
            nothing
        )..., weight, tmp
    )
    linsolve = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        Pl = Pl, Pr = Pr,
        assumptions = LinearSolve.OperatorAssumptions(true),
        verbose = verbose.linear_verbosity
    )

    # Detect algebraic variables from mass matrix
    mass_matrix = f.mass_matrix
    if mass_matrix === I
        diff_vars = collect(1:length(u))
        alg_vars = Int[]
        g_z = zeros(eltype(u), 0, 0)
        g_y = zeros(eltype(u), 0, 0)
        W_z = zeros(eltype(u), 0, 0)
        linsolve_tmp_z = zeros(eltype(u), 0)
        linsolve_z = nothing
    else
        n = length(u)
        diff_vars = findall(i -> mass_matrix[i, i] != 0, 1:n)
        alg_vars = findall(i -> mass_matrix[i, i] == 0, 1:n)
        n_g = length(alg_vars)
        n_f = length(diff_vars)
        g_z = zeros(eltype(u), n_g, n_g)
        g_y = zeros(eltype(u), n_g, n_f)
        W_z = zeros(eltype(u), n_g, n_g)
        linsolve_tmp_z = zeros(eltype(u), n_g)
    end

    return HybridExplicitImplicitCache(
        u, uprev, dense, du, du1, du2, ks,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf,
        linsolve_tmp, linsolve, jac_config, grad_config, reltol, alg,
        alg.step_limiter!, alg.stage_limiter!, interp_order,
        diff_vars, alg_vars, g_z, g_y, W_z, linsolve_tmp_z
    )
end
