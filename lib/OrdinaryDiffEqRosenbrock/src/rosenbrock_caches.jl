abstract type RosenbrockMutableCache <: OrdinaryDiffEqMutableCache end
abstract type GenericRosenbrockMutableCache <: RosenbrockMutableCache end
abstract type RosenbrockConstantCache <: OrdinaryDiffEqConstantCache end

# Fake values since non-FSAL
get_fsalfirstlast(cache::RosenbrockMutableCache, u) = (nothing, nothing)
function get_fsalfirstlast(cache::GenericRosenbrockMutableCache, u)
    (cache.fsalfirst, cache.fsallast)
end

@cache mutable struct Rosenbrock23Cache{uType, rateType, uNoUnitsType, JType, WType,
    TabType, TFType, UFType, F, JCType, GCType,
    RTolType, A, AV, StepLimiter, StageLimiter} <: RosenbrockMutableCache
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
end

@cache mutable struct Rosenbrock32Cache{uType, rateType, uNoUnitsType, JType, WType,
    TabType, TFType, UFType, F, JCType, GCType,
    RTolType, A, AV, StepLimiter, StageLimiter} <: RosenbrockMutableCache
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
end

function alg_cache(alg::Rosenbrock23, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
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
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rosenbrock23Tableau(constvalue(uBottomEltypeNoUnits))
    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)
    linsolve_tmp = zero(rate_prototype)

    linprob = LinearProblem(W, _vec(linsolve_tmp); u0 = _vec(tmp))
    Pl, Pr = wrapprecs(
        alg.precs(W, nothing, u, p, t, nothing, nothing, nothing,
            nothing)..., weight, tmp)
    linsolve = init(linprob, alg.linsolve, alias_A = true, alias_b = true,
        Pl = Pl, Pr = Pr,
        assumptions = LinearSolve.OperatorAssumptions(true))

    grad_config = build_grad_config(alg, f, tf, du1, t)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, du2)
    algebraic_vars = f.mass_matrix === I ? nothing :
                     [all(iszero, x) for x in eachcol(f.mass_matrix)]

    Rosenbrock23Cache(u, uprev, k₁, k₂, k₃, du1, du2, f₁,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf,
        linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, algebraic_vars, alg.step_limiter!,
        alg.stage_limiter!)
end

function alg_cache(alg::Rosenbrock32, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
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
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rosenbrock32Tableau(constvalue(uBottomEltypeNoUnits))

    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)
    linsolve_tmp = zero(rate_prototype)
    linprob = LinearProblem(W, _vec(linsolve_tmp); u0 = _vec(tmp))

    Pl, Pr = wrapprecs(
        alg.precs(W, nothing, u, p, t, nothing, nothing, nothing,
            nothing)..., weight, tmp)
    linsolve = init(linprob, alg.linsolve, alias_A = true, alias_b = true,
        Pl = Pl, Pr = Pr,
        assumptions = LinearSolve.OperatorAssumptions(true))
    grad_config = build_grad_config(alg, f, tf, du1, t)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, du2)
    algebraic_vars = f.mass_matrix === I ? nothing :
                     [all(iszero, x) for x in eachcol(f.mass_matrix)]

    Rosenbrock32Cache(u, uprev, k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W,
        tmp, atmp, weight, tab, tf, uf, linsolve_tmp, linsolve, jac_config,
        grad_config, reltol, alg, algebraic_vars, alg.step_limiter!, alg.stage_limiter!)
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
end

function Rosenbrock23ConstantCache(::Type{T}, tf, uf, J, W, linsolve, autodiff) where {T}
    tab = Rosenbrock23Tableau(T)
    Rosenbrock23ConstantCache(tab.c₃₂, tab.d, tf, uf, J, W, linsolve, autodiff)
end

function alg_cache(alg::Rosenbrock23, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rosenbrock23ConstantCache(constvalue(uBottomEltypeNoUnits), tf, uf, J, W, linsolve,
        alg_autodiff(alg))
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
end

function Rosenbrock32ConstantCache(::Type{T}, tf, uf, J, W, linsolve, autodiff) where {T}
    tab = Rosenbrock32Tableau(T)
    Rosenbrock32ConstantCache(tab.c₃₂, tab.d, tf, uf, J, W, linsolve, autodiff)
end

function alg_cache(alg::Rosenbrock32, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rosenbrock32ConstantCache(constvalue(uBottomEltypeNoUnits), tf, uf, J, W, linsolve,
        alg_autodiff(alg))
end

################################################################################

# Shampine's Low-order Rosenbrocks
mutable struct RosenbrockCache{uType, rateType, uNoUnitsType, JType, WType, TabType,
    TFType, UFType, F, JCType, GCType, RTolType, A, AV, StepLimiter, StageLimiter} <: RosenbrockMutableCache
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
    algebraic_vars::AV
    step_limiter!::StepLimiter
    stage_limiter!::StageLimiter
    interp_order::Int
end

function full_cache(c::RosenbrockCache)
    return [c.u, c.uprev, c.dense..., c.du, c.du1, c.du2,
        c.ks..., c.fsalfirst, c.fsallast, c.dT, c.tmp, c.atmp, c.weight, c.linsolve_tmp]
end

struct RosenbrockCombinedConstantCache{TF, UF, Tab, JType, WType, F, AD} <: RosenbrockConstantCache
    tf::TF
    uf::UF
    tab::Tab
    J::JType
    W::WType
    linsolve::F
    autodiff::AD
    interp_order::Int
end

@ROS2(:cache)

################################################################################

### ROS23 methods

@ROS23(:cache)

################################################################################

### ROS34PW methods

@ROS34PW(:cache)

################################################################################

### ROS4 methods

@Rosenbrock4(:cache)
jac_cache(c::Rosenbrock4Cache) = (c.J, c.W)

###############################################################################

### Rodas methods
tabtype(::Rosenbrock23) = Rosenbrock23Tableau
tabtype(::Rosenbrock32) = Rosenbrock32Tableau
tabtype(::Rodas23W) = Rodas23WTableau
tabtype(::ROS3P) = ROS3PTableau
tabtype(::Rodas3) = Rodas3Tableau
tabtype(::Rodas4) = Rodas4Tableau
tabtype(::Rodas42) = Rodas42Tableau
tabtype(::Rodas4P) = Rodas4PTableau
tabtype(::Rodas4P2) = Rodas4P2Tableau
tabtype(::Rodas5) = Rodas5Tableau
tabtype(::Rodas5P) = Rodas5PTableau
tabtype(::Rodas5Pr) = Rodas5PTableau
tabtype(::Rodas5Pe) = Rodas5PTableau

function alg_cache(
        alg::Union{ROS3P, Rodas3, Rodas23W, Rodas4, Rodas42, Rodas4P, Rodas4P2, Rodas5, Rodas5P, Rodas5Pe, Rodas5Pr},
        u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = tabtype(alg)(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    dense = [zero(rate_prototype) for _ in 1:size(tab.H, 1)]
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    ks = [zero(rate_prototype) for _ in 1:size(tab.A, 1)]

    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)

    tf = TimeGradientWrapper(f, uprev, p)
    uf = UJacobianWrapper(f, t, p)
    linsolve_tmp = zero(rate_prototype)
    linprob = LinearProblem(W, _vec(linsolve_tmp); u0 = _vec(tmp))
    Pl, Pr = wrapprecs(
        alg.precs(W, nothing, u, p, t, nothing, nothing, nothing,
            nothing)..., weight, tmp)
    linsolve = init(linprob, alg.linsolve, alias_A = true, alias_b = true,
        Pl = Pl, Pr = Pr,
        assumptions = LinearSolve.OperatorAssumptions(true))
    grad_config = build_grad_config(alg, f, tf, du1, t)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, du2)
    algebraic_vars = f.mass_matrix === I ? nothing :
                     [all(iszero, x) for x in eachcol(f.mass_matrix)]
    RosenbrockCache(u, uprev, dense, du, du1, du2, ks, fsalfirst, fsallast,
                    dT, J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp,
                    linsolve, jac_config, grad_config, reltol, alg, algebraic_vars,
                    alg.step_limiter!, alg.stage_limiter!, size(tab.H, 1))
end

function alg_cache(
        alg::Union{ROS3P, Rodas3, Rodas23W, Rodas4, Rodas42, Rodas4P, Rodas4P2, Rodas5, Rodas5P, Rodas5Pe, Rodas5Pr},
        u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    tab = tabtype(alg)(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    RosenbrockCombinedConstantCache(tf, uf, tab, J, W, linsolve, alg_autodiff(alg), size(tab.H, 1))
end

function get_fsalfirstlast(
        cache::Union{RosenbrockCache,
            Rosenbrock4Cache},
        u)
    (cache.fsalfirst, cache.fsallast)
end

################################################################################

### RosenbrockW6S4O

@RosenbrockW6S4OS(:cache)
