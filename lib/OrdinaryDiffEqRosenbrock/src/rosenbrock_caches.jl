abstract type RosenbrockMutableCache <: OrdinaryDiffEqMutableCache end
################################################################################

# Shampine's Low-order Rosenbrocks

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

TruncatedStacktraces.@truncate_stacktrace Rosenbrock23Cache 1

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
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, du2, Val(false))
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
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, du2, Val(false))
    algebraic_vars = f.mass_matrix === I ? nothing :
                     [all(iszero, x) for x in eachcol(f.mass_matrix)]

    Rosenbrock32Cache(u, uprev, k₁, k₂, k₃, du1, du2, f₁, fsalfirst, fsallast, dT, J, W,
        tmp, atmp, weight, tab, tf, uf, linsolve_tmp, linsolve, jac_config,
        grad_config, reltol, alg, algebraic_vars, alg.step_limiter!, alg.stage_limiter!)
end

struct Rosenbrock23ConstantCache{T, TF, UF, JType, WType, F, AD} <:
       OrdinaryDiffEqConstantCache
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
       OrdinaryDiffEqConstantCache
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

### 3rd order specialized Rosenbrocks

struct Rosenbrock33ConstantCache{TF, UF, Tab, JType, WType, F} <:
       OrdinaryDiffEqConstantCache
    tf::TF
    uf::UF
    tab::Tab
    J::JType
    W::WType
    linsolve::F
end

@cache mutable struct Rosenbrock33Cache{uType, rateType, uNoUnitsType, JType, WType,
    TabType, TFType, UFType, F, JCType, GCType,
    RTolType, A, StepLimiter, StageLimiter} <: RosenbrockMutableCache
    u::uType
    uprev::uType
    du::rateType
    du1::rateType
    du2::rateType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
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
end

function alg_cache(alg::ROS3P, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = ROS3PTableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
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
    Rosenbrock33Cache(u, uprev, du, du1, du2, k1, k2, k3, k4,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf,
        linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

function alg_cache(alg::ROS3P, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rosenbrock33ConstantCache(tf, uf,
        ROS3PTableau(constvalue(uBottomEltypeNoUnits),
            constvalue(tTypeNoUnits)), J, W, linsolve)
end

@cache mutable struct Rosenbrock34Cache{uType, rateType, uNoUnitsType, JType, WType,
    TabType, TFType, UFType, F, JCType, GCType, StepLimiter, StageLimiter} <:
                      RosenbrockMutableCache
    u::uType
    uprev::uType
    du::rateType
    du1::rateType
    du2::rateType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
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
    step_limiter!::StepLimiter
    stage_limiter!::StageLimiter
end

function alg_cache(alg::Rodas3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rodas3Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))

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
    Rosenbrock34Cache(u, uprev, du, du1, du2, k1, k2, k3, k4,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf,
        linsolve_tmp,
        linsolve, jac_config, grad_config, alg.step_limiter!,
        alg.stage_limiter!)
end

struct Rosenbrock34ConstantCache{TF, UF, Tab, JType, WType, F} <:
       OrdinaryDiffEqConstantCache
    tf::TF
    uf::UF
    tab::Tab
    J::JType
    W::WType
    linsolve::F
end

function alg_cache(alg::Rodas3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rosenbrock34ConstantCache(tf, uf,
        Rodas3Tableau(constvalue(uBottomEltypeNoUnits),
            constvalue(tTypeNoUnits)), J, W, linsolve)
end

################################################################################

### ROS2 methods

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

struct Rodas23WConstantCache{TF, UF, Tab, JType, WType, F, AD} <:
       OrdinaryDiffEqConstantCache
    tf::TF
    uf::UF
    tab::Tab
    J::JType
    W::WType
    linsolve::F
    autodiff::AD
end

struct Rodas3PConstantCache{TF, UF, Tab, JType, WType, F, AD} <: OrdinaryDiffEqConstantCache
    tf::TF
    uf::UF
    tab::Tab
    J::JType
    W::WType
    linsolve::F
    autodiff::AD
end

@cache mutable struct Rodas23WCache{uType, rateType, uNoUnitsType, JType, WType, TabType,
    TFType, UFType, F, JCType, GCType, RTolType, A, StepLimiter, StageLimiter} <:
                      RosenbrockMutableCache
    u::uType
    uprev::uType
    dense1::rateType
    dense2::rateType
    dense3::rateType
    du::rateType
    du1::rateType
    du2::rateType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
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
end

@cache mutable struct Rodas3PCache{uType, rateType, uNoUnitsType, JType, WType, TabType,
    TFType, UFType, F, JCType, GCType, RTolType, A, StepLimiter, StageLimiter} <:
                      RosenbrockMutableCache
    u::uType
    uprev::uType
    dense1::rateType
    dense2::rateType
    dense3::rateType
    du::rateType
    du1::rateType
    du2::rateType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
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
end

function alg_cache(alg::Rodas23W, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    dense3 = zero(rate_prototype)
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rodas3PTableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))

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
    Rodas23WCache(u, uprev, dense1, dense2, dense3, du, du1, du2, k1, k2, k3, k4, k5,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

TruncatedStacktraces.@truncate_stacktrace Rodas23WCache 1
function alg_cache(alg::Rodas3P, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    dense3 = zero(rate_prototype)
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rodas3PTableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))

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
    Rodas3PCache(u, uprev, dense1, dense2, dense3, du, du1, du2, k1, k2, k3, k4, k5,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

TruncatedStacktraces.@truncate_stacktrace Rodas3PCache 1

function alg_cache(alg::Rodas23W, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rodas23WConstantCache(tf, uf,
        Rodas3PTableau(constvalue(uBottomEltypeNoUnits),
            constvalue(tTypeNoUnits)), J, W, linsolve,
        alg_autodiff(alg))
end

function alg_cache(alg::Rodas3P, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rodas3PConstantCache(tf, uf,
        Rodas3PTableau(constvalue(uBottomEltypeNoUnits),
            constvalue(tTypeNoUnits)), J, W, linsolve,
        alg_autodiff(alg))
end

### Rodas4 methods

struct Rodas4ConstantCache{TF, UF, Tab, JType, WType, F, AD} <: OrdinaryDiffEqConstantCache
    tf::TF
    uf::UF
    tab::Tab
    J::JType
    W::WType
    linsolve::F
    autodiff::AD
end

@cache mutable struct Rodas4Cache{uType, rateType, uNoUnitsType, JType, WType, TabType,
    TFType, UFType, F, JCType, GCType, RTolType, A, StepLimiter, StageLimiter} <:
                      RosenbrockMutableCache
    u::uType
    uprev::uType
    dense1::rateType
    dense2::rateType
    du::rateType
    du1::rateType
    du2::rateType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
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
end

function alg_cache(alg::Rodas4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rodas4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))

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
    Rodas4Cache(u, uprev, dense1, dense2, du, du1, du2, k1, k2, k3, k4,
        k5, k6,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

TruncatedStacktraces.@truncate_stacktrace Rodas4Cache 1

function alg_cache(alg::Rodas4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rodas4ConstantCache(tf, uf,
        Rodas4Tableau(constvalue(uBottomEltypeNoUnits),
            constvalue(tTypeNoUnits)), J, W, linsolve,
        alg_autodiff(alg))
end

function alg_cache(alg::Rodas42, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rodas42Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))

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
    Rodas4Cache(u, uprev, dense1, dense2, du, du1, du2, k1, k2, k3, k4,
        k5, k6,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

function alg_cache(alg::Rodas42, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rodas4ConstantCache(tf, uf,
        Rodas42Tableau(constvalue(uBottomEltypeNoUnits),
            constvalue(tTypeNoUnits)), J, W, linsolve,
        alg_autodiff(alg))
end

function alg_cache(alg::Rodas4P, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rodas4PTableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))

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
    Rodas4Cache(u, uprev, dense1, dense2, du, du1, du2, k1, k2, k3, k4,
        k5, k6,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

function alg_cache(alg::Rodas4P, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rodas4ConstantCache(tf, uf,
        Rodas4PTableau(constvalue(uBottomEltypeNoUnits),
            constvalue(tTypeNoUnits)), J, W, linsolve,
        alg_autodiff(alg))
end

function alg_cache(alg::Rodas4P2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rodas4P2Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))

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
    Rodas4Cache(u, uprev, dense1, dense2, du, du1, du2, k1, k2, k3, k4,
        k5, k6,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

function alg_cache(alg::Rodas4P2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rodas4ConstantCache(tf, uf,
        Rodas4P2Tableau(constvalue(uBottomEltypeNoUnits),
            constvalue(tTypeNoUnits)), J, W, linsolve,
        alg_autodiff(alg))
end

################################################################################

### Rosenbrock5

struct Rosenbrock5ConstantCache{TF, UF, Tab, JType, WType, F} <: OrdinaryDiffEqConstantCache
    tf::TF
    uf::UF
    tab::Tab
    J::JType
    W::WType
    linsolve::F
end

@cache mutable struct Rosenbrock5Cache{
    uType, rateType, uNoUnitsType, JType, WType, TabType,
    TFType, UFType, F, JCType, GCType, RTolType, A, StepLimiter, StageLimiter} <:
                      RosenbrockMutableCache
    u::uType
    uprev::uType
    dense1::rateType
    dense2::rateType
    dense3::rateType
    du::rateType
    du1::rateType
    du2::rateType
    k1::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k6::rateType
    k7::rateType
    k8::rateType
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
end

TruncatedStacktraces.@truncate_stacktrace Rosenbrock5Cache 1

function alg_cache(alg::Rodas5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    dense3 = zero(rate_prototype)
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rodas5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))

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
    Rosenbrock5Cache(u, uprev, dense1, dense2, dense3, du, du1, du2, k1, k2, k3, k4,
        k5, k6, k7, k8,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf,
        linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

function alg_cache(alg::Rodas5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rosenbrock5ConstantCache(tf, uf,
        Rodas5Tableau(constvalue(uBottomEltypeNoUnits),
            constvalue(tTypeNoUnits)), J, W, linsolve)
end

function alg_cache(
        alg::Union{Rodas5P, Rodas5Pe, Rodas5Pr}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    dense3 = zero(rate_prototype)
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    fsallast = zero(rate_prototype)
    dT = zero(rate_prototype)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    tmp = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    tab = Rodas5PTableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))

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
    Rosenbrock5Cache(u, uprev, dense1, dense2, dense3, du, du1, du2, k1, k2, k3, k4,
        k5, k6, k7, k8,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf,
        linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

function alg_cache(
        alg::Union{Rodas5P, Rodas5Pe, Rodas5Pr}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tf = TimeDerivativeWrapper(f, u, p)
    uf = UDerivativeWrapper(f, t, p)
    J, W = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(false))
    linprob = nothing #LinearProblem(W,copy(u); u0=copy(u))
    linsolve = nothing #init(linprob,alg.linsolve,alias_A=true,alias_b=true)
    Rosenbrock5ConstantCache(tf, uf,
        Rodas5PTableau(constvalue(uBottomEltypeNoUnits),
            constvalue(tTypeNoUnits)), J, W, linsolve)
end

################################################################################

### RosenbrockW6S4O

@RosenbrockW6S4OS(:cache)