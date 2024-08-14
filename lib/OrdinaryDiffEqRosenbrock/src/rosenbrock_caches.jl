abstract type RosenbrockMutableCache <: OrdinaryDiffEqMutableCache end
################################################################################

# Shampine's Low-order Rosenbrocks

@cache mutable struct RosenbrockConstantCache{T, TF, UF, Tab, JType, WType, F, AD} <: OrdinaryDiffEqConstantCache
    c₃₂::T
    d::T
    tf::TF
    uf::UF
    J::JType
    W::WType
    linsolve::F
    autodiff::AD
    tab::Tab
end

@cache mutable struct RosenbrockCache{uType, rateType, uNoUnitsType, JType, WType,
    TabType, TFType, UFType, F, JCType, GCType,
    RTolType, A, AV, AD, TF, UF, T, StepLimiter, StageLimiter} <: RosenbrockMutableCache
    u::uType
    uprev::uType
    dus::Vector{rateType}
    ks::Array{rateType, 1}
    f1::rateType
    dense::Array{rateType, 1}
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
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    ks = [k1, k2, k3]
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du1, du2]
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

    grad_config = build_grad_config(alg, f, tf, dus[1], t)
    jac_config = build_jac_config(alg, f, uf, dus[1], uprev, u, tmp, du2, Val(false))
    algebraic_vars = f.mass_matrix === I ? nothing :
                     [all(iszero, x) for x in eachcol(f.mass_matrix)]

    RosenbrockCache(u, uprev, ks[1], ks[2], ks[3], dus[1], dus[2], f₁,
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf,
        linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, algebraic_vars, alg.step_limiter!,
        alg.stage_limiter!)
end

function alg_cache(alg::Rosenbrock32, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    ks = [k1, k2, k3]
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du1, du2]
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

    RosenbrockCache(u, uprev, ks[1], ks[2], ks[3], dus[1], dus[2], f₁, fsalfirst, fsallast, dT, J, W,
        tmp, atmp, weight, tab, tf, uf, linsolve_tmp, linsolve, jac_config,
        grad_config, reltol, alg, algebraic_vars, alg.step_limiter!, alg.stage_limiter!)
end

function RosenbrockConstantCache(::Type{T}, tf, uf, J, W, linsolve, autodiff) where {T}
    tab = Rosenbrock23Tableau(T)
    RosenbrockConstantCache(tab.c₃₂, tab.d, tf, uf, J, W, linsolve, autodiff)
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
    RosenbrockConstantCache(constvalue(uBottomEltypeNoUnits), tf, uf, J, W, linsolve,
        alg_autodiff(alg))
end

function Rosenbrock32ConstantCache(::Type{T}, tf, uf, J, W, linsolve, autodiff) where {T}
    tab = Rosenbrock32Tableau(T)
    RosenbrockConstantCache(tab.c₃₂, tab.d, tf, uf, J, W, linsolve, autodiff)
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
    RosenbrockConstantCache(constvalue(uBottomEltypeNoUnits), tf, uf, J, W, linsolve,
        alg_autodiff(alg))
end

################################################################################

### 3rd order specialized Rosenbrocks

function alg_cache(alg::ROS3P, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du, du1, du2]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    ks = [k1, k2, k3, k4]
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
    RosenbrockCache(u, uprev, dus[1], dus[2], dus[3], ks[1], ks[2], ks[3], ks[4],
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
    RosenbrockConstantCache(tf, uf,
        ROS3PTableau(constvalue(uBottomEltypeNoUnits),
            constvalue(tTypeNoUnits)), J, W, linsolve)
end

function alg_cache(alg::Rodas3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du, du1, du2]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    ks = [k1, k2, k3, k4]
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
    RosenbrockCache(u, uprev, dus[1], dus[2], dus[3], ks[1], ks[2], ks[3], ks[4],
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf,
        linsolve_tmp,
        linsolve, jac_config, grad_config, alg.step_limiter!,
        alg.stage_limiter!)
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
    RosenbrockConstantCache(tf, uf,
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

function alg_cache(alg::Rodas23W, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    dense3 = zero(rate_prototype)
    dense = [dense1, dense2, dense3]
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du, du1, du2]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    ks = [k1, k2, k3, k4, k5]
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
    Rodas23WCache(u, uprev, dense[1], dense[2], dense[3], dus[1], dus[2], dus[3], ks[1], ks[2], ks[3], ks[4], ks[5],
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

function alg_cache(alg::Rodas3P, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    dense3 = zero(rate_prototype)
    dense = [dense1, dense2, dense3]
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du, du1, du2]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    ks = [k1, k2, k3, k4, k5]
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
    Rodas3PCache(u, uprev, dense[1], dense[2], dense[3], dus[1], dus[2], dus[3], ks[1], ks[2], ks[3], ks[4], ks[5],
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

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

function alg_cache(alg::Rodas4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    dense = [dense1, dense2]
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du, du1, du2]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    ks = [k1, k2, k3, k4, k5, k6]
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
    Rodas4Cache(u, uprev, dense[1], dense[2], dus[1], dus[2], dus[3], ks[1], ks[2], ks[3], ks[4],
        ks[5], ks[6],
        fsalfirst, fsallast, dT, J, W, tmp, atmp, weight, tab, tf, uf, linsolve_tmp,
        linsolve, jac_config, grad_config, reltol, alg, alg.step_limiter!,
        alg.stage_limiter!)
end

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
    dense = [dense1, dense2]
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du, du1, du2]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    ks = [k1, k2, k3, k4, k5, k6]
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
    Rodas4Cache(u, uprev, dense[1], dense[2], dus[1], dus[2], dus[3], ks[1], ks[2], ks[3], ks[4],
        ks[5], ks[6],
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
    dense = [dense1, dense2]
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du, du1, du2]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    ks = [k1, k2, k3, k4, k5, k6]
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
    Rodas4Cache(u, uprev, dense[1], dense[2], dus[1], dus[2], dus[3], ks[1], ks[2], ks[3], ks[4],
        ks[5], ks[6],
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
    dense = [dense1, dense2]
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du, du1, du2]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    ks = [k1, k2, k3, k4, k5, k6]
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
    grad_config = build_grad_config(alg, f, tf, dus[1], t)
    jac_config = build_jac_config(alg, f, uf, dus[2], uprev, u, tmp, dus[3])
    Rodas4Cache(u, uprev, dense[1], dense[2], dus[1], dus[2], dus[3], ks[1], ks[2], ks[3], ks[4],
        ks[5], ks[6],
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
function alg_cache(alg::Rodas5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dense1 = zero(rate_prototype)
    dense2 = zero(rate_prototype)
    dense3 = zero(rate_prototype)
    dense = [dense1, dense2, dense3]
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du, du1, du2]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    ks = [k1, k2, k3, k4, k5, k6, k7, k8]
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
    Rosenbrock5Cache(u, uprev, dense[1], dense[2], dense[3], dus[1], dus[2], dus[3], ks[1], ks[2], ks[3], ks[4],
        ks[5], ks[6], ks[7], ks[8],
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
    dense = [dense1, dense2, dense3]
    du = zero(rate_prototype)
    du1 = zero(rate_prototype)
    du2 = zero(rate_prototype)
    dus = [du, du1, du2]
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    k8 = zero(rate_prototype)
    ks = [k1, k2, k3, k4, k5, k6, k7, k8]
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
    Rosenbrock5Cache(u, uprev, dense[1], dense[2], dense[3], dus[1], dus[2], dus[3], ks[1], ks[2], ks[3], ks[4],
        ks[5], ks[6], ks[7], ks[8],
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