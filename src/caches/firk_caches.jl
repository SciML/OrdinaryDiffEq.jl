mutable struct RadauIIA3ConstantCache{F, Tab, Tol, Dt, U, JType} <:
               OrdinaryDiffEqConstantCache
    uf::F
    tab::Tab
    κ::Tol
    ηold::Tol
    iter::Int
    cont1::U
    cont2::U
    cont3::U
    dtprev::Dt
    W_γdt::Dt
    status::NLStatus
    J::JType
end

function alg_cache(alg::RadauIIA3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uf = UDerivativeWrapper(f, t, p)
    uToltype = constvalue(uBottomEltypeNoUnits)
    tab = RadauIIA3Tableau(uToltype, constvalue(tTypeNoUnits))

    κ = convert(uToltype, 1 // 100)
    J = false .* _vec(rate_prototype) .* _vec(rate_prototype)'

    RadauIIA3ConstantCache(uf, tab, κ, one(uToltype), 10000, u, u, u, dt, dt,
        Convergence, J)
end

mutable struct RadauIIA3Cache{uType, cuType, uNoUnitsType, rateType, JType, W1Type, UF, JC,
    F1, Tab, Tol, Dt, rTol, aTol} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    z1::uType
    z2::uType
    w1::uType
    w2::uType
    dw12::cuType
    cubuff::cuType
    cont1::uType
    cont2::uType
    du1::rateType
    fsalfirst::rateType
    k::rateType
    k2::rateType
    fw1::rateType
    fw2::rateType
    J::JType
    W1::W1Type
    uf::UF
    tab::Tab
    κ::Tol
    ηold::Tol
    iter::Int
    tmp::uType
    atmp::uNoUnitsType
    jac_config::JC
    linsolve::F1
    rtol::rTol
    atol::aTol
    dtprev::Dt
    W_γdt::Dt
    status::NLStatus
end

function alg_cache(alg::RadauIIA3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uf = UJacobianWrapper(f, t, p)
    uToltype = constvalue(uBottomEltypeNoUnits)
    tab = RadauIIA3Tableau(uToltype, constvalue(tTypeNoUnits))

    κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1 // 100)

    z1 = zero(u)
    z2 = zero(u)
    w1 = zero(u)
    w2 = zero(u)
    dw12 = similar(u, Complex{eltype(u)})
    recursivefill!(dw12, false)
    cubuff = similar(u, Complex{eltype(u)})
    recursivefill!(cubuff, false)
    cont1 = zero(u)
    cont2 = zero(u)

    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    k2 = zero(rate_prototype)
    fw1 = zero(rate_prototype)
    fw2 = zero(rate_prototype)

    J, W1 = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    W1 = similar(J, Complex{eltype(W1)})
    recursivefill!(W1, false)

    du1 = zero(rate_prototype)

    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    jac_config = jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dw12)

    linprob = LinearProblem(W1, _vec(cubuff); u0 = _vec(dw12))
    linsolve = init(linprob, alg.linsolve, alias_A = true, alias_b = true,
        assumptions = LinearSolve.OperatorAssumptions(true))
    #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
    #Pr = Diagonal(_vec(weight)))

    rtol = reltol isa Number ? reltol : zero(reltol)
    atol = reltol isa Number ? reltol : zero(reltol)

    RadauIIA3Cache(u, uprev,
        z1, z2, w1, w2,
        dw12, cubuff, cont1, cont2,
        du1, fsalfirst, k, k2, fw1, fw2,
        J, W1,
        uf, tab, κ, one(uToltype), 10000,
        tmp, atmp, jac_config, linsolve, rtol, atol, dt, dt,
        Convergence)
end

mutable struct RadauIIA5ConstantCache{F, Tab, Tol, Dt, U, JType} <:
               OrdinaryDiffEqConstantCache
    uf::F
    tab::Tab
    κ::Tol
    ηold::Tol
    iter::Int
    cont1::U
    cont2::U
    cont3::U
    dtprev::Dt
    W_γdt::Dt
    status::NLStatus
    J::JType
end

function alg_cache(alg::RadauIIA5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uf = UDerivativeWrapper(f, t, p)
    uToltype = constvalue(uBottomEltypeNoUnits)
    tab = RadauIIA5Tableau(uToltype, constvalue(tTypeNoUnits))

    κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1 // 100)
    J = false .* _vec(rate_prototype) .* _vec(rate_prototype)'

    RadauIIA5ConstantCache(uf, tab, κ, one(uToltype), 10000, u, u, u, dt, dt,
        Convergence, J)
end

mutable struct RadauIIA5Cache{uType, cuType, uNoUnitsType, rateType, JType, W1Type, W2Type,
    UF, JC, F1, F2, Tab, Tol, Dt, rTol, aTol} <:
               OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    z1::uType
    z2::uType
    z3::uType
    w1::uType
    w2::uType
    w3::uType
    dw1::uType
    ubuff::uType
    dw23::cuType
    cubuff::cuType
    cont1::uType
    cont2::uType
    cont3::uType
    du1::rateType
    fsalfirst::rateType
    k::rateType
    k2::rateType
    k3::rateType
    fw1::rateType
    fw2::rateType
    fw3::rateType
    J::JType
    W1::W1Type
    W2::W2Type # complex
    uf::UF
    tab::Tab
    κ::Tol
    ηold::Tol
    iter::Int
    tmp::uType
    atmp::uNoUnitsType
    jac_config::JC
    linsolve1::F1
    linsolve2::F2
    rtol::rTol
    atol::aTol
    dtprev::Dt
    W_γdt::Dt
    status::NLStatus
end
TruncatedStacktraces.@truncate_stacktrace RadauIIA5Cache 1

function alg_cache(alg::RadauIIA5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uf = UJacobianWrapper(f, t, p)
    uToltype = constvalue(uBottomEltypeNoUnits)
    tab = RadauIIA5Tableau(uToltype, constvalue(tTypeNoUnits))

    κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1 // 100)

    z1 = zero(u)
    z2 = zero(u)
    z3 = zero(u)
    w1 = zero(u)
    w2 = zero(u)
    w3 = zero(u)
    dw1 = zero(u)
    ubuff = zero(u)
    dw23 = similar(u, Complex{eltype(u)})
    recursivefill!(dw23, false)
    cubuff = similar(u, Complex{eltype(u)})
    recursivefill!(cubuff, false)
    cont1 = zero(u)
    cont2 = zero(u)
    cont3 = zero(u)

    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    fw1 = zero(rate_prototype)
    fw2 = zero(rate_prototype)
    fw3 = zero(rate_prototype)

    J, W1 = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    if J isa AbstractSciMLOperator
        error("Non-concrete Jacobian not yet supported by RadauIIA5.")
    end
    W2 = similar(J, Complex{eltype(W1)})
    recursivefill!(W2, false)

    du1 = zero(rate_prototype)

    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dw1)

    linprob = LinearProblem(W1, _vec(ubuff); u0 = _vec(dw1))
    linsolve1 = init(linprob, alg.linsolve, alias_A = true, alias_b = true,
        assumptions = LinearSolve.OperatorAssumptions(true))
    #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
    #Pr = Diagonal(_vec(weight)))
    linprob = LinearProblem(W2, _vec(cubuff); u0 = _vec(dw23))
    linsolve2 = init(linprob, alg.linsolve, alias_A = true, alias_b = true,
        assumptions = LinearSolve.OperatorAssumptions(true))
    #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
    #Pr = Diagonal(_vec(weight)))

    rtol = reltol isa Number ? reltol : zero(reltol)
    atol = reltol isa Number ? reltol : zero(reltol)

    RadauIIA5Cache(u, uprev,
        z1, z2, z3, w1, w2, w3,
        dw1, ubuff, dw23, cubuff, cont1, cont2, cont3,
        du1, fsalfirst, k, k2, k3, fw1, fw2, fw3,
        J, W1, W2,
        uf, tab, κ, one(uToltype), 10000,
        tmp, atmp, jac_config, linsolve1, linsolve2, rtol, atol, dt, dt,
        Convergence)
end
