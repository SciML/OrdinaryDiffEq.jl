abstract type FIRKMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::FIRKMutableCache, u) = (cache.fsalfirst, cache.k)

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
    F1, Tab, Tol, Dt, rTol, aTol, StepLimiter} <: FIRKMutableCache
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
    step_limiter!::StepLimiter
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
    linsolve = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
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
        Convergence, alg.step_limiter!)
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
    UF, JC, F1, F2, Tab, Tol, Dt, rTol, aTol, StepLimiter} <:
               FIRKMutableCache
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
    step_limiter!::StepLimiter
end

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
    linsolve1 = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        assumptions = LinearSolve.OperatorAssumptions(true))
    #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
    #Pr = Diagonal(_vec(weight)))
    linprob = LinearProblem(W2, _vec(cubuff); u0 = _vec(dw23))
    linsolve2 = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
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
        Convergence, alg.step_limiter!)
end

mutable struct RadauIIA9ConstantCache{F, Tab, Tol, Dt, U, JType} <:
               OrdinaryDiffEqConstantCache
    uf::F
    tab::Tab
    κ::Tol
    ηold::Tol
    iter::Int
    cont1::U
    cont2::U
    cont3::U
    cont4::U
    cont5::U
    dtprev::Dt
    W_γdt::Dt
    status::NLStatus
    J::JType
end

function alg_cache(alg::RadauIIA9, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uf = UDerivativeWrapper(f, t, p)
    uToltype = constvalue(uBottomEltypeNoUnits)
    tab = RadauIIA9Tableau(uToltype, constvalue(tTypeNoUnits))

    κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1 // 100)
    J = false .* _vec(rate_prototype) .* _vec(rate_prototype)'

    RadauIIA9ConstantCache(uf, tab, κ, one(uToltype), 10000, u, u, u, u, u, dt, dt,
        Convergence, J)
end

mutable struct RadauIIA9Cache{uType, cuType, uNoUnitsType, rateType, JType, W1Type, W2Type,
    UF, JC, F1, F2, Tab, Tol, Dt, rTol, aTol, StepLimiter} <:
               FIRKMutableCache
    u::uType
    uprev::uType
    z1::uType
    z2::uType
    z3::uType
    z4::uType
    z5::uType
    w1::uType
    w2::uType
    w3::uType
    w4::uType
    w5::uType
    dw1::uType
    ubuff::uType
    dw23::cuType
    dw45::cuType
    cubuff1::cuType
    cubuff2::cuType
    cont1::uType
    cont2::uType
    cont3::uType
    cont4::uType
    cont5::uType
    du1::rateType
    fsalfirst::rateType
    k::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    fw1::rateType
    fw2::rateType
    fw3::rateType
    fw4::rateType
    fw5::rateType
    J::JType
    W1::W1Type
    W2::W2Type # complex
    W3::W2Type
    uf::UF
    tab::Tab
    κ::Tol
    ηold::Tol
    iter::Int
    tmp::uType
    tmp2::uType
    tmp3::uType
    tmp4::uType
    tmp5::uType
    tmp6::uType
    tmp7::uType
    tmp8::uType
    tmp9::uType
    tmp10::uType
    atmp::uNoUnitsType
    jac_config::JC
    linsolve1::F1
    linsolve2::F2
    linsolve3::F2
    rtol::rTol
    atol::aTol
    dtprev::Dt
    W_γdt::Dt
    status::NLStatus
    step_limiter!::StepLimiter
end

function alg_cache(alg::RadauIIA9, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uf = UJacobianWrapper(f, t, p)
    uToltype = constvalue(uBottomEltypeNoUnits)
    tab = RadauIIA9Tableau(uToltype, constvalue(tTypeNoUnits))

    κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1 // 100)

    z1 = zero(u)
    z2 = zero(u)
    z3 = zero(u)
    z4 = zero(u)
    z5 = zero(u)
    w1 = zero(u)
    w2 = zero(u)
    w3 = zero(u)
    w4 = zero(u)
    w5 = zero(u)
    dw1 = zero(u)
    ubuff = zero(u)
    dw23 = similar(u, Complex{eltype(u)})
    dw45 = similar(u, Complex{eltype(u)})
    recursivefill!(dw23, false)
    recursivefill!(dw45, false)
    cubuff1 = similar(u, Complex{eltype(u)})
    cubuff2 = similar(u, Complex{eltype(u)})
    recursivefill!(cubuff1, false)
    recursivefill!(cubuff2, false)
    cont1 = zero(u)
    cont2 = zero(u)
    cont3 = zero(u)
    cont4 = zero(u)
    cont5 = zero(u)

    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    fw1 = zero(rate_prototype)
    fw2 = zero(rate_prototype)
    fw3 = zero(rate_prototype)
    fw4 = zero(rate_prototype)
    fw5 = zero(rate_prototype)

    J, W1 = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    if J isa AbstractSciMLOperator
        error("Non-concrete Jacobian not yet supported by RadauIIA5.")
    end
    W2 = similar(J, Complex{eltype(W1)})
    W3 = similar(J, Complex{eltype(W1)})
    recursivefill!(W2, false)
    recursivefill!(W3, false)

    du1 = zero(rate_prototype)

    tmp = zero(u)
    tmp2 = zero(u)
    tmp3 = zero(u)
    tmp4 = zero(u)
    tmp5 = zero(u)
    tmp6 = zero(u)
    tmp7 = zero(u)
    tmp8 = zero(u)
    tmp9 = zero(u)
    tmp10 = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, tmp, dw1)

    linprob = LinearProblem(W1, _vec(ubuff); u0 = _vec(dw1))
    linsolve1 = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        assumptions = LinearSolve.OperatorAssumptions(true))
    #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
    #Pr = Diagonal(_vec(weight)))
    linprob = LinearProblem(W2, _vec(cubuff1); u0 = _vec(dw23))
    linsolve2 = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        assumptions = LinearSolve.OperatorAssumptions(true))
    #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
    #Pr = Diagonal(_vec(weight)))
    linprob = LinearProblem(W3, _vec(cubuff2); u0 = _vec(dw45))
    linsolve3 = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        assumptions = LinearSolve.OperatorAssumptions(true))
    #Pl = LinearSolve.InvPreconditioner(Diagonal(_vec(weight))),
    #Pr = Diagonal(_vec(weight)))

    rtol = reltol isa Number ? reltol : zero(reltol)
    atol = reltol isa Number ? reltol : zero(reltol)

    RadauIIA9Cache(u, uprev,
        z1, z2, z3, z4, z5, w1, w2, w3, w4, w5,
        dw1, ubuff, dw23, dw45, cubuff1, cubuff2, cont1, cont2, cont3, cont4, cont5,
        du1, fsalfirst, k, k2, k3, k4, k5, fw1, fw2, fw3, fw4, fw5,
        J, W1, W2, W3,
        uf, tab, κ, one(uToltype), 10000,
        tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, atmp, jac_config,
        linsolve1, linsolve2, linsolve3, rtol, atol, dt, dt,
        Convergence, alg.step_limiter!)
end

mutable struct AdaptiveRadauConstantCache{F, Tab, Tol, Dt, U, JType} <:
               OrdinaryDiffEqConstantCache
    uf::F
    tabs::Vector{Tab}
    κ::Tol
    ηold::Tol
    iter::Int
    cont::Vector{U}
    dtprev::Dt
    W_γdt::Dt
    status::NLStatus
    J::JType
    num_stages::Int
    step::Int
    hist_iter::Float64
    index::Int
end

function alg_cache(alg::AdaptiveRadau, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uf = UDerivativeWrapper(f, t, p)
    uToltype = constvalue(uBottomEltypeNoUnits)

    max_order = alg.max_order
    min_order = alg.min_order
    max_stages = (max_order - 1) ÷ 4 * 2 + 1
    min_stages = (min_order - 1) ÷ 4 * 2 + 1
    if (alg.min_order < 5)
        error("min_order choice $min_order below 5 is not compatible with the algorithm")
    elseif (max_stages < min_stages)
        error("max_order $max_order is below min_order $min_order")
    end
    num_stages = min_stages

    tabs = [RadauIIATableau(uToltype, constvalue(tTypeNoUnits), i)
            for i in min_stages:2:max_stages]
    cont = Vector{typeof(u)}(undef, max_stages)
    for i in 1:max_stages
        cont[i] = zero(u)
    end

    index = 1

    κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1 // 100)
    J = false .* _vec(rate_prototype) .* _vec(rate_prototype)'
    AdaptiveRadauConstantCache(uf, tabs, κ, one(uToltype), 10000, cont, dt, dt,
        Convergence, J, num_stages, 1, 0.0, index)
end

mutable struct AdaptiveRadauCache{
    uType, cuType, tType, uNoUnitsType, rateType, JType, W1Type, W2Type,
    UF, JC, F1, F2, Tab, Tol, Dt, rTol, aTol, StepLimiter} <:
               FIRKMutableCache
    u::uType
    uprev::uType
    z::Vector{uType}
    w::Vector{uType}
    c_prime::Vector{tType}
    αdt::Vector{tType}
    βdt::Vector{tType}
    dw1::uType
    ubuff::uType
    dw2::Vector{cuType}
    cubuff::Vector{cuType}
    dw::Vector{uType}
    cont::Vector{uType}
    derivatives::Matrix{uType}
    du1::rateType
    fsalfirst::rateType
    ks::Vector{rateType}
    k::rateType
    fw::Vector{rateType}
    J::JType
    W1::W1Type #real
    W2::Vector{W2Type} #complex
    uf::UF
    tabs::Vector{Tab}
    κ::Tol
    ηold::Tol
    iter::Int
    tmp::uType
    atmp::uNoUnitsType
    jac_config::JC
    linsolve1::F1 #real
    linsolve2::Vector{F2} #complex
    rtol::rTol
    atol::aTol
    dtprev::Dt
    W_γdt::Dt
    status::NLStatus
    step_limiter!::StepLimiter
    num_stages::Int
    step::Int
    hist_iter::Float64
    index::Int
end

function alg_cache(alg::AdaptiveRadau, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uf = UJacobianWrapper(f, t, p)
    uToltype = constvalue(uBottomEltypeNoUnits)

    max_order = alg.max_order
    min_order = alg.min_order
    max_stages = (max_order - 1) ÷ 4 * 2 + 1
    min_stages = (min_order - 1) ÷ 4 * 2 + 1
    if (alg.min_order < 5)
        error("min_order choice $min_order below 5 is not compatible with the algorithm")
    elseif (max_stages < min_stages)
        error("max_order $max_order is below min_order $min_order")
    end
    num_stages = min_stages

    tabs = [RadauIIATableau(uToltype, constvalue(tTypeNoUnits), i)
            for i in min_stages:2:max_stages]

    index = 1

    κ = alg.κ !== nothing ? convert(uToltype, alg.κ) : convert(uToltype, 1 // 100)

    z = Vector{typeof(u)}(undef, max_stages)
    w = Vector{typeof(u)}(undef, max_stages)
    for i in 1:max_stages
        z[i] = zero(u)
        w[i] = zero(u)
    end

    αdt = [zero(t) for i in 1:max_stages]
    βdt = [zero(t) for i in 1:max_stages]
    c_prime = Vector{typeof(t)}(undef, max_stages) #time stepping
    for i in 1:max_stages
        c_prime[i] = zero(t)
    end

    dw1 = zero(u)
    ubuff = zero(u)
    dw2 = [similar(u, Complex{eltype(u)}) for _ in 1:((max_stages - 1) ÷ 2)]
    recursivefill!.(dw2, false)
    cubuff = [similar(u, Complex{eltype(u)}) for _ in 1:((max_stages - 1) ÷ 2)]
    recursivefill!.(cubuff, false)
    dw = [zero(u) for i in 1:max_stages]

    cont = [zero(u) for i in 1:max_stages]

    derivatives = Matrix{typeof(u)}(undef, max_stages, max_stages)
    for i in 1:max_stages, j in 1:max_stages
        derivatives[i, j] = zero(u)
    end

    fsalfirst = zero(rate_prototype)
    fw = [zero(rate_prototype) for i in 1:max_stages]
    ks = [zero(rate_prototype) for i in 1:max_stages]

    k = ks[1]

    J, W1 = build_J_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits, Val(true))
    if J isa AbstractSciMLOperator
        error("Non-concrete Jacobian not yet supported by AdaptiveRadau.")
    end

    W2 = [similar(J, Complex{eltype(W1)}) for _ in 1:((max_stages - 1) ÷ 2)]
    recursivefill!.(W2, false)

    du1 = zero(rate_prototype)

    tmp = zero(u)

    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    jac_config = build_jac_config(alg, f, uf, du1, uprev, u, zero(u), dw1)

    linprob = LinearProblem(W1, _vec(ubuff); u0 = _vec(dw1))
    linsolve1 = init(
        linprob, alg.linsolve, alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        assumptions = LinearSolve.OperatorAssumptions(true))

    linsolve2 = [init(LinearProblem(W2[i], _vec(cubuff[i]); u0 = _vec(dw2[i])),
                     alg.linsolve, alias = LinearAliasSpecifier(
                         alias_A = true, alias_b = true),
                     assumptions = LinearSolve.OperatorAssumptions(true))
                 for i in 1:((max_stages - 1) ÷ 2)]

    rtol = reltol isa Number ? reltol : zero(reltol)
    atol = reltol isa Number ? reltol : zero(reltol)

    AdaptiveRadauCache(u, uprev,
        z, w, c_prime, αdt, βdt, dw1, ubuff, dw2, cubuff, dw, cont, derivatives,
        du1, fsalfirst, ks, k, fw,
        J, W1, W2,
        uf, tabs, κ, one(uToltype), 10000, tmp,
        atmp, jac_config,
        linsolve1, linsolve2, rtol, atol, dt, dt,
        Convergence, alg.step_limiter!, num_stages, 1, 0.0, index)
end
