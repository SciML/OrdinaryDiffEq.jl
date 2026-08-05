struct MREILConstantCache{T, UF, JFType} <: OrdinaryDiffEqConstantCache
    T::T  # pre-allocated extrapolation table: Vector of length `order`
    uf::UF
    jac_f::JFType
end

@cache mutable struct MREILCache{
        uType, rateType, uNoUnitsType, JType, WType, UFType, JCType, F, JFType,
    } <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    atmp::uNoUnitsType
    weight::uNoUnitsType
    k_slow::rateType
    k_fast::rateType
    du1::rateType
    linsolve_tmp::rateType
    T::Array{uType, 1}
    J::JType
    W::WType
    uf::UFType
    jac_config::JCType
    linsolve::F
    fsalfirst::rateType
    k::rateType
    jac_f::JFType         # the resolved fast component J and W were built from
    J_iter::Int           # `integrator.iter` when J was taken …
    J_success_iter::Int   # … and `success_iter`, so a rejected attempt can reuse J
    J_isvalid::Bool
end

get_fsalfirstlast(cache::MREILCache, u) = (cache.fsalfirst, cache.k)

# `build_J_W` reads `jac_prototype`/`sparsity`, decides whether the right-hand side is
# a linear operator, and builds any matrix-free Jacobian-vector product against the
# function it is handed. All of those must describe the *fast* component for MREIL, so
# hand it `f1` — otherwise a matrix-free `linsolve` silently differentiates `f1 + f2`
# and solves a different method than the factorizing one. The mass matrix lives on the
# outer function, as do `jac`/`jac_prototype`/`sparsity` when the user set them there
# rather than on `f1`.
#
# A `SplitFunction`'s `jac` is `df1/du`: SciMLBase documents it that way ("the
# derivatives … are only defined on the `f1` portion") and its constructor defaults the
# field *from* `f1`, so merging the outer one onto `f1` is the documented reading rather
# than a guess. Resolving it here is what makes reading it safe later: `jac_f` becomes
# the single description of the fast component, `has_jac(jac_f)` agrees with
# `has_jac(f)` by construction, and the shared `calc_J!`/`calc_J` fall-through in the
# step is then reachable only when there is no user Jacobian anywhere.
#
# Leaving `jac` unresolved instead sent that fall-through through `calc_J!`, which pairs
# the outer `jac` with the outer `jac_prototype` while `J` was built from the merged
# one. Those two disagree as soon as a user writes them separately, and `calc_J!` then
# skips the sparse-structure reset a `jac` that only accumulates into stored entries
# needs. With a dense prototype `J[1, 2]` came out 11, 21, 31, 41, 51 over five steps
# of a problem whose fast Jacobian is constant at 10; a sparse prototype resets
# correctly at every placement.
function _mreil_jac_function(f::SplitFunction, nf)
    return SciMLBase.remake(
        nf;
        mass_matrix = f.mass_matrix,
        jac = something(nf.jac, f.jac, Some(nothing)),
        jac_prototype = something(nf.jac_prototype, f.jac_prototype, Some(nothing)),
        sparsity = something(nf.sparsity, f.sparsity, Some(nothing)),
    )
end
_mreil_jac_function(f, nf) = f

# A sparse matrix states its sparsity pattern through its stored structure, so a sparse
# `jac_prototype` whose stored values are all zero still names a pattern. A dense matrix
# has no structure to read, so the pattern has to come from the values, and an all-zero
# dense one therefore says the Jacobian has no nonzeros anywhere:
# `prepare_user_sparsity` builds a `KnownJacobianSparsityDetector` over it, the coloring
# comes back empty, and AD fills in nothing — `J` stays identically zero.
#
# What that costs is *stability*, not order. Extrapolated linearly implicit Euler is a
# W-method: the Richardson order survives an arbitrary frozen `J`, measured at 2.00,
# 2.97 and 4.24 for `order` 2, 3 and 4 with `J = 0` against 2.01, 3.00 and 3.99 with the
# true one. With `J = 0` the substep `(M - h*J)Δ = h*(f1 + f2)` is explicit Euler and
# MREIL reproduces `MREEF`, so the user silently gets the explicit method they
# were avoiding, together with its step size restriction — on a fast block with
# `λ = 500` at `dt = 0.05` that is an error of 6.8e60 where the true `J` gives 7.1e-7.
# `zeros(n, n)` is a natural way to write "a dense n-by-n buffer", so refuse it rather
# than hand back a different method under MREIL's name.
#
# This looks at the outer `f`, because `prepare_user_sparsity` reads `prob.f.sparsity`;
# `ODEFunction` and `SplitFunction` both default `sparsity` from `jac_prototype`. A user
# `jac` fills the prototype itself and never consults the pattern, so it is exempt.
function _mreil_check_ad_sparsity(f)
    SciMLBase.has_jac(f) && return nothing
    hasproperty(f, :sparsity) || return nothing
    sparsity = f.sparsity
    (sparsity isa AbstractMatrix && !issparse(sparsity)) || return nothing
    all(iszero, sparsity) && throw(
        ArgumentError(
            "MREIL: the fast component's `jac_prototype`/`sparsity` is a dense " *
                "all-zero matrix, which declares that the Jacobian has no nonzero " *
                "entries; automatic differentiation would then leave it identically " *
                "zero and MREIL would silently degrade to explicit Euler substeps. " *
                "Drop the prototype to get a dense AD Jacobian, or pass one whose " *
                "nonzero entries are the actual sparsity pattern (`ones(n, n)` for a " *
                "dense Jacobian, a `SparseMatrixCSC` with the right structure " *
                "otherwise), or supply `jac`."
        )
    )
    return nothing
end

function alg_cache(
        alg::MREIL, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    _mreil_check_ad_sparsity(f)
    tmp = zero(u)
    dz = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    weight = similar(u, uEltypeNoUnits)
    recursivefill!(weight, false)
    k_slow = zero(rate_prototype)
    k_fast = zero(rate_prototype)
    du1 = zero(rate_prototype)
    linsolve_tmp = zero(rate_prototype)
    T = [zero(u) for _ in 1:(alg.order)]
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)

    nf = nlsolve_f(f, alg)
    jac_f = _mreil_jac_function(f, nf)
    uf = build_uf(alg, nf, t, p, Val(true))
    jac_config = build_jac_config(alg, nf, uf, du1, uprev, u, tmp, dz)
    J, W = build_J_W(
        alg, u, uprev, p, t, dt, jac_f, jac_config, uEltypeNoUnits, Val(true)
    )

    linprob = LinearProblem(W, _vec(linsolve_tmp), (nothing, u, p, t); u0 = _vec(dz))
    linsolve = init(
        linprob, wrapprecs(alg.linsolve, W, weight),
        alias = LinearAliasSpecifier(alias_A = true, alias_b = true),
        abstol = reltol, reltol = reltol,
        assumptions = LinearSolve.OperatorAssumptions(true),
        verbose = verbose.linear_verbosity
    )

    return MREILCache(
        u, uprev, tmp, dz, atmp, weight, k_slow, k_fast, du1, linsolve_tmp,
        T, J, W, uf, jac_config, linsolve, fsalfirst, k, jac_f, 0, 0, false
    )
end

function alg_cache(
        alg::MREIL, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    _mreil_check_ad_sparsity(f)
    T = Vector{typeof(u)}(undef, alg.order)
    nf = nlsolve_f(f, alg)
    uf = build_uf(alg, nf, t, p, Val(false))
    return MREILConstantCache(T, uf, _mreil_jac_function(f, nf))
end

mutable struct MRIGARKImplicitConstantCache{N, TabType} <: OrdinaryDiffEqConstantCache
    nlsolver::N
    tab::TabType
end

@cache mutable struct MRIGARKImplicitCache{uType, rateType, uNoUnitsType, N, TabType} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    atmp::uNoUnitsType
    v::uType
    vtmp::uType
    f1eval::rateType
    kk::Vector{uType}
    z::Vector{uType}
    fS::Vector{rateType}
    zemb::uType
    fsalfirst::rateType
    k::rateType
    nlsolver::N
    tab::TabType
end

get_fsalfirstlast(cache::MRIGARKImplicitCache, u) = (cache.fsalfirst, cache.k)

_mrigark_impl_γ(tab) = tab.γ0[findfirst(!iszero, tab.γ0)]

function alg_cache(
        alg::MRIGARKImplicitAlg, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = mri_gark_tableau(alg, eltype(u))
    s = length(tab.Δc)
    γ = _mrigark_impl_γ(tab)
    c = one(tTypeNoUnits)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    v = zero(u)
    vtmp = zero(u)
    f1eval = zero(rate_prototype)
    kk = [zero(u) for _ in 1:(tab.q)]
    z = [zero(u) for _ in 1:(s + 1)]
    fS = [zero(rate_prototype) for _ in 1:s]
    zemb = zero(u)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    return MRIGARKImplicitCache(
        u, uprev, tmp, atmp, v, vtmp, f1eval, kk, z, fS, zemb,
        fsalfirst, k, nlsolver, tab
    )
end

function alg_cache(
        alg::MRIGARKImplicitAlg, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = mri_gark_tableau(alg, eltype(u))
    γ = _mrigark_impl_γ(tab)
    c = one(tTypeNoUnits)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    return MRIGARKImplicitConstantCache(nlsolver, tab)
end
