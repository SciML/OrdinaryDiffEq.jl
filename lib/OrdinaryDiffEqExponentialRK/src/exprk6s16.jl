const _EXPRK6S16_C = (
    1 / 2,
    1 / 2,
    1 / 3,
    1 / 2,
    1 / 5,
    1 / 4,
    1 / 2,
    1 / 5,
    1 / 4,
    1 / 3,
    1 / 2,
    1 / 5,
    1 / 4,
    1 / 3,
    1.0,
)

function expRK_operators(::ExpRK6s16, dt, A)
    return ntuple(i -> phi((i == 16 ? one(dt) : _EXPRK6S16_C[i]) * dt * A, 6), 16)
end

struct ExpRK6s16ConstantCache{OpType,FType} <: ExpRKConstantCache
    ops::OpType
    uf::FType
end

function alg_cache(
    alg::ExpRK6s16,
    u,
    rate_prototype,
    ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits},
    ::Type{tTypeNoUnits},
    uprev,
    uprev2,
    f,
    t,
    dt,
    reltol,
    p,
    calck,
    ::Val{false},
    verbose,
) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    ops =
        alg.krylov ? nothing :
        begin
            isa(f, SplitFunction) ||
                throw(ArgumentError("Caching can only be used with SplitFunction"))
            A =
                size(f.f1.f) == () ? convert(Number, f.f1.f) :
                convert(AbstractMatrix, f.f1.f)
            expRK_operators(alg, dt, A)
        end
    uf =
        isa(f, SplitFunction) || SciMLBase.has_jac(f) ? nothing :
        UDerivativeWrapper(f, t, p)
    return ExpRK6s16ConstantCache(ops, uf)
end

@cache struct ExpRK6s16Cache{uType,rateType,JCType,FType,JType,OpType,KsType} <: ExpRKCache
    u::uType
    uprev::uType
    tmp::uType
    dz::uType
    rtmp::rateType
    g1::rateType
    Au::rateType
    du1::rateType
    Ds::Vector{rateType}
    jac_config::JCType
    uf::FType
    J::JType
    ops::OpType
    KsCache::KsType
end

function alg_cache(
    alg::ExpRK6s16,
    u,
    rate_prototype,
    ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits},
    ::Type{tTypeNoUnits},
    uprev,
    uprev2,
    f,
    t,
    dt,
    reltol,
    p,
    calck,
    ::Val{true},
    verbose,
) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp, dz = (zero(u) for _ = 1:2)
    rtmp, g1, Au, du1 = (zero(rate_prototype) for _ = 1:4)
    Ds = [zero(rate_prototype) for _ = 1:15]
    uf, jac_config, J, ops, KsCache =
        alg_cache_expRK(alg, u, uEltypeNoUnits, uprev, f, t, dt, p, du1, tmp, dz, (6,))
    return ExpRK6s16Cache(
        u,
        uprev,
        tmp,
        dz,
        rtmp,
        g1,
        Au,
        du1,
        Ds,
        jac_config,
        uf,
        J,
        ops,
        KsCache,
    )
end

@inline function _exprk6s16_coefficients(c, others, order; scaled = one(c))
    denominator = c * prod(c - other for other in others)
    numerator = if order == 4
        (prod(others), -2sum(others), 6)
    elseif order == 5
        (
            -prod(others),
            2sum(
                others[i] * others[j] for i in eachindex(others) for
                j = (i+1):length(others)
            ),
            -6sum(others),
            24,
        )
    else
        (
            prod(others),
            -2sum(
                others[i] * others[j] * others[k] for i in eachindex(others) for
                j = (i+1):length(others) for k = (j+1):length(others)
            ),
            6sum(
                others[i] * others[j] for i in eachindex(others) for
                j = (i+1):length(others)
            ),
            -24sum(others),
            120,
        )
    end
    return [scaled^q * numerator[q-1] / denominator for q = 2:order]
end

function _exprk6s16_stage(uprev, F1, Ds, stage, terms, phi_action, dt)
    c = stage == 16 ? one(dt) : _EXPRK6S16_C[stage]
    u = uprev + dt * c * phi_action(F1, stage, 1)
    for (D, coefficients) in terms
        for (index, coefficient) in enumerate(coefficients)
            u += dt * coefficient * phi_action(D, stage, index + 1)
        end
    end
    return u
end

function _exprk6s16_terms(Ds, stage)
    c = _EXPRK6S16_C
    if stage == 1
        return ()
    elseif stage in 2:3
        return [(Ds[1], [c[stage]^2 / c[1]])]
    elseif stage in 4:6
        sources = (2, 3)
        return [
            (
                Ds[j],
                [
                    -c[stage]^2 * c[k] / (c[j] * (c[j] - c[k])),
                    2c[stage]^3 / (c[j] * (c[j] - c[k])),
                ],
            ) for j in sources for k in sources if j != k
        ]
    elseif stage in 7:10
        sources = 4:6
        return [
            (
                Ds[j],
                _exprk6s16_coefficients(
                    c[j],
                    [c[k] for k in sources if k != j],
                    4;
                    scaled = c[stage],
                ),
            ) for j in sources
        ]
    elseif stage in 11:15
        sources = 7:10
        return [
            (
                Ds[j],
                _exprk6s16_coefficients(
                    c[j],
                    [c[k] for k in sources if k != j],
                    5;
                    scaled = c[stage],
                ),
            ) for j in sources
        ]
    end
    return [
        (Ds[j], _exprk6s16_coefficients(c[j], [c[k] for k = 11:15 if k != j], 6)) for
        j = 11:15
    ]
end

function _exprk6s16_step(uprev, F1, nonlinear, A, alg, integrator, dt, t, Ds, g1, ops)
    phi_action = if alg.krylov
        kwargs = (
            m = min(alg.m, size(A, 1)),
            opnorm = integrator.opts.internalopnorm,
            iop = alg.iop,
        )
        (v, stage, q) -> @view(
            phiv((stage == 16 ? one(dt) : _EXPRK6S16_C[stage]) * dt, A, v, q; kwargs...)[
                :,
                q+1,
            ]
        )
    else
        (v, stage, q) -> ops[stage][q+1] * v
    end
    for stage = 1:15
        U = _exprk6s16_stage(
            uprev,
            F1,
            Ds,
            stage,
            _exprk6s16_terms(Ds, stage),
            phi_action,
            dt,
        )
        Ds[stage] = nonlinear(U, t + _EXPRK6S16_C[stage] * dt) - g1
    end
    return _exprk6s16_stage(uprev, F1, Ds, 16, _exprk6s16_terms(Ds, 16), phi_action, dt)
end

function perform_step!(integrator, cache::ExpRK6s16ConstantCache, repeat_step = false)
    (; t, dt, uprev, f, p) = integrator
    A = isa(f, SplitFunction) ? f.f1.f : calc_J(integrator, cache)
    alg = unwrap_alg(integrator, true)
    F1 = integrator.fsalfirst
    nonlinear = (u, time) -> _compute_nl(f, u, p, time, A)
    g1 = nonlinear(uprev, t)
    Ds = Vector{typeof(g1)}(undef, 15)
    u = _exprk6s16_step(uprev, F1, nonlinear, A, alg, integrator, dt, t, Ds, g1, cache.ops)
    if isa(f, SplitFunction)
        integrator.stats.nf2 += 15
    else
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 15)
    end
    integrator.fsallast = f(u, p, t + dt)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    return integrator.u = u
end

function perform_step!(integrator, cache::ExpRK6s16Cache, repeat_step = false)
    (; t, dt, uprev, u, f, p) = integrator
    (; g1, Au, J, Ds, ops) = cache
    A = isa(f, SplitFunction) ? f.f1.f : (calc_J!(J, integrator, cache); J)
    alg = unwrap_alg(integrator, true)
    F1 = integrator.fsalfirst
    _compute_nl!(g1, f, uprev, p, t, A, Au)
    nonlinear = function (state, time)
        G = zero(g1)
        _compute_nl!(G, f, state, p, time, A, Au)
        return G
    end
    result = _exprk6s16_step(uprev, F1, nonlinear, A, alg, integrator, dt, t, Ds, g1, ops)
    u .= result
    if isa(f, SplitFunction)
        integrator.stats.nf2 += 15
    else
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 15)
    end
    f(integrator.fsallast, u, p, t + dt)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end
