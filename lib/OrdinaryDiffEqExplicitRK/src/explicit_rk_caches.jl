@cache struct ExplicitRKCache{uType, rateType, TabType, TmpC <: TmpCache} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    fsallast::rateType
    kk::Vector{rateType}
    tab::TabType
    # Unified scratch: `tmp` (stage scratch), `tmp2` (replaces the former inline
    # `utilde` accumulation buffer) and `atmp` (error-norm scaling). Default
    # layout is `TmpCache{uType, Nothing, uNoUnitsType, Nothing}` (rate slots
    # skipped).
    tmp_cache::TmpC
end

get_fsalfirstlast(cache::ExplicitRKCache, u) = (cache.kk[1], cache.fsallast)

function alg_cache(
        alg::ExplicitRK, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    kk = Vector{typeof(rate_prototype)}(undef, 0)
    for i in 1:(alg.tableau.stages)
        push!(kk, zero(rate_prototype))
    end
    fsalfirst = kk[1]
    if isfsal(alg.tableau)
        fsallast = kk[end]
    else
        fsallast = zero(rate_prototype)
    end
    # The historical inline scratch (`tmp`, `utilde`, `atmp`) migrates into the
    # unified `TmpCache`, so the array count matches the historical cache
    # exactly (the former rate-typed `utilde` allocation becomes the state-typed
    # `tmp2`; the >17-stage fallback accumulators fold `dt` in so the buffer
    # always carries state units). There are no safe rate donors — `kk` feeds
    # the generic RK interpolant via `_ode_addsteps!` and `fsalfirst`/`fsallast`
    # feed the Hermite interpolant — and `ExplicitRK` is not a `@kwdef`
    # algorithm (no `preallocate_initdt_buffers` field), so the rate slots stay
    # `nothing` and `initdt` allocates its rate temporaries at call time.
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    tab = ExplicitRKConstantCache(alg.tableau, rate_prototype, typeof(dt))
    return ExplicitRKCache(u, uprev, fsalfirst, fsallast, kk, tab, tmp_cache)
end

struct ExplicitRKConstantCache{MType, VType, CType, KType, BType, BiType} <:
    OrdinaryDiffEqConstantCache
    A::MType
    c::CType
    α::VType
    αEEst::VType
    stages::Int
    kk::KType
    B_interp::BType
    bi::BiType  # Pre-allocated buffer for interpolation polynomial weights
end

function ExplicitRKConstantCache(tableau, rate_prototype, ::Type{tType} = Float64) where {tType}
    (; A, c, α, αEEst, stages) = tableau
    A = copy(A') # Transpose A to column major looping
    # Convert c to match the dimensionless numeric type of dt so that
    # t + c[i]*dt doesn't promote t beyond what FunctionWrapper signatures
    # expect under AutoSpecialize. `one(tType)` strips units for Unitful
    # quantities while preserving precision for BigFloat/Float32 time types.
    # Use promote_type so that when `c`'s eltype is already wider than the
    # dimensionless dt type (e.g. BigFloat `c` with Rational{Int} dt during
    # convergence testing), we don't attempt a lossy narrowing that could
    # throw InexactError for non-exactly-representable coefficients.
    cType = promote_type(eltype(c), typeof(one(tType)))
    c = cType.(c)
    kk = Array{typeof(rate_prototype)}(undef, stages) # Not ks since that's for integrator.opts.dense
    αEEst = isempty(αEEst) ? αEEst : α .- αEEst
    B_interp = hasproperty(tableau, :B_interp) ? tableau.B_interp : nothing
    bi = if isnothing(B_interp)
        nothing
    else
        Vector{eltype(B_interp)}(undef, size(B_interp, 1))
    end
    return ExplicitRKConstantCache(A, c, α, αEEst, stages, kk, B_interp, bi)
end

function alg_cache(
        alg::ExplicitRK, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return ExplicitRKConstantCache(alg.tableau, rate_prototype, typeof(dt))
end
