@cache struct ExplicitRKCache{uType, rateType, uNoUnitsType, TabType} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    utilde::rateType
    atmp::uNoUnitsType
    fsalfirst::rateType
    fsallast::rateType
    kk::Vector{rateType}
    tab::TabType
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
    utilde = zero(rate_prototype)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tab = ExplicitRKConstantCache(alg.tableau, rate_prototype, typeof(dt))
    return ExplicitRKCache(u, uprev, tmp, utilde, atmp, fsalfirst, fsallast, kk, tab)
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
