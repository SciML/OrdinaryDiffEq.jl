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

# Prefer the dimensionless time type for tableau `c` when conversion is exact
# (and for IEEE float narrowing), so `t + c[i]*dt` keeps FunctionWrapper types.
# Fall back to promote_type to preserve wider coefficient types (e.g. BigFloat).
_is_ieee_float(::Type{T}) where {T} = false
_is_ieee_float(::Type{Float16}) = true
_is_ieee_float(::Type{Float32}) = true
_is_ieee_float(::Type{Float64}) = true

function _explicit_rk_c_eltype(::Type{T}, c) where {T}
    C = eltype(c)
    preferred = T
    fallback = promote_type(C, T)
    preferred === fallback && return preferred
    if _is_ieee_float(preferred) && _is_ieee_float(C)
        return preferred
    elseif _coefficients_exactly_convertible(preferred, c)
        return preferred
    else
        return fallback
    end
end

_exactly_convertible(::Type{T}, ::T) where {T} = true

function _exactly_convertible(::Type{T}, x::AbstractFloat) where {T <: AbstractFloat}
    y = T(x)
    return convert(typeof(x), y) == x
end

function _exactly_convertible(::Type{Rational{BigInt}}, x::AbstractFloat)
    return isfinite(x)
end

function _exactly_convertible(::Type{Rational{I}}, x::AbstractFloat) where {I <: Integer}
    isfinite(x) || return false
    r = Rational{BigInt}(x)
    return (typemin(I) <= r.num <= typemax(I)) & (0 < r.den <= typemax(I))
end

function _exactly_convertible(::Type{Rational{BigInt}}, ::Rational)
    return true
end

function _exactly_convertible(::Type{Rational{I}}, x::Rational) where {I <: Integer}
    return (typemin(I) <= numerator(x) <= typemax(I)) &
        (0 < denominator(x) <= typemax(I))
end

_exactly_convertible(::Type{T}, x) where {T} = false

function _coefficients_exactly_convertible(::Type{T}, c) where {T}
    @inbounds for i in eachindex(c)
        _exactly_convertible(T, c[i]) || return false
    end
    return true
end

function ExplicitRKConstantCache(tableau, rate_prototype, ::Type{tType} = Float64) where {tType}
    (; A, c, α, αEEst, stages) = tableau
    A = copy(A') # Transpose A to column major looping
    # `one(tType)` strips units for Unitful while preserving BigFloat/Float32.
    cType = _explicit_rk_c_eltype(typeof(one(tType)), c)
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
