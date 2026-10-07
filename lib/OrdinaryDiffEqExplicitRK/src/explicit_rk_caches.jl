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

# When both the dimensionless time type and tableau `c` are IEEE floats, match
# `c` to the time type so `t + c[i]*dt` keeps FunctionWrapper signatures under
# AutoSpecialize. Otherwise use promote_type (preserves BigFloat, Rational, etc.).
const _IEEEFloat = Union{Float16, Float32, Float64}

function _explicit_rk_c_eltype(::Type{T}, ::Type{C}) where {T, C}
    return (T <: _IEEEFloat && C <: _IEEEFloat) ? T : promote_type(C, T)
end

function ExplicitRKConstantCache(tableau, rate_prototype, ::Type{tType} = Float64) where {tType}
    (; A, c, α, αEEst, stages) = tableau
    A = copy(A') # Transpose A to column major looping
    # `one(tType)` strips units for Unitful while preserving BigFloat/Float32.
    cType = _explicit_rk_c_eltype(typeof(one(tType)), eltype(c))
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
