abstract type LinearMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::LinearMutableCache, u) = (cache.fsalfirst, cache.k)

@cache struct MagnusMidpointCache{uType, rateType, WType, expType} <:
    LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::MagnusMidpoint, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return MagnusMidpointCache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct MagnusMidpointConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::MagnusMidpoint, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MagnusMidpointConstantCache()
end

@cache struct RKMK2Cache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::RKMK2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return RKMK2Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct RKMK2ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::RKMK2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RKMK2ConstantCache()
end

@cache struct LieRK4Cache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::LieRK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return LieRK4Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct LieRK4ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::LieRK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return LieRK4ConstantCache()
end

@cache struct CG3Cache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::CG3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return CG3Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct CG3ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::CG3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return CG3ConstantCache()
end

@cache struct CG2Cache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::CG2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return CG2Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct CG2ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::CG2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return CG2ConstantCache()
end

@cache struct CG4aCache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::CG4a, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, ::Val{true}, verbose
    )
    W = false .* vec(rate_prototype) .* vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return CG4aCache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct CG4aConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::CG4a, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, ::Val{false}, verbose
    )
    return CG4aConstantCache()
end

@cache struct RKMK4Cache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::RKMK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return RKMK4Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct RKMK4ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::RKMK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return RKMK4ConstantCache()
end

@cache struct MagnusAdapt4Cache{uType, rateType, WType, uNoUnitsType, expType} <:
    LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    utilde::uType
    atmp::uNoUnitsType
    exp_cache::expType
end

function alg_cache(
        alg::MagnusAdapt4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return MagnusAdapt4Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, utilde, atmp, exp_cache)
end

struct MagnusAdapt4ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::MagnusAdapt4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MagnusAdapt4ConstantCache()
end

@cache struct MagnusNC8Cache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::MagnusNC8, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return MagnusNC8Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct MagnusNC8ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::MagnusNC8, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MagnusNC8ConstantCache()
end

@cache struct MagnusGL4Cache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::MagnusGL4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return MagnusGL4Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct MagnusGL4ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::MagnusGL4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MagnusGL4ConstantCache()
end

@cache struct MagnusGL8Cache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::MagnusGL8, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return MagnusGL8Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct MagnusGL8ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::MagnusGL8, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MagnusGL8ConstantCache()
end

@cache struct MagnusNC6Cache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::MagnusNC6, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return MagnusNC6Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct MagnusNC6ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::MagnusNC6, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MagnusNC6ConstantCache()
end

@cache struct MagnusGL6Cache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::MagnusGL6, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return MagnusGL6Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct MagnusGL6ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::MagnusGL6, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MagnusGL6ConstantCache()
end
@cache struct MagnusGauss4Cache{uType, rateType, WType, expType} <:
    LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::MagnusGauss4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return MagnusGauss4Cache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct MagnusGauss4ConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::MagnusGauss4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MagnusGauss4ConstantCache()
end

@cache struct LieEulerCache{uType, rateType, WType, expType} <: LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::LieEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return LieEulerCache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct LieEulerConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::LieEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return LieEulerConstantCache()
end

@cache struct CayleyEulerCache{uType, rateType} <: LinearMutableCache
    u::uType
    uprev::uType
    tmp::uType
    V::uType
    fsalfirst::rateType
    k::rateType
end

function alg_cache(
        alg::CayleyEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    return CayleyEulerCache(u, uprev, zero(u), zero(u), fsalfirst, k)
end

struct CayleyEulerConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::CayleyEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return CayleyEulerConstantCache()
end

@cache struct MagnusLeapfrogCache{uType, rateType, WType, expType} <:
    LinearMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    tmp::uType
    fsalfirst::rateType
    W::WType
    k::rateType
    exp_cache::expType
end

function alg_cache(
        alg::MagnusLeapfrog, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    W = false .* _vec(rate_prototype) .* _vec(rate_prototype)' # uEltype?
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return MagnusLeapfrogCache(u, uprev, uprev2, zero(u), fsalfirst, W, k, exp_cache)
end

struct MagnusLeapfrogConstantCache <: OrdinaryDiffEqConstantCache
end

function alg_cache(
        alg::MagnusLeapfrog, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return MagnusLeapfrogConstantCache()
end

struct LinearExponentialConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(
        alg::LinearExponential, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return LinearExponentialConstantCache()
end

@cache struct LinearExponentialCache{uType, rateType, KsType, expType} <:
    LinearMutableCache
    u::uType
    uprev::uType
    tmp::uType
    rtmp::rateType
    KsCache::KsType # different depending on alg.krylov
    exp_cache::expType
end

get_fsalfirstlast(cache::LinearExponentialCache, u) = (zero(u), zero(u))

function _phiv_timestep_caches(u_prototype, maxiter::Int, p::Int)
    n = length(u_prototype)
    T = eltype(u_prototype)
    u = zero(u_prototype)                         # stores the current state
    W = similar(u_prototype, n, p + 1)              # stores the w vectors
    P = similar(u_prototype, n, p + 2)              # stores output from phiv!
    Ks = KrylovSubspace{T, T, typeof(similar(u_prototype, size(u_prototype, 1), 2))}(
        n, maxiter
    ) # stores output from arnoldi!
    phiv_cache = PhivCache(u_prototype, maxiter, p + 1) # cache used by phiv! (need +1 for error estimation)
    return u, W, P, Ks, phiv_cache
end

function alg_cache(
        alg::LinearExponential, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    rtmp = zero(rate_prototype)
    n = length(u)
    T = eltype(u)
    m = min(alg.m, n)

    if alg.krylov == :off
        KsCache = nothing
    elseif alg.krylov == :simple
        Ks = KrylovSubspace{T, T, typeof(similar(u, size(u, 1), 2))}(n, m)
        expv_cache = ExpvCache{T}(m)
        KsCache = (Ks, expv_cache)
    elseif alg.krylov == :adaptive
        KsCache = _phiv_timestep_caches(u, m, 0)
    else
        throw(ArgumentError("Unknown krylov setting $(alg.krylov). Can be :off, :simple or :adaptive."))
    end
    exp_cache = ExponentialUtilities.alloc_mem(f, ExpMethodGeneric())
    return LinearExponentialCache(u, uprev, tmp, rtmp, KsCache, exp_cache)
end
