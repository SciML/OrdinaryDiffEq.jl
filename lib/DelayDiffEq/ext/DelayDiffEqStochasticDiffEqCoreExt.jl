module DelayDiffEqStochasticDiffEqCoreExt

using DelayDiffEq
import DelayDiffEq: _sde_alg_order, _sde_isadaptive, _sde_alg_cache, _create_sdde_noise,
    SDEAlgUnion

import StochasticDiffEqCore
using StochasticDiffEqCore: alg_cache as sde_alg_cache_impl

using DiffEqNoiseProcess: WienerProcess, WienerProcess!, RSWM
using DiffEqBase: is_diagonal_noise
using Random: Random

using OrdinaryDiffEqCore: DEVerbosity
using SciMLBase: isinplace

# ── Trait bridges ──────────────────────────────────────────────────────

function DelayDiffEq._sde_alg_order(alg::SDEAlgUnion)
    return StochasticDiffEqCore.alg_order(alg)
end

function DelayDiffEq._sde_isadaptive(alg::SDEAlgUnion)
    return StochasticDiffEqCore.isadaptive(alg)
end

# ── SDE alg_cache wrapper ─────────────────────────────────────────────

function DelayDiffEq._sde_alg_cache(
        alg, prob, u, dW, dZ, p,
        rate_prototype, noise_rate_prototype,
        jump_rate_prototype,
        uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, f, t, dt,
        iip_val, verbose
    )
    return sde_alg_cache_impl(
        alg, prob, u, dW, dZ, p,
        rate_prototype, noise_rate_prototype,
        jump_rate_prototype,
        uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, f, t, dt,
        iip_val, verbose
    )
end

# ── Noise creation ─────────────────────────────────────────────────────

function DelayDiffEq._create_sdde_noise(
        prob, u0, t0, dt, tdir, noise_rate_prototype,
        save_noise, seed, iip, adaptive
    )
    _seed = seed == 0 ? rand(UInt64) : seed
    _rng = Random.Xoshiro(_seed)

    if is_diagonal_noise(prob)
        rand_prototype = u0 isa Number ? u0 : zero(u0)
    elseif noise_rate_prototype !== nothing
        rand_prototype = false .* noise_rate_prototype[1, :]
    else
        rand_prototype = u0 isa Number ? u0 : zero(u0)
    end

    if prob.noise === nothing
        if iip
            W = WienerProcess!(t0, rand_prototype, save_everystep = save_noise, rng = _rng)
        else
            W = WienerProcess(t0, rand_prototype, save_everystep = save_noise, rng = _rng)
        end
    else
        W = prob.noise
        if W.curt != t0
            reinit!(W, t0, t0 = t0)
        end
    end

    P = nothing
    tType = typeof(dt)
    sqdt = tdir * sqrt(abs(dt))

    return W, P, sqdt
end

end # module
