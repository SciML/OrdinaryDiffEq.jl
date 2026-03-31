module DelayDiffEqStochasticDiffEqCoreExt

using DelayDiffEq
import DelayDiffEq: _sde_alg_cache, _create_sdde_noise

import StochasticDiffEqCore
using StochasticDiffEqCore: alg_cache as sde_alg_cache_impl,
    alg_needs_extra_process, _z_prototype

using DiffEqNoiseProcess: WienerProcess, WienerProcess!, RSWM
using DiffEqBase: is_diagonal_noise
using Random: Random

using SciMLBase: isinplace

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
        prob, alg, u0, t0, dt, tdir, noise_rate_prototype,
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

    needs_dZ = alg_needs_extra_process(alg)

    if prob.noise === nothing
        if iip
            if needs_dZ
                rand_prototype2 = _z_prototype(alg, rand_prototype, true)
                W = WienerProcess!(
                    t0, rand_prototype, rand_prototype2,
                    save_everystep = save_noise, rng = _rng
                )
            else
                W = WienerProcess!(
                    t0, rand_prototype,
                    save_everystep = save_noise, rng = _rng
                )
            end
        else
            if needs_dZ
                rand_prototype2 = _z_prototype(alg, rand_prototype, false)
                W = WienerProcess(
                    t0, rand_prototype, rand_prototype2,
                    save_everystep = save_noise, rng = _rng
                )
            else
                W = WienerProcess(
                    t0, rand_prototype,
                    save_everystep = save_noise, rng = _rng
                )
            end
        end
    else
        W = prob.noise
        if needs_dZ && (!hasproperty(W, :dZ) || W.dZ === nothing)
            error("Higher order SDE solver requires extra Brownian process Z. Use `WienerProcess(t, W0, Z0)` instead of `WienerProcess(t, W0)`.")
        end
        if W.curt != t0
            reinit!(W, t0, t0 = t0)
        end
    end

    P = nothing
    sqdt = tdir * sqrt(abs(dt))

    return W, P, sqdt
end

end # module
