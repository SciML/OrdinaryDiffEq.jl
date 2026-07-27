# Newton adaptive get_tmp_cache (SDE-specific nlsolver fields):
@inline DiffEqBase.get_tmp_cache(
    integrator, alg::StochasticDiffEqNewtonAdaptiveAlgorithm,
    cache::StochasticDiffEqMutableCache
) = (cache.nlsolver.tmp, cache.nlsolver.ztmp)
@inline DiffEqBase.get_tmp_cache(
    integrator, alg::StochasticCompositeAlgorithm,
    cache::StochasticCompositeCache
) = get_tmp_cache(integrator, alg.algs[1], cache.caches[1])

function full_cache(integrator::StochasticCompositeCache)
    return Iterators.flatten(full_cache(c) for c in integrator.caches)
end

ratenoise_cache(integrator::SDEIntegrator) = ratenoise_cache(integrator.cache)
function ratenoise_cache(integrator::StochasticCompositeCache)
    return Iterators.flatten(ratenoise_cache(c) for c in integrator.caches)
end

rand_cache(integrator::SDEIntegrator) = rand_cache(integrator.cache)
function rand_cache(integrator::StochasticCompositeCache)
    return Iterators.flatten(rand_cache(c) for c in integrator.caches)
end

jac_iter(integrator::SDEIntegrator) = jac_iter(integrator.cache)
function jac_iter(integrator::StochasticCompositeCache)
    return Iterators.flatten(jac_iter(c) for c in integrator.caches)
end

"""
    resize_noise!(integrator, cache, bot_idx, i) -> nothing

Resize the noise process of `integrator` to length `i` after the state was resized.

All of the noise process's buffers (`dW`, `dWtilde`, `dWtmp`, `curW`, the `dZ` family
when [`alg_needs_extra_process`](@ref) holds, and the cached increments in the
interpolation stacks `S₁`/`S₂`) are grown or shrunk together. New entries from
`bot_idx` up to `i` are filled with freshly sampled increments via
[`fill_new_noise_caches!`](@ref) so the process stays consistent with the already
generated path.

This is called from the `resize!` integrator interface; see
[`deleteat_noise!`](@ref) and [`addat_noise!`](@ref) for the index-wise variants.
"""
function resize_noise!(integrator, cache, bot_idx, i)
    for c in integrator.W.S₁
        resize!(c[2], i)
        if alg_needs_extra_process(integrator.alg)
            resize!(c[3], i)
        end
        if i >= bot_idx # fill in rands
            fill_new_noise_caches!(integrator, c, c[1], bot_idx:i)
        end
    end
    for c in integrator.W.S₂
        resize!(c[2], i)
        if alg_needs_extra_process(integrator.alg)
            resize!(c[3], i)
        end
        if i >= bot_idx # fill in rands
            fill_new_noise_caches!(integrator, c, c[1], bot_idx:i)
        end
    end
    resize!(integrator.W.dW, i)
    integrator.W.dW[end] = zero(eltype(integrator.u))
    resize!(integrator.W.dWtilde, i)
    integrator.W.dWtilde[end] = zero(eltype(integrator.u))
    resize!(integrator.W.dWtmp, i)
    integrator.W.dWtmp[end] = zero(eltype(integrator.u))
    resize!(integrator.W.curW, i)
    integrator.W.curW[end] = zero(eltype(integrator.u))
    DiffEqNoiseProcess.resize_stack!(integrator.W, i)

    if alg_needs_extra_process(integrator.alg)
        resize!(integrator.W.dZ, i)
        integrator.W.dZ[end] = zero(eltype(integrator.u))
        resize!(integrator.W.dZtilde, i)
        integrator.W.dZtilde[end] = zero(eltype(integrator.u))
        resize!(integrator.W.dZtmp, i)
        integrator.W.dZtmp[end] = zero(eltype(integrator.u))
        resize!(integrator.W.curZ, i)
        integrator.W.curZ[end] = zero(eltype(integrator.u))
    end
    return if i >= bot_idx # fill in rands
        fill!(@view(integrator.W.curW[bot_idx:i]), zero(eltype(integrator.u)))
        if alg_needs_extra_process(integrator.alg)
            fill!(@view(integrator.W.curZ[bot_idx:i]), zero(eltype(integrator.u)))
        end
    end
end

"""
    fill_new_noise_caches!(integrator, c, scaling_factor, idxs) -> nothing

Sample fresh noise increments into positions `idxs` of a cached noise entry `c`.

`c` is one entry of the noise process's interpolation stacks (`W.S₁`/`W.S₂`), a tuple
whose first element is the step's scaling factor and whose remaining elements hold the
cached `ΔW` and, when [`alg_needs_extra_process`](@ref) holds, `ΔZ` values. The new
values are drawn from the process's own distribution so that the extended path has
the correct law.

Used by [`resize_noise!`](@ref) and [`addat_noise!`](@ref) whenever the state grows.
"""
@inline function fill_new_noise_caches!(integrator, c, scaling_factor, idxs)
    return if isinplace(integrator.W)
        integrator.W.dist(
            @view(c[2][idxs]), integrator.W, scaling_factor,
            integrator.u, integrator.p, integrator.t, integrator.W.rng
        )
        if alg_needs_extra_process(integrator.alg)
            integrator.W.dist(
                @view(c[3][idxs]), integrator.W, scaling_factor,
                integrator.u, integrator.p, integrator.t, integrator.W.rng
            )
        end
    else
        c[2][idxs] .= integrator.noise(length(idxs), integrator, scaling_factor)
        if alg_needs_extra_process(integrator.alg)
            c[3][idxs] .= integrator.noise(length(idxs), integrator, scaling_factor)
        end
    end
end

function resize_non_user_cache!(integrator::SDEIntegrator, cache, i)
    bot_idx = length(integrator.u) + 1
    if is_diagonal_noise(integrator.sol.prob)
        resize_noise!(integrator, cache, bot_idx, i)
        for c in rand_cache(integrator)
            resize!(c, i)
        end
    end
    for c in ratenoise_cache(integrator)
        resize!(c, i)
    end
    return
end

function resize_non_user_cache!(integrator::SDEIntegrator, cache::CompositeCache, i)
    for _cache in cache.caches
        resize_non_user_cache!(integrator, _cache, i)
    end
    return
end

function deleteat_non_user_cache!(integrator::SDEIntegrator, cache, idxs)
    if is_diagonal_noise(integrator.sol.prob)
        deleteat_noise!(integrator, cache, idxs)
        for c in rand_cache(integrator)
            deleteat!(c, idxs)
        end
    end
    for c in ratenoise_cache(integrator)
        deleteat!(c, idxs)
    end
    return
end

function deleteat_non_user_cache!(integrator::SDEIntegrator, cache::CompositeCache, idxs)
    for _cache in cache.caches
        deleteat_non_user_cache!(integrator, _cache, idxs)
    end
    return
end

function addat_non_user_cache!(integrator::SDEIntegrator, cache, idxs)
    if is_diagonal_noise(integrator.sol.prob)
        addat_noise!(integrator, cache, idxs)
        for c in rand_cache(integrator)
            addat!(c, idxs)
        end
    end
    for c in ratenoise_cache(integrator)
        addat!(c, idxs)
    end
    return
end

function addat_non_user_cache!(integrator::SDEIntegrator, cache::CompositeCache, idxs)
    for _cache in cache.caches
        addat_non_user_cache!(integrator, _cache, idxs)
    end
    return
end

"""
    deleteat_noise!(integrator, cache, idxs) -> nothing

Delete the components `idxs` from the noise process of `integrator`.

The counterpart of [`addat_noise!`](@ref), called from `deleteat!` on the integrator:
the same components are removed from every noise buffer and from the cached
increments in the interpolation stacks, keeping the process dimension in step with
the shrunken state.
"""
function deleteat_noise!(integrator, cache, idxs)
    for c in integrator.W.S₁
        deleteat!(c[2], idxs)
        if alg_needs_extra_process(integrator.alg)
            deleteat!(c[3], idxs)
        end
    end
    for c in integrator.W.S₂
        deleteat!(c[2], idxs)
        if alg_needs_extra_process(integrator.alg)
            deleteat!(c[3], idxs)
        end
    end
    deleteat!(integrator.W.dW, idxs)
    deleteat!(integrator.W.dWtilde, idxs)
    deleteat!(integrator.W.dWtmp, idxs)
    deleteat!(integrator.W.curW, idxs)
    DiffEqNoiseProcess.resize_stack!(integrator.W, length(integrator.u))

    return if alg_needs_extra_process(integrator.alg)
        deleteat!(integrator.W.curZ, idxs)
        deleteat!(integrator.W.dZtmp, idxs)
        deleteat!(integrator.W.dZtilde, idxs)
        deleteat!(integrator.W.dZ, idxs)
    end
end

"""
    addat_noise!(integrator, cache, idxs) -> nothing

Insert new components at positions `idxs` in the noise process of `integrator`.

The counterpart of [`deleteat_noise!`](@ref), called from `addat!` on the integrator.
Space is made in every noise buffer and in the cached increments of the interpolation
stacks, and the new slots are filled with freshly sampled increments through
[`fill_new_noise_caches!`](@ref).
"""
function addat_noise!(integrator, cache, idxs)
    for c in integrator.W.S₁
        addat!(c[2], idxs)
        if alg_needs_extra_process(integrator.alg)
            addat!(c[3], idxs)
        end
        fill_new_noise_caches!(integrator, c, c[1], idxs)
    end
    for c in integrator.W.S₂
        addat!(c[2], idxs)
        if alg_needs_extra_process(integrator.alg)
            addat!(c[3], idxs)
        end
        fill_new_noise_caches!(integrator, c, c[1], idxs)
    end

    addat!(integrator.W.dW, idxs)
    integrator.W.dW[idxs] .= zero(eltype(integrator.u))
    addat!(integrator.W.curW, idxs)
    integrator.W.curW[idxs] .= zero(eltype(integrator.u))
    if alg_needs_extra_process(integrator.alg)
        addat!(integrator.W.dZ, idxs)
        integrator.W.dZ[idxs] .= zero(eltype(integrator.u))
        addat!(integrator.W.curZ, idxs)
        integrator.W.curZ[idxs] .= zero(eltype(integrator.u))
    end

    i = length(integrator.u)
    resize!(integrator.W.dWtilde, i)
    resize!(integrator.W.dWtmp, i)
    DiffEqNoiseProcess.resize_stack!(integrator.W, i)
    if alg_needs_extra_process(integrator.alg)
        resize!(integrator.W.dZtmp, i)
        resize!(integrator.W.dZtilde, i)
    end

    # fill in rands
    fill!(@view(integrator.W.curW[idxs]), zero(eltype(integrator.u)))
    return if alg_needs_extra_process(integrator.alg)
        fill!(@view(integrator.W.curZ[idxs]), zero(eltype(integrator.u)))
    end
end
