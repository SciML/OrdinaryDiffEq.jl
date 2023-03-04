## initialize!

@muladd function initialize!(nlsolver::NLSolver{<:NLFunctional},
                             integrator::DiffEqBase.DEIntegrator)
    nlsolver.cache.tstep = integrator.t + nlsolver.c * integrator.dt

    nothing
end

@muladd function initialize!(nlsolver::NLSolver{<:NLAnderson},
                             integrator::DiffEqBase.DEIntegrator)
    @static if VERSION >= 1.8
        (; cache) = nlsolver
    else
        @unpack cache = nlsolver
    end

    cache.history = 0
    cache.tstep = integrator.t + nlsolver.c * integrator.dt

    nothing
end

## initial_η

function initial_η(nlsolver::NLSolver{<:Union{NLFunctional, NLAnderson}}, integrator)
    nlsolver.ηold
end

## compute_step!

"""
    compute_step!(nlsolver::NLSolver{<:Union{NLFunctional,NLAnderson}}, integrator)

Compute the next step of the fixed-point iteration
```math
g(z) = dt⋅f(tmp + γ⋅z, p, t + c⋅dt),
```
and return the norm of ``g(z) - z``.

# References

Ernst Hairer and Gerhard Wanner, "Solving Ordinary Differential
Equations II, Springer Series in Computational Mathematics. ISBN
978-3-642-05221-7. Section IV.8.
[doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7).
"""
function compute_step!(nlsolver::NLSolver{<:NLFunctional}, integrator)
    compute_step_fixedpoint!(nlsolver, integrator)
end

@muladd function compute_step!(nlsolver::NLSolver{<:NLAnderson, false}, integrator)
    @static if VERSION >= 1.8
        (; cache) = nlsolver
    else
        @unpack cache = nlsolver
    end
    @static if VERSION >= 1.8
        (; aa_start) = cache
    else
        @unpack aa_start = cache
    end

    # perform Anderson acceleration
    previter = nlsolver.iter - 1
    if previter == aa_start
        # update cached values for next step of Anderson acceleration
        cache.dzold = cache.dz
        cache.z₊old = nlsolver.z
    elseif previter > aa_start
        # actually perform Anderson acceleration
        nlsolver.z = anderson(nlsolver.z, cache)
        if DiffEqBase.has_destats(integrator)
            integrator.destats.nsolve += 1
        end
    end

    # compute next step
    compute_step_fixedpoint!(nlsolver, integrator)
end

@muladd function compute_step!(nlsolver::NLSolver{<:NLAnderson, true}, integrator)
    @static if VERSION >= 1.8
        (; cache) = nlsolver
    else
        @unpack cache = nlsolver
    end
    @static if VERSION >= 1.8
        (; aa_start) = cache
    else
        @unpack aa_start = cache
    end

    # perform Anderson acceleration
    previter = nlsolver.iter - 1
    if previter == aa_start
        # update cached values for next step of Anderson acceleration
        @.. broadcast=false cache.dzold=cache.dz
        @.. broadcast=false cache.z₊old=nlsolver.z
    elseif previter > aa_start
        # actually perform Anderson acceleration
        anderson!(nlsolver.z, cache)
        if DiffEqBase.has_destats(integrator)
            integrator.destats.nsolve += 1
        end
    end

    # compute next step
    compute_step_fixedpoint!(nlsolver, integrator)
end

@muladd function compute_step_fixedpoint!(nlsolver::NLSolver{
                                                             <:Union{NLFunctional,
                                                                     NLAnderson}, false},
                                          integrator)
    @static if VERSION >= 1.8
        (; uprev, t, p, dt, opts) = integrator
    else
        @unpack uprev, t, p, dt, opts = integrator
    end
    @static if VERSION >= 1.8
        (; z, γ, α, cache, tmp) = nlsolver
    else
        @unpack z, γ, α, cache, tmp = nlsolver
    end
    @static if VERSION >= 1.8
        (; tstep) = cache
    else
        @unpack tstep = cache
    end

    f = nlsolve_f(integrator)
    isdae = f isa DAEFunction

    γdt = γ * dt
    if isdae
        ustep = @.. broadcast=false uprev+z
        invγdt = inv(γdt)
        dustep = @.. broadcast=false (tmp + α * z)*invγdt
        dz = f(dustep, ustep, p, t)
        ztmp = @.. broadcast=false z+dz
    else
        mass_matrix = integrator.f.mass_matrix
        if nlsolver.method === COEFFICIENT_MULTISTEP
            ustep = z
            if mass_matrix === I
                ztmp = (tmp .+ f(z, p, tstep)) * (γdt / α)
                dz = ztmp .- z
            else
                ztmp = _reshape(mass_matrix * _vec(z), axes(z))
                dz = (tmp .+ f(z, p, tstep)) * γdt - α .* ztmp
                ztmp = dz .+ z
            end
        else
            ustep = @.. broadcast=false tmp+γ * z
            if mass_matrix === I
                ztmp = dt .* f(ustep, p, tstep)
                dz = ztmp .- z
            else
                ztmp = _reshape(mass_matrix * _vec(z), axes(z))
                dz = dt .* f(ustep, p, tstep) .- ztmp
                ztmp = z .+ dz
            end
        end
    end
    if DiffEqBase.has_destats(integrator)
        integrator.destats.nf += 1
    end

    # compute norm of residuals
    atmp = calculate_residuals(dz, uprev, ustep, opts.abstol, opts.reltol,
                               opts.internalnorm, t)
    ndz = opts.internalnorm(atmp, t)

    # cache results
    nlsolver.ztmp = ztmp
    if isdefined(cache, :dz)
        cache.dz = dz
    end

    ndz
end

@muladd function compute_step_fixedpoint!(nlsolver::NLSolver{
                                                             <:Union{NLFunctional,
                                                                     NLAnderson}, true},
                                          integrator)
    @static if VERSION >= 1.8
        (; uprev, t, p, dt, opts) = integrator
    else
        @unpack uprev, t, p, dt, opts = integrator
    end
    @static if VERSION >= 1.8
        (; z, tmp, ztmp, γ, α, cache) = nlsolver
    else
        @unpack z, tmp, ztmp, γ, α, cache = nlsolver
    end
    @static if VERSION >= 1.8
        (; ustep, tstep, k, atmp, dz) = cache
    else
        @unpack ustep, tstep, k, atmp, dz = cache
    end

    f = nlsolve_f(integrator)
    isdae = f isa DAEFunction

    γdt = γ * dt
    if isdae
        @.. broadcast=false ustep=uprev + z
        @.. broadcast=false ztmp=(tmp + α * z) * inv(γdt)
        f(k, ztmp, ustep, p, tstep)
        @.. broadcast=false dz=k
        @.. broadcast=false ztmp=z + dz
    else
        mass_matrix = integrator.f.mass_matrix
        if nlsolver.method === COEFFICIENT_MULTISTEP
            ustep = z
            f(k, ustep, p, tstep)
            if mass_matrix === I
                @.. broadcast=false ztmp=(tmp + k) * (γdt / α)
                @.. broadcast=false dz=ztmp - z
            else
                update_coefficients!(mass_matrix, ustep, p, tstep)
                mul!(_vec(ztmp), mass_matrix, _vec(z))
                @.. broadcast=false dz=(tmp + k) * γdt - α * ztmp
                @.. broadcast=false ztmp=dz + z
            end
        else
            @.. broadcast=false ustep=tmp + γ * z
            f(k, ustep, p, tstep)
            if mass_matrix === I
                @.. broadcast=false ztmp=dt * k
                @.. broadcast=false dz=ztmp - z
            else
                update_coefficients!(mass_matrix, ustep, p, tstep)
                mul!(_vec(ztmp), mass_matrix, _vec(z))
                @.. broadcast=false dz=dt * k - ztmp
                @.. broadcast=false ztmp=z + dz
            end
        end
    end

    if DiffEqBase.has_destats(integrator)
        integrator.destats.nf += 1
    end

    # compute norm of residuals
    calculate_residuals!(atmp, dz, uprev, ustep, opts.abstol, opts.reltol,
                         opts.internalnorm, t)
    ndz = opts.internalnorm(atmp, t)

    ndz
end

## resize!

function Base.resize!(nlcache::NLFunctionalCache, i::Int)
    resize!(nlcache.ustep, i)
    resize!(nlcache.k, i)
    resize!(nlcache.atmp, i)
    resize!(nlcache.dz, i)
    nothing
end

function Base.resize!(nlcache::NLAndersonCache, nlsolver::NLSolver{<:NLAnderson},
                      integrator, i::Int)
    resize!(nlcache, nlsolver.alg, i)
end

function Base.resize!(nlcache::NLAndersonCache, nlalg::NLAnderson, i::Int)
    @static if VERSION >= 1.8
        (; z₊old, Δz₊s) = nlcache
    else
        @unpack z₊old, Δz₊s = nlcache
    end

    resize!(nlcache.ustep, i)
    resize!(nlcache.k, i)
    resize!(nlcache.atmp, i)
    resize!(nlcache.dz, i)
    resize!(nlcache.dzold, i)
    resize!(z₊old, i)

    # update history of Anderson cache
    max_history_old = length(Δz₊s)
    max_history = min(nlalg.max_history, nlalg.max_iter, i)

    resize!(nlcache.γs, max_history)
    resize!(nlcache.Δz₊s, max_history)

    if max_history != max_history_old
        nlcache.Q = typeof(nlcache.Q)(undef, i, max_history)
        nlcache.R = typeof(nlcache.R)(undef, max_history, max_history)
    end

    max_history = length(Δz₊s)
    if max_history > max_history_old
        for i in (max_history_old + 1):max_history
            Δz₊s[i] = zero(z₊old)
        end
    end

    nothing
end
