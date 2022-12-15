## initialize!

@muladd function initialize!(nlsolver::NLSolver{<:NLNewton, false},
                             integrator::DiffEqBase.DEIntegrator)
    @unpack dt = integrator
    @unpack cache = nlsolver

    cache.invγdt = inv(dt * nlsolver.γ)
    cache.tstep = integrator.t + nlsolver.c * dt

    nothing
end

@muladd function initialize!(nlsolver::NLSolver{<:NLNewton, true},
                             integrator::DiffEqBase.DEIntegrator)
    @unpack u, uprev, t, dt, opts = integrator
    @unpack cache = nlsolver
    @unpack weight = cache

    cache.invγdt = inv(dt * nlsolver.γ)
    cache.tstep = integrator.t + nlsolver.c * dt
    calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, u,
                         opts.abstol, opts.reltol, opts.internalnorm, t)

    nothing
end

## compute_step!

"""
    compute_step!(nlsolver::NLSolver{<:NLNewton}, integrator)

Compute next iterate of numerically stable modified Newton iteration
that is specialized for implicit methods.

Please check
https://github.com/SciML/DiffEqDevMaterials/blob/master/newton/output/main.pdf
for more details.

# References

M.E.Hoseaa and L.F.Shampine, "Analysis and implementation of TR-BDF2",
Applied Numerical Mathematics, Volume 20, Issues 1–2, February 1996, Pages
21-37.
[doi:10.1016/0168-9274(95)00115-8](https://doi.org/10.1016/0168-9274(95)00115-8).

Ernst Hairer and Gerhard Wanner, "Solving Ordinary Differential
Equations II, Springer Series in Computational Mathematics. ISBN
978-3-642-05221-7. Section IV.8.
[doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7).
"""
@muladd function compute_step!(nlsolver::NLSolver{<:NLNewton, false}, integrator)
    @unpack uprev, t, p, dt, opts = integrator
    @unpack z, tmp, γ, α, cache = nlsolver
    @unpack tstep, W, invγdt = cache

    f = nlsolve_f(integrator)

    ztmp, ustep = _compute_rhs(nlsolver, integrator, f, z)

    if DiffEqBase.has_destats(integrator)
        integrator.destats.nf += 1
    end

    # update W
    if W isa DiffEqBase.AbstractDiffEqLinearOperator
        W = update_coefficients!(W, ustep, p, tstep)
    end

    dz = _reshape(W \ _vec(ztmp), axes(ztmp))
    dz = relax(dz, nlsolver, integrator, f)
    if DiffEqBase.has_destats(integrator)
        integrator.destats.nsolve += 1
    end

    atmp = calculate_residuals(dz, uprev, ustep, opts.abstol, opts.reltol,
                               opts.internalnorm, t)
    ndz = opts.internalnorm(atmp, t)
    # NDF and BDF are special because the truncation error is directly
    # proportional to the total displacement.
    if integrator.alg isa QNDF
        ndz *= error_constant(integrator, alg_order(integrator.alg))
    end

    # compute next iterate
    nlsolver.ztmp = z .- dz

    ndz
end

@muladd function compute_step!(nlsolver::NLSolver{<:NLNewton, true}, integrator)
    @unpack uprev, t, p, dt, opts = integrator
    @unpack z, tmp, ztmp, γ, α, iter, cache = nlsolver
    @unpack W_γdt, ustep, tstep, k, atmp, dz, W, new_W, invγdt, linsolve, weight = cache

    f = nlsolve_f(integrator)
    isdae = f isa DAEFunction

    if DiffEqBase.has_destats(integrator)
        integrator.destats.nf += 1
    end

    b, ustep = _compute_rhs!(nlsolver, integrator, f, z)

    # update W
    if W isa DiffEqBase.AbstractDiffEqLinearOperator
        update_coefficients!(W, ustep, p, tstep)
    end

    if integrator.opts.adaptive
        reltol = integrator.opts.reltol
    else
        reltol = eps(eltype(dz))
    end

    if is_always_new(nlsolver) || (iter == 1 && new_W)
        linres = dolinsolve(integrator, linsolve; A = W, b = _vec(b), linu = _vec(dz),
                            reltol = reltol)
    else
        linres = dolinsolve(integrator, linsolve; A = nothing, b = _vec(b), linu = _vec(dz),
                            reltol = reltol)
    end

    cache.linsolve = linres.cache

    if DiffEqBase.has_destats(integrator)
        integrator.destats.nsolve += 1
    end

    # relaxed Newton
    # Diagonally Implicit Runge-Kutta Methods for Ordinary Differential
    # Equations. A Review, by Christopher A. Kennedy and Mark H. Carpenter
    # page 54.
    if isdae
        γdt = α * invγdt
    else
        γdt = γ * dt
    end

    !(W_γdt ≈ γdt) && (rmul!(dz, 2 / (1 + γdt / W_γdt)))
    relax!(dz, nlsolver, integrator, f)

    calculate_residuals!(atmp, dz, uprev, ustep, opts.abstol, opts.reltol,
                         opts.internalnorm, t)
    ndz = opts.internalnorm(atmp, t)
    # NDF and BDF are special because the truncation error is directly
    # proportional to the total displacement.
    if integrator.alg isa QNDF
        ndz *= error_constant(integrator, alg_order(integrator.alg))
    end

    # compute next iterate
    @.. broadcast=false ztmp=z - dz

    ndz
end

@muladd function compute_step!(nlsolver::NLSolver{<:NLNewton, true, <:Array}, integrator)
    @unpack uprev, t, p, dt, opts = integrator
    @unpack z, tmp, ztmp, γ, α, iter, cache = nlsolver
    @unpack W_γdt, ustep, tstep, k, atmp, dz, W, new_W, invγdt, linsolve, weight = cache
    f = nlsolve_f(integrator)
    isdae = f isa DAEFunction

    if DiffEqBase.has_destats(integrator)
        integrator.destats.nf += 1
    end

    b, ustep = _compute_rhs!(nlsolver, integrator, f, z)

    # update W
    if W isa DiffEqBase.AbstractDiffEqLinearOperator
        update_coefficients!(W, ustep, p, tstep)
    end

    if integrator.opts.adaptive
        reltol = integrator.opts.reltol
    else
        reltol = eps(eltype(dz))
    end

    if is_always_new(nlsolver) || (iter == 1 && new_W)
        linres = dolinsolve(integrator, linsolve; A = W, b = _vec(b), linu = _vec(dz),
                            reltol = reltol)
    else
        linres = dolinsolve(integrator, linsolve; A = nothing, b = _vec(b), linu = _vec(dz),
                            reltol = reltol)
    end

    cache.linsolve = linres.cache

    if DiffEqBase.has_destats(integrator)
        integrator.destats.nsolve += 1
    end

    # relaxed Newton
    # Diagonally Implicit Runge-Kutta Methods for Ordinary Differential
    # Equations. A Review, by Christopher A. Kennedy and Mark H. Carpenter
    # page 54.
    if isdae
        γdt = α * invγdt
    else
        γdt = γ * dt
    end

    !(W_γdt ≈ γdt) && (rmul!(dz, 2 / (1 + γdt / W_γdt)))
    relax!(dz, nlsolver, integrator, f)

    calculate_residuals!(atmp, dz, uprev, ustep, opts.abstol, opts.reltol,
                         opts.internalnorm, t)
    ndz = opts.internalnorm(atmp, t)
    # NDF and BDF are special because the truncation error is directly
    # proportional to the total displacement.
    if integrator.alg isa QNDF
        ndz *= error_constant(integrator, alg_order(integrator.alg))
    end

    # compute next iterate
    @inbounds @simd ivdep for i in eachindex(z)
        ztmp[i] = z[i] - dz[i]
    end

    ndz
end

@inline function _compute_rhs(nlsolver::NLSolver{<:NLNewton, false}, integrator, f::TF, z) where {TF}
    @unpack uprev, t, p, dt = integrator
    @unpack tmp, ztmp, γ, α, cache = nlsolver
    @unpack tstep, invγdt = cache
    isdae = TF <: DAEFunction
    if isdae
        # not all predictors are uprev, for other forms of predictors, defined in u₀
        if isdefined(integrator.cache, :u₀)
            ustep = @.. broadcast=false integrator.cache.u₀+z
        else
            ustep = @.. broadcast=false uprev+z
        end
        dustep = @. (tmp + α * z) * invγdt
        ztmp = f(dustep, ustep, p, t)
    else
        mass_matrix = integrator.f.mass_matrix
        if nlsolver.method === COEFFICIENT_MULTISTEP
            ustep = z
            # tmp = outertmp ./ hγ
            if mass_matrix === I
                ztmp = tmp .+ f(z, p, tstep) .- (α * invγdt) .* z
            else
                update_coefficients!(mass_matrix, ustep, p, tstep)
                ztmp = tmp .+ f(z, p, tstep) .- (mass_matrix * z) .* (α * invγdt)
            end
        else
            ustep = @. tmp + γ * z
            if mass_matrix === I
                ztmp = (dt .* f(ustep, p, tstep) .- z) .* invγdt
            else
                update_coefficients!(mass_matrix, ustep, p, tstep)
                ztmp = (dt .* f(ustep, p, tstep) .- mass_matrix * z) .* invγdt
            end
        end
    end
    return ztmp, ustep
end

@inline function _compute_rhs!(nlsolver::NLSolver{<:NLNewton, true}, integrator, f::TF, z) where {TF}
    @unpack uprev, t, p, dt = integrator
    @unpack tmp, ztmp, γ, α, cache = nlsolver
    @unpack ustep, tstep, k, invγdt = cache
    isdae = TF <: DAEFunction
    b = if isdae
        @.. broadcast=false ztmp=(tmp + α * z) * invγdt
        # not all predictors are uprev, for other forms of predictors, defined in u₀
        if isdefined(integrator.cache, :u₀)
            @.. broadcast=false ustep=integrator.cache.u₀ + z
        else
            @.. broadcast=false ustep=uprev + z
        end
        f(k, ztmp, ustep, p, tstep)
        _vec(k)
    else
        mass_matrix = integrator.f.mass_matrix
        if nlsolver.method === COEFFICIENT_MULTISTEP
            ustep = z
            f(k, z, p, tstep)
            if mass_matrix === I
                @.. broadcast=false ztmp=tmp + k - (α * invγdt) * z
            else
                update_coefficients!(mass_matrix, ustep, p, tstep)
                mul!(_vec(ztmp), mass_matrix, _vec(z))
                @.. broadcast=false ztmp=tmp + k - (α * invγdt) * ztmp
            end
        else
            @.. broadcast=false ustep=tmp + γ * z
            f(k, ustep, p, tstep)
            if mass_matrix === I
                @.. broadcast=false ztmp=(dt * k - z) * invγdt
            else
                update_coefficients!(mass_matrix, ustep, p, tstep)
                mul!(_vec(ztmp), mass_matrix, _vec(z))
                @.. broadcast=false ztmp=(dt * k - ztmp) * invγdt
            end
        end
        _vec(ztmp)
    end
end

@inline function _compute_rhs!(nlsolver::NLSolver{<:NLNewton, true, <:Array}, integrator, f::TF, z) where {TF}
    @unpack uprev, t, p, dt = integrator
    @unpack tmp, ztmp, γ, α, cache = nlsolver
    @unpack ustep, tstep, k, invγdt = cache
    isdae = TF <: DAEFunction
    b = if isdae
        @inbounds @simd ivdep for i in eachindex(z)
            ztmp[i] = (tmp[i] + α * z[i]) * invγdt
        end
        if isdefined(integrator.cache, :u₀)
            @inbounds @simd ivdep for i in eachindex(z)
                ustep[i] = integrator.cache.u₀[i] + z[i]
            end
            #@.. broadcast=false ustep = integrator.cache.u₀ + z
        else
            @inbounds @simd ivdep for i in eachindex(z)
                ustep[i] = uprev[i] + z[i]
            end
        end
        f(k, ztmp, ustep, p, tstep)
        _vec(k)
    else
        mass_matrix = integrator.f.mass_matrix
        if nlsolver.method === COEFFICIENT_MULTISTEP
            ustep = z
            f(k, z, p, tstep)
            if mass_matrix === I
                @inbounds @simd ivdep for i in eachindex(z)
                    ztmp[i] = tmp[i] + k[i] - (α * invγdt) * z[i]
                end
            else
                update_coefficients!(mass_matrix, ustep, p, tstep)
                mul!(_vec(ztmp), mass_matrix, _vec(z))

                @inbounds @simd ivdep for i in eachindex(z)
                    ztmp[i] = tmp[i] + k[i] - (α * invγdt) * ztmp[i]
                end
            end
        else
            @inbounds @simd ivdep for i in eachindex(z)
                ustep[i] = tmp[i] + γ * z[i]
            end
            f(k, ustep, p, tstep)
            if mass_matrix === I
                @inbounds @simd ivdep for i in eachindex(z)
                    ztmp[i] = (dt * k[i] - z[i]) * invγdt
                end
            else
                update_coefficients!(mass_matrix, ustep, p, tstep)
                mul!(_vec(ztmp), mass_matrix, _vec(z))
                @inbounds @simd ivdep for i in eachindex(z)
                    ztmp[i] = (dt * k[i] - ztmp[i]) * invγdt
                end
            end
        end
        _vec(ztmp)
    end
    return b, ustep
end

## relax!
relax!(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF) where {TF} = relax!(dz, nlsolver, integrator, f, relax(nlsolver))
relax(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF) where {TF} = relax(dz, nlsolver, integrator, f, relax(nlsolver))
function relax!(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF, r::Number) where {TF}
    if !iszero(r)
        rmul!(dz, 1 - r)
    end
end

function relax!(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF, linesearch) where {TF}
    let dz=dz,
        integrator=integrator,
        nlsolver=nlsolver,
        f=f,
        linesearch=linesearch;
        @unpack uprev, t, p, dt, opts = integrator
        @unpack z, tmp, ztmp, γ, iter, cache = nlsolver
        @unpack ustep, atmp = cache
        function resid(z)
            # recompute residual (rhs)
            b, ustep = _compute_rhs!(nlsolver, integrator, f, z)
            calculate_residuals!(atmp, b, uprev, ustep, opts.abstol, opts.reltol, opts.internalnorm, t)
            ndz = opts.internalnorm(atmp, t)
            return ndz
        end
        function ϕ(α)
            z = @.. atmp = nlsolver.z - dz * α
            return resid(z)
        end
        function dϕ(α)
            ϵ = sqrt(eps())
            return (ϕ(α+ϵ) - ϕ(α)) / ϵ
        end
        function ϕdϕ(α)
            ϵ = sqrt(eps())
            ϕ_1 = ϕ(α)
            ϕ_2 = ϕ(α + ϵ)
            dϕ = (ϕ_2 - ϕ_1) / ϵ
            return ϕ_1, dϕ
        end
        α0 = one(eltype(dz))
        ϕ0, dϕ0 = ϕdϕ(zero(α0))
        α, _ = linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
        @.. dz = dz * α
        return dz
    end
end

function relax(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF, r::Number) where {TF}
    if !iszero(r)
        dz = (1 - r)*dz
    end
    return dz
end

function relax(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF, linesearch) where {TF}
    let dz=dz,
        integrator=integrator,
        nlsolver=nlsolver,
        f=f,
        linesearch=linesearch;
        @unpack uprev, t, p, dt, opts = integrator
        @unpack z, tmp, ztmp, γ, iter, cache = nlsolver
        function resid(z)
            # recompute residual (rhs)
            b, ustep = _compute_rhs(nlsolver, integrator, f, z)
            atmp = calculate_residuals(b, uprev, ustep, opts.abstol, opts.reltol, opts.internalnorm, t)
            ndz = opts.internalnorm(atmp, t)
            return ndz
        end
        function ϕ(α)
            z = @.. nlsolver.z - dz * α
            return resid(z)
        end
        function dϕ(α)
            ϵ = sqrt(eps())
            return (ϕ(α+ϵ) - ϕ(α)) / ϵ
        end
        function ϕdϕ(α)
            ϵ = sqrt(eps())
            ϕ_1 = ϕ(α)
            ϕ_2 = ϕ(α + ϵ)
            dϕ = (ϕ_2 - ϕ_1) / ϵ
            return ϕ_1, dϕ
        end
        α0 = one(eltype(dz))
        ϕ0, dϕ0 = ϕdϕ(zero(α0))
        α, _ = linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
        dz = dz * α
        return dz
    end
end

## resize!

function Base.resize!(nlcache::NLNewtonCache, ::AbstractNLSolver, integrator, i::Int)
    resize!(nlcache.ustep, i)
    resize!(nlcache.k, i)
    resize!(nlcache.atmp, i)
    resize!(nlcache.dz, i)
    resize!(nlcache.du1, i)
    if nlcache.jac_config !== nothing
        resize_jac_config!(nlcache.jac_config, i)
    end
    resize!(nlcache.weight, i)

    # resize J and W (or rather create new ones of appropriate size and type)
    resize_J_W!(nlcache, integrator, i)

    nothing
end
