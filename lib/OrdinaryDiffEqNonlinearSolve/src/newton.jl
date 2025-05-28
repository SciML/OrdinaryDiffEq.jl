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

function initialize!(nlsolver::NLSolver{<:NonlinearSolveAlg, false},
        integrator::DiffEqBase.DEIntegrator)
    @unpack uprev, t, p, dt, opts, f = integrator
    @unpack z, tmp, ztmp, γ, α, iter, cache, method, alg = nlsolver
    cache.invγdt = inv(dt * nlsolver.γ)
    cache.tstep = integrator.t + nlsolver.c * dt

    @unpack ustep, tstep, k, invγdt = cache
    if DiffEqBase.has_stats(integrator)
        integrator.stats.nf += cache.cache.stats.nf
        integrator.stats.nnonliniter += cache.cache.stats.nsteps
        integrator.stats.njacs += cache.cache.stats.njacs
    end
    if f isa DAEFunction
        nlp_params = (tmp, α, tstep, invγdt, p, dt, uprev, f)
    else
        nlp_params = (tmp, γ, α, tstep, invγdt, method, p, dt, f)
    end
    new_prob = remake(cache.prob, p = nlp_params, u0 = z)
    cache.cache = init(new_prob, alg.alg)
    nothing
end

function initialize!(nlsolver::NLSolver{<:NonlinearSolveAlg, true},
        integrator::DiffEqBase.DEIntegrator)
    @unpack uprev, t, p, dt, opts, f = integrator
    @unpack z, tmp, ztmp, γ, α, iter, cache, method, alg = nlsolver

    cache.invγdt = inv(dt * nlsolver.γ)
    cache.tstep = integrator.t + nlsolver.c * dt

    @unpack ustep, tstep, k, invγdt = cache

    if DiffEqBase.has_stats(integrator)
        integrator.stats.nf += cache.cache.stats.nf
        integrator.stats.nnonliniter += cache.cache.stats.nsteps
        integrator.stats.njacs += cache.cache.stats.njacs
    end
    if f isa DAEFunction
        nlp_params = (tmp, ztmp, ustep, γ, α, tstep, k, invγdt, p, dt, f)
    else
        nlp_params = (tmp, ustep, γ, α, tstep, k, invγdt, method, p, dt, f)
    end
    new_prob = remake(cache.prob, p = nlp_params, u0 = z)
    cache.cache = init(new_prob, alg.alg)
    nothing
end

## compute_step!

@muladd function compute_step!(nlsolver::NLSolver{<:NonlinearSolveAlg, false}, integrator)
    @unpack uprev, t, p, dt, opts = integrator
    @unpack z, tmp, ztmp, γ, α, cache, method = nlsolver
    @unpack tstep, invγdt = cache

    nlcache = nlsolver.cache.cache

    if is_always_new(nlsolver) || new_jac || new_W
        recompute_jacobian = true
    else
        recompute_jacobian = false
    end

    step!(nlcache; recompute_jacobian)
    nlsolver.ztmp = nlcache.u

    ustep = compute_ustep(tmp, γ, z, method)
    atmp = calculate_residuals(nlcache.fu, uprev, ustep, opts.abstol, opts.reltol,
        opts.internalnorm, t)
    ndz = opts.internalnorm(atmp, t)
    #ndz = opts.internalnorm(nlcache.fu, t)
    # NDF and BDF are special because the truncation error is directly
    # proportional to the total displacement.
    if has_special_newton_error(integrator.alg)
        ndz *= error_constant(integrator, alg_order(integrator.alg))
    end
    return ndz
end

@muladd function compute_step!(nlsolver::NLSolver{<:NonlinearSolveAlg, true}, integrator)
    @unpack uprev, t, p, dt, opts = integrator
    @unpack z, tmp, ztmp, γ, α, cache, method = nlsolver
    @unpack tstep, invγdt, atmp, ustep = cache

    nlcache = nlsolver.cache.cache
    new_jac, new_W = do_newJW(integrator, integrator.alg, nlsolver, false)

    if is_always_new(nlsolver) || new_jac || new_W
        recompute_jacobian = true
    else
        recompute_jacobian = false
    end

    step!(nlcache; recompute_jacobian)
    @.. broadcast=false ztmp=nlcache.u

    ustep = compute_ustep!(ustep, tmp, γ, z, method)
    calculate_residuals!(atmp, nlcache.fu, uprev, ustep, opts.abstol, opts.reltol,
        opts.internalnorm, t)
    ndz = opts.internalnorm(atmp, t)
    #ndz = opts.internalnorm(nlcache.fu, t)
    # NDF and BDF are special because the truncation error is directly
    # proportional to the total displacement.
    if has_special_newton_error(integrator.alg)
        ndz *= error_constant(integrator, alg_order(integrator.alg))
    end
    ndz
end

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
@muladd function compute_step!(nlsolver::NLSolver{<:NLNewton, false}, integrator, γW)
    @unpack uprev, t, p, dt, opts = integrator
    @unpack z, tmp, ztmp, γ, α, cache, method = nlsolver
    @unpack tstep, W, invγdt = cache

    f = nlsolve_f(integrator)

    if f isa DAEFunction
        _uprev = get_dae_uprev(integrator, uprev)
        ztmp, ustep = _compute_rhs(tmp, α, tstep, invγdt, p, _uprev, f, z)
    else
        ztmp, ustep = _compute_rhs(tmp, γ, α, tstep, invγdt, method, p, dt, f, z)
    end

    if DiffEqBase.has_stats(integrator)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    # update W
    if W isa Union{WOperator, StaticWOperator}
        update_coefficients!(W, ustep, p, tstep)
    elseif W isa AbstractSciMLOperator
        error("Non-concrete Jacobian not yet supported by out-of-place Newton solve.")
    end
    dz = _reshape(W \ _vec(ztmp), axes(ztmp))
    dz = relax(dz, nlsolver, integrator, f)
    if DiffEqBase.has_stats(integrator)
        integrator.stats.nsolve += 1
    end

    atmp = calculate_residuals(dz, uprev, ustep, opts.abstol, opts.reltol,
        opts.internalnorm, t)
    ndz = opts.internalnorm(atmp, t)
    # NDF and BDF are special because the truncation error is directly
    # proportional to the total displacement.
    if has_special_newton_error(integrator.alg)
        ndz *= error_constant(integrator, alg_order(integrator.alg))
    end

    # compute next iterate
    nlsolver.ztmp = z .- dz

    ndz
end

@muladd function compute_step!(nlsolver::NLSolver{<:NLNewton, true}, integrator, γW)
    @unpack uprev, t, p, dt, opts = integrator
    @unpack z, tmp, ztmp, γ, α, iter, cache, method = nlsolver
    @unpack W_γdt, ustep, tstep, k, atmp, dz, W, new_W, invγdt, linsolve, weight = cache

    f = nlsolve_f(integrator)
    isdae = f isa DAEFunction

    if DiffEqBase.has_stats(integrator)
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end

    if isdae
        _uprev = get_dae_uprev(integrator, uprev)
        b, ustep = _compute_rhs!(tmp, ztmp, ustep, α, tstep, k, invγdt, p, _uprev, f, z)
    else
        b, ustep = _compute_rhs!(
            tmp, ztmp, ustep, γ, α, tstep, k, invγdt, method, p, dt, f, z)
    end

    # update W
    if W isa Union{WOperator, StaticWOperator}
        update_coefficients!(W, ustep, p, tstep)
    elseif W isa AbstractSciMLOperator
        # logic for generic AbstractSciMLOperator does not yet support partial state updates, so provide full state
        update_coefficients!(W, ustep, p, tstep; dtgamma = γW, transform = true)
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
        linres = dolinsolve(
            integrator, linsolve; A = nothing, b = _vec(b), linu = _vec(dz),
            reltol = reltol)
    end

    if !SciMLBase.successful_retcode(linres.retcode) &&
       linres.retcode != SciMLBase.ReturnCode.Default
        return convert(eltype(atmp,), Inf)
    end

    cache.linsolve = linres.cache

    if DiffEqBase.has_stats(integrator)
        integrator.stats.nsolve += 1
    end

    # relaxed Newton
    # Diagonally Implicit Runge-Kutta Methods for Ordinary Differential
    # Equations. A Review, by Christopher A. Kennedy and Mark H. Carpenter
    # page 54.
    γdt = isdae ? α * invγdt : γ * dt

    !(W_γdt ≈ γdt) && (rmul!(dz, 2 / (1 + γdt / W_γdt)))
    relax!(dz, nlsolver, integrator, f)

    calculate_residuals!(atmp, dz, uprev, ustep, opts.abstol, opts.reltol,
        opts.internalnorm, t)
    ndz = opts.internalnorm(atmp, t)

    # NDF and BDF are special because the truncation error is directly
    # proportional to the total displacement.
    if has_special_newton_error(integrator.alg)
        ndz *= error_constant(integrator, alg_order(integrator.alg))
    end

    # compute next iterate
    @.. broadcast=false ztmp=z - dz

    ndz
end

function get_dae_uprev(integrator, uprev)
    # not all predictors are uprev, for other forms of predictors, defined in u₀
    if isdefined(integrator.cache, :u₀)
        integrator.cache.u₀
    else
        uprev
    end
end

function _compute_rhs(tmp, α, tstep, invγdt, p, uprev, f::TF, z) where {TF <: DAEFunction}
    ustep = @.. uprev + z
    dustep = @. (tmp + α * z) * invγdt
    ztmp = f(dustep, ustep, p, tstep)
    return ztmp, ustep
end

function compute_ustep(tmp, γ, z, method)
    if method === COEFFICIENT_MULTISTEP
        z
    else
        @. tmp + γ * z
    end
end

function compute_ustep!(ustep, tmp, γ, z, method)
    if method === COEFFICIENT_MULTISTEP
        ustep = z
    else
        @.. ustep = tmp + γ * z
    end
    ustep
end

function _compute_rhs(tmp, γ, α, tstep, invγdt, method::MethodType, p, dt, f, z)
    mass_matrix = f.mass_matrix
    ustep = compute_ustep(tmp, γ, z, method)
    if method === COEFFICIENT_MULTISTEP
        # tmp = outertmp ./ hγ
        if mass_matrix === I
            ztmp = tmp .+ f(z, p, tstep) .- (α * invγdt) .* z
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            ztmp = tmp .+ f(z, p, tstep) .- (mass_matrix * z) .* (α * invγdt)
        end
    else
        if mass_matrix === I
            ztmp = (dt * f(ustep, p, tstep) - z) * invγdt
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            ztmp = (dt * f(ustep, p, tstep) - mass_matrix * z) * invγdt
        end
    end
    return ztmp, ustep
end

function _compute_rhs!(tmp, ztmp, ustep, α, tstep, k,
        invγdt, p, uprev, f::TF, z) where {TF <: DAEFunction}
    @.. broadcast=false ztmp=(tmp + α * z) * invγdt
    @.. ustep = uprev + z
    f(k, ztmp, ustep, p, tstep)
    return _vec(k), ustep
end

function _compute_rhs!(tmp, ztmp, ustep, γ, α, tstep, k,
        invγdt, method::MethodType, p, dt, f, z)
    mass_matrix = f.mass_matrix
    ustep = compute_ustep!(ustep, tmp, γ, z, method)
    if method === COEFFICIENT_MULTISTEP
        f(k, z, p, tstep)
        if mass_matrix === I
            @.. broadcast=false ztmp=tmp + k - (α * invγdt) * z
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            mul!(_vec(ztmp), mass_matrix, _vec(z))
            @.. broadcast=false ztmp=tmp + k - (α * invγdt) * ztmp
        end
    else
        f(k, ustep, p, tstep)
        if mass_matrix === I
            @.. ztmp = (dt * k - z) * invγdt
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            mul!(_vec(ztmp), mass_matrix, _vec(z))
            @.. ztmp = (dt * k - ztmp) * invγdt
        end
    end
    return _vec(ztmp), ustep
end

function _compute_rhs!(tmp::Array, ztmp::Array, ustep::Array, α, tstep, k,
        invγdt, p, uprev, f::TF, z) where {TF <: DAEFunction}
    @inbounds @simd ivdep for i in eachindex(z)
        ztmp[i] = (tmp[i] + α * z[i]) * invγdt
    end
    @inbounds @simd ivdep for i in eachindex(z)
        ustep[i] = uprev[i] + z[i]
    end
    f(k, ztmp, ustep, p, tstep)

    return _vec(k), ustep
end

function _compute_rhs!(tmp::Array, ztmp::Array, ustep::Array, γ, α, tstep, k,
        invγdt, method::MethodType, p, dt, f, z)
    mass_matrix = f.mass_matrix
    ustep = compute_ustep!(ustep, tmp, γ, z, method)
    if method === COEFFICIENT_MULTISTEP
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

    return _vec(ztmp), ustep
end

## relax!
function relax!(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF) where {TF}
    relax!(dz, nlsolver, integrator, f, relax(nlsolver))
end
function relax(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF) where {TF}
    relax(dz, nlsolver, integrator, f, relax(nlsolver))
end
function relax!(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        r::Nothing) where {TF}
        dz
end
function relax!(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        r::Number) where {TF}
    if !iszero(r)
        rmul!(dz, 1 - r)
    end
    dz
end

function relax!(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        linesearch) where {TF}
    let dz = dz,
        integrator = integrator,
        nlsolver = nlsolver,
        f = f,
        linesearch = linesearch

        @unpack uprev, t, p, dt, opts, isdae = integrator
        @unpack z, tmp, ztmp, γ, iter, α, cache, method = nlsolver
        @unpack ustep, atmp, tstep, k, invγdt = cache
        function resid(z)
            # recompute residual (rhs)
            if isdae
                _uprev = get_dae_uprev(integrator, uprev)
                b, ustep2 = _compute_rhs!(
                    tmp, ztmp, ustep, α, tstep, k, invγdt, p, _uprev, f::TF, z)
            else
                b, ustep2 = _compute_rhs!(
                    tmp, ztmp, ustep, γ, α, tstep, k, invγdt, method, p, dt, f, z)
            end
            calculate_residuals!(atmp, b, uprev, ustep2, opts.abstol, opts.reltol,
                opts.internalnorm, t)
            ndz = opts.internalnorm(atmp, t)
            return ndz
        end
        function ϕ(α)
            local z = @.. atmp = nlsolver.z - dz * α
            res = resid(z)
            return res
        end
        function dϕ(α)
            ϵ = sqrt(eps())
            return (ϕ(α + ϵ) - ϕ(α)) / ϵ
        end
        function ϕdϕ(α)
            ϵ = sqrt(eps())
            ϕ_1 = ϕ(α)
            ϕ_2 = ϕ(α + ϵ)
            ∂ϕ∂α = (ϕ_2 - ϕ_1) / ϵ
            return ϕ_1, ∂ϕ∂α
        end
        α0 = one(eltype(ustep))
        ϕ0, dϕ0 = ϕdϕ(zero(α0))
        α, _ = linesearch(ϕ, dϕ, ϕdϕ, α0, ϕ0, dϕ0)
        @.. dz = dz * α
        return dz
    end
end

function relax(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        r::Number) where {TF}
    if !iszero(r)
        dz = (1 - r) * dz
    end
    return dz
end

function relax(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        r::Nothing) where {TF}
    return dz
end

function relax(dz, nlsolver::AbstractNLSolver, integrator::DEIntegrator, f::TF,
        linesearch) where {TF}
    let dz = dz,
        integrator = integrator,
        nlsolver = nlsolver,
        f = f,
        linesearch = linesearch

        @unpack uprev, t, p, dt, opts = integrator
        @unpack z, tmp, ztmp, γ, iter, cache, method = nlsolver
        @unpack ustep, atmp, tstep, k, invγdt = cache
        function resid(z)
            # recompute residual (rhs)
            if f isa DAEFunction
                _uprev = get_dae_uprev(integrator, uprev)
                ztmp, ustep2 = _compute_rhs(tmp, α, tstep, invγdt, p, dt, _uprev, f, z)
            else
                ztmp, ustep2 = _compute_rhs(tmp, γ, α, tstep, invγdt, method, p, f, z)
            end
            atmp = calculate_residuals(b, uprev, ustep2, opts.abstol, opts.reltol,
                opts.internalnorm, t)
            ndz = opts.internalnorm(atmp, t)
            return ndz
        end
        function ϕ(α)
            local z = @.. nlsolver.z - dz * α
            return resid(z)
        end
        function dϕ(α)
            ϵ = sqrt(eps())
            return (ϕ(α + ϵ) - ϕ(α)) / ϵ
        end
        function ϕdϕ(α)
            ϵ = sqrt(eps())
            ϕ_1 = ϕ(α)
            ϕ_2 = ϕ(α + ϵ)
            ∂ϕ∂α = (ϕ_2 - ϕ_1) / ϵ
            return ϕ_1, ∂ϕ∂α
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

    resize_jac_config!(nlcache, integrator)
    resize!(nlcache.weight, i)

    # resize J and W (or rather create new ones of appropriate size and type)
    resize_J_W!(nlcache, integrator, i)

    nothing
end

