# Step-size homotopy for the implicit stage equations.
#
# The stage residual of `nlsolve!` (see nlsolve.jl) is embedded into the one-parameter
# family obtained by scaling the `f` evaluation with `λ ∈ [0, 1]`:
#
#   DIRK form:                H(z, λ) = λ⋅dt⋅f(tmp + γ⋅z, p, tstep) - M z
#   COEFFICIENT_MULTISTEP:    H(z, λ) = (γ dt / α)⋅(tmp + λ⋅f(z, p, tstep)) - M z
#
# `λ` is a dimensionless multiplier on the effective step size, so `λspan = (0, 1)`
# regardless of the magnitude of `dt` (a `(0, dt)` span would interact badly with the
# absolute `min_dλ`/`min_ds` floors of the continuation solvers for small `dt`). At
# `λ = 0` the DIRK-form solution is exactly `z = 0` and the multistep-form problem is
# linear, so the continuation always starts from a solvable anchor. Everything except
# the leading `λ` (history `tmp`, stage time `tstep`, the residual scaling `1/(γ dt)`)
# is frozen at the target step size: intermediate `λ` solutions are continuation
# waypoints, not method steps, which keeps the embedding valid for methods whose
# coefficients depend on `dt` (BDF step-ratio history, SDIRK tableaus).

_scalar_tol(x::Number) = x
_scalar_tol(x) = minimum(x)

function _init_homotopy_nonlinear_cache(prob, alg, abstol, verbose)
    alg.alg isa Union{HomotopySweep, KantorovichHomotopy} || return nothing
    if alg.reltol === nothing
        return init(prob, alg.alg; abstol, verbose)
    end
    return init(prob, alg.alg; abstol, reltol = alg.reltol, verbose)
end

function homotopy_odenlf(ztmp, z, p, λ)
    tmp, ustep, γ, α, tstep, k, invγdt, method, _p, dt, f, nf = p
    nf[] += 1
    return _compute_homotopy_rhs!(
        tmp, ztmp, ustep, γ, α, tstep, k, invγdt, method, _p, dt, λ, f, z
    )[1]
end

function homotopy_oopodenlf(z, p, λ)
    tmp, γ, α, tstep, invγdt, method, _p, dt, f, nf = p
    nf[] += 1
    return _compute_homotopy_rhs(tmp, γ, α, tstep, invγdt, method, _p, dt, λ, f, z)[1]
end

# Unlike `_compute_rhs`, the homotopy residuals are kept in `u` units (the `1/(γ dt)`
# scaling of the Newton residual is dropped, i.e. the DIRK form is multiplied by `γ dt`
# and the multistep form by `γ dt / α`). The continuation solver terminates on the
# residual norm with an absolute tolerance, and the `1/(γ dt)`-scaled residuals sit at
# `~|u|/(γ dt)`, far above any fixed tolerance floating point can reach for small `dt`.
function _compute_homotopy_rhs(
        tmp, γ, α, tstep, invγdt, method::MethodType, p, dt, λ, f::F, z
    ) where {F}
    mass_matrix = f.mass_matrix
    if method === COEFFICIENT_MULTISTEP
        # The continuation solver may hand `z` as a view into its augmented `(u, λ)`
        # vector; materialize it so wrapped functions (AutoSpecialize) and user code
        # always see a plain array.
        ustep = copy(z)
        s = inv(α * invγdt)
        if mass_matrix === I
            ztmp = s .* (tmp .+ λ .* f(ustep, p, tstep)) .- ustep
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            ztmp = s .* (tmp .+ λ .* f(ustep, p, tstep)) .- mass_matrix * ustep
        end
    else
        ustep = compute_ustep(tmp, γ, z, method)
        if mass_matrix === I
            ztmp = (λ * dt) * f(ustep, p, tstep) - z
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            ztmp = (λ * dt) * f(ustep, p, tstep) - mass_matrix * z
        end
    end
    return ztmp, ustep
end

function _compute_homotopy_rhs!(
        tmp, ztmp, ustep, γ, α, tstep, k, invγdt, method::MethodType, p, dt, λ, f, z
    )
    mass_matrix = f.mass_matrix
    if method === COEFFICIENT_MULTISTEP
        # See the out-of-place method: `z` can be a view into the continuation
        # solver's augmented vector, so evaluate `f` on the dense `ustep` buffer.
        copyto!(ustep, z)
        f(k, ustep, p, tstep)
        s = inv(α * invγdt)
        if mass_matrix === I
            @.. broadcast = false ztmp = s * (tmp + λ * k) - z
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            mul!(_vec(ztmp), mass_matrix, _vec(ustep))
            @.. broadcast = false ztmp = s * (tmp + λ * k) - ztmp
        end
    else
        ustep = compute_ustep!(ustep, tmp, γ, z, method)
        f(k, ustep, p, tstep)
        if mass_matrix === I
            @.. broadcast = false ztmp = (λ * dt) * k - z
        else
            update_coefficients!(mass_matrix, ustep, p, tstep)
            mul!(_vec(ztmp), mass_matrix, _vec(z))
            @.. broadcast = false ztmp = (λ * dt) * k - ztmp
        end
    end
    return _vec(ztmp), ustep
end

# One full continuation solve replaces the shared Newton iteration loop: the generic
# `nlsolve!` drives `compute_step!` iteratively with η-convergence control, which does
# not fit a solver whose single `solve` already runs an entire predictor-corrector
# sweep, so this method bypasses it.
function nlsolve!(
        nlsolver::NLSolver{<:HomotopyNonlinearSolveAlg, iip},
        integrator::SciMLBase.DEIntegrator,
        cache = nothing, repeat_step = false
    ) where {iip}
    (; t, p, dt) = integrator
    (; z, tmp, ztmp, γ, α, method, alg) = nlsolver
    nlcache = nlsolver.cache
    f = nlsolve_f(integrator)

    nlcache.invγdt = inv(dt * γ)
    nlcache.tstep = t + nlsolver.c * dt
    (; tstep, invγdt, nlfunc) = nlcache

    nlcache.nf[] = 0
    λspan = (zero(γ), one(γ))
    if iip
        (; ustep, k) = nlcache
        # For the DIRK form the λ = 0 anchor solution is exactly `z = 0`; the multistep
        # form is linear at λ = 0 but not trivial, so keep the method's predictor there.
        u0 = method === COEFFICIENT_MULTISTEP ? copyto!(ztmp, z) : fill!(ztmp, false)
        nlp_params = (tmp, ustep, γ, α, tstep, k, invγdt, method, p, dt, f, nlcache.nf)
    else
        u0 = method === COEFFICIENT_MULTISTEP ? z : zero(z)
        nlp_params = (tmp, γ, α, tstep, invγdt, method, p, dt, f, nlcache.nf)
    end
    # The residual is in `u` units, so the integrator's absolute tolerance sets its
    # natural convergence scale; κ tightens it like the Newton κ⋅tol criterion.
    abstol = alg.abstol === nothing ?
        nlsolver.κ * _scalar_tol(integrator.opts.abstol) : alg.abstol
    kwargs = (; abstol = abstol)
    alg.reltol === nothing || (kwargs = merge(kwargs, (; reltol = alg.reltol)))
    if nlcache.needs_rebuild
        prob = SciMLBase.HomotopyProblem{iip}(nlfunc, u0, nlp_params; λspan = λspan)
        nlcache.continuation_cache = _init_homotopy_nonlinear_cache(
            prob, alg, abstol, nlcache.verbose
        )
        nlcache.needs_rebuild = false
    end
    sol = if nlcache.continuation_cache === nothing
        prob = SciMLBase.HomotopyProblem{iip}(nlfunc, u0, nlp_params; λspan = λspan)
        solve(prob, alg.alg; kwargs...)
    else
        SciMLBase.reinit!(nlcache.continuation_cache, u0; p = nlp_params, kwargs...)
        solve!(nlcache.continuation_cache)
    end

    if iip
        copyto!(z, sol.u)
    else
        nlsolver.z = sol.u
    end

    if SciMLBase.successful_retcode(sol.retcode)
        nlsolver.status = Convergence
        nlsolver.nfails = 0
    else
        @SciMLMessage(
            lazy"Homotopy continuation of the stage equation failed to reach λ = 1 (retcode = $(sol.retcode)); rejecting the step",
            integrator.opts.verbose, :newton_convergence
        )
        nlsolver.status = Divergence
        nlsolver.nfails += 1
    end
    nlsolver.iter = 1
    # The continuation solutions do not carry stats, so `f` calls are counted by the
    # residual itself; inner Jacobian/linear-solve counts are not tracked.
    if SciMLBase.has_stats(integrator)
        integrator.stats.nf += nlcache.nf[]
    end

    return postamble!(nlsolver, integrator)
end

function du_cache(nlcache::HomotopyNonlinearSolveCache)
    return nlcache.k === nothing ? nothing : (nlcache.k,)
end

function Base.resize!(
        nlcache::HomotopyNonlinearSolveCache, ::AbstractNLSolver, integrator, i::Int
    )
    nlcache.ustep === nothing || resize!(nlcache.ustep, i)
    nlcache.k === nothing || resize!(nlcache.k, i)
    nlcache.needs_rebuild = nlcache.continuation_cache !== nothing
    return nothing
end
