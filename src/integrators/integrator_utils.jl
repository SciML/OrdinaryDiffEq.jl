@inline function DiffEqNoiseProcess.setup_next_step!(integrator::SDEIntegrator)
    !isnothing(integrator.W) &&
        DiffEqNoiseProcess.setup_next_step!(integrator.W, integrator.u, integrator.p)
    return !isnothing(integrator.P) &&
        DiffEqNoiseProcess.setup_next_step!(integrator.P, integrator.u, integrator.p)
end

@inline function DiffEqNoiseProcess.reject_step!(integrator::SDEIntegrator, dtnew = integrator.dt)
    !isnothing(integrator.W) &&
        reject_step!(integrator.W, dtnew, integrator.u, integrator.p)
    return !isnothing(integrator.P) &&
        reject_step!(integrator.P, dtnew, integrator.u, integrator.p)
end

@inline function DiffEqNoiseProcess.accept_step!(integrator::SDEIntegrator, setup)
    !isnothing(integrator.W) &&
        accept_step!(integrator.W, integrator.dt, integrator.u, integrator.p, setup)
    return !isnothing(integrator.P) &&
        accept_step!(integrator.P, integrator.dt, integrator.u, integrator.p, setup)
end

@inline function DiffEqNoiseProcess.save_noise!(integrator::SDEIntegrator)
    !isnothing(integrator.W) && DiffEqNoiseProcess.save_noise!(integrator.W)
    return !isnothing(integrator.P) && DiffEqNoiseProcess.save_noise!(integrator.P)
end

# SDE's loopfooter! dispatches to ODE's shared _loopfooter!.
@inline loopfooter!(integrator::SDEIntegrator) = _loopfooter!(integrator)

function last_step_failed(integrator::SDEIntegrator)
    return integrator.last_stepfail && !integrator.opts.adaptive
end

# Use OrdinaryDiffEqCore's _savevalues! via hooks (interp_at_saveat, is_composite_algorithm, post_savevalues!)
@inline function DiffEqBase.savevalues!(
        integrator::SDEIntegrator, force_save = false, reduce_size = true,
    )::Tuple{Bool, Bool}
    return _savevalues!(integrator, force_save, reduce_size)
end

# solution_endpoint_match_cur_integrator! is now imported from OrdinaryDiffEqCore.
# SDE-specific behavior (noise acceptance, noise saving) handled via
# finalize_endpoint! and finalize_solution_storage! hooks.

# Use OrdinaryDiffEqCore's _postamble! via hooks (finalize_solution_storage!, finalize_endpoint!)
@inline function DiffEqBase.postamble!(integrator::SDEIntegrator)
    return _postamble!(integrator)
end

# handle_callbacks! is now imported from OrdinaryDiffEqCore.
# SDE-specific behavior handled via on_callbacks_complete! hook.

@inline function handle_callback_modifiers!(integrator::SDEIntegrator)
    #integrator.reeval_fsal = true
    return if integrator.P !== nothing && integrator.opts.adaptive
        if integrator.cache isa StochasticDiffEqMutableCache
            oldrate = integrator.P.cache.currate
            integrator.P.cache.rate(oldrate, integrator.u, integrator.p, integrator.t)
        else
            integrator.P.cache.currate = integrator.P.cache.rate(integrator.u, integrator.p, integrator.t)
        end
    end
end

# handle_tstop! is now imported from OrdinaryDiffEqCore.
# SDE benefits from ODE's while-loop fix for multiple matching tstops.

@inline initialize!(integrator, cache::StochasticDiffEqCache, f = integrator.f) = nothing

function nlsolve!(integrator, cache)
    return DiffEqBase.nlsolve!(cache.nlsolver, cache.nlsolver.cache, integrator)
end

function OrdinaryDiffEqCore.nlsolve_f(f, alg::StochasticDiffEqAlgorithm)
    return f isa SplitSDEFunction && issplit(alg) ? f.f1 : f
end
function OrdinaryDiffEqCore.nlsolve_f(integrator::SDEIntegrator)
    return nlsolve_f(integrator.f, unwrap_alg(integrator, true))
end

# TauLeapingDrift: wrapper for tau-leaping drift function used by nlsolver
# Computes drift(u, p, t) = c(u, p, t, rate(u, p, t), nothing)
# where c is the stoichiometry function and rate is the propensity function
struct TauLeapingDrift{C, R, RateCache, IIP}
    c::C              # Stoichiometry function (from integrator.c)
    rate::R           # Rate function (from integrator.P.cache.rate)
    rate_cache::RateCache  # Cache for rate values (for in-place version)
end

# Out-of-place: drift(u, p, t)
function (td::TauLeapingDrift{C, R, Nothing, false})(u, p, t) where {C, R}
    rates = td.rate(u, p, t)
    return td.c(u, p, t, rates, nothing)
end

# In-place: drift(du, u, p, t)
function (td::TauLeapingDrift{C, R, RateCache, true})(du, u, p, t) where {C, R, RateCache}
    td.rate(td.rate_cache, u, p, t)
    td.c(du, u, p, t, td.rate_cache, nothing)
    return nothing
end

# nlsolve_f override for ImplicitTauLeaping: return TauLeapingDrift wrapper
function OrdinaryDiffEqCore.nlsolve_f(integrator::SDEIntegrator{A}) where {A <: ImplicitTauLeaping}
    # Determine if the cache is in-place or constant (out-of-place) based on cache type
    cache = integrator.cache
    if cache isa ImplicitTauLeapingCache
        # In-place version - use rate cache from the integrator's cache
        rate_cache = cache.rate_at_uprev
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), typeof(rate_cache), true}(
            integrator.c, integrator.P.cache.rate, rate_cache
        )
    else
        # Out-of-place (constant cache) version
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), Nothing, false}(
            integrator.c, integrator.P.cache.rate, nothing
        )
    end
end

# nlsolve_f override for ThetaTrapezoidalTauLeaping: return TauLeapingDrift wrapper
function OrdinaryDiffEqCore.nlsolve_f(integrator::SDEIntegrator{A}) where {A <: ThetaTrapezoidalTauLeaping}
    # Determine if the cache is in-place or constant (out-of-place) based on cache type
    cache = integrator.cache
    if cache isa ThetaTrapezoidalTauLeapingCache
        # In-place version - use rate cache from the integrator's cache
        rate_cache = cache.rate_at_uprev
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), typeof(rate_cache), true}(
            integrator.c, integrator.P.cache.rate, rate_cache
        )
    else
        # Out-of-place (constant cache) version
        return TauLeapingDrift{typeof(integrator.c), typeof(integrator.P.cache.rate), Nothing, false}(
            integrator.c, integrator.P.cache.rate, nothing
        )
    end
end

function iip_generate_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits)
    if alg.nlsolve isa NLNewton
        nf = nlsolve_f(f, alg)
        islin = f isa Union{SDEFunction, SplitSDEFunction} && islinear(nf.f)
        if islin
            J = nf.f
            W = WOperator{true}(f.mass_matrix, dt, J, _vec(u))
        else
            if ArrayInterface.isstructured(f.jac_prototype) ||
                    f.jac_prototype isa SparseMatrixCSC
                J = similar(f.jac_prototype)
                W = similar(J)
            elseif DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) &&
                    f.jac_prototype !== nothing
                J = deepcopy(f.jac_prototype)
                if J isa AbstractMatrix
                    J = MatrixOperator(J; update_func! = f.jac)
                end
                W = WOperator{true}(f.mass_matrix, dt, J, _vec(u))
            else
                J = false .* vec(u) .* vec(u)'
                W = similar(J)
            end
        end
    else
        J = nothing
        W = nothing
    end
    return J, W
end

function oop_generate_W(alg, u, uprev, p, t, dt, f, uEltypeNoUnits)
    nf = nlsolve_f(f, alg)
    islin = f isa Union{SDEFunction, SplitSDEFunction} && islinear(nf.f)
    if islin || DiffEqBase.has_jac(f)
        # get the operator
        J = islin ? nf.f : f.jac(uprev, p, t)
        if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
            J = MatrixOperator(J)
        end
        W = WOperator{false}(f.mass_matrix, dt, J, _vec(u))
    else
        if u isa StaticArray
            # get a "fake" `J`
            J = if u isa AbstractMatrix && size(u, 1) > 1 # `u` is already a matrix
                u
            elseif size(u, 1) == 1 # `u` is a row vector
                vcat(u, u)
            else # `u` is a column vector
                hcat(u, u)
            end
            W = lu(J)
        else
            W = u isa Number ? u :
                LU{LinearAlgebra.lutype(uEltypeNoUnits)}(
                    Matrix{uEltypeNoUnits}(undef, 0, 0),
                    Vector{LinearAlgebra.BlasInt}(undef, 0),
                    zero(LinearAlgebra.BlasInt)
                )
            J = u isa Number ? u : (false .* vec(u) .* vec(u)')
        end
    end
    return J, W
end

# ============================================================================
# Traits for SDE composite algorithms/caches
# ============================================================================

# Trait: SDE composite algorithm types
function OrdinaryDiffEqCore.is_composite_algorithm(
        alg::Union{StochasticDiffEqCompositeAlgorithm, StochasticDiffEqRODECompositeAlgorithm},
    )
    return true
end

# Trait: SDE composite cache type
OrdinaryDiffEqCore.is_composite_cache(cache::StochasticCompositeCache) = true

# ============================================================================
# Noise interface methods for OrdinaryDiffEqCore's noise functions
# These implement the accept_noise!/reject_noise!/save_noise!/noise_curt/
# is_noise_saveable interface that ODE's unified loop functions call.
# ============================================================================

function OrdinaryDiffEqCore.accept_noise!(W::DiffEqNoiseProcess.AbstractNoiseProcess, dt, u, p, setup)
    return DiffEqNoiseProcess.accept_step!(W, dt, u, p, setup)
end

function OrdinaryDiffEqCore.reject_noise!(W::DiffEqNoiseProcess.AbstractNoiseProcess, dt, u, p)
    return DiffEqNoiseProcess.reject_step!(W, dt, u, p)
end

function OrdinaryDiffEqCore.save_noise!(W::DiffEqNoiseProcess.AbstractNoiseProcess)
    return DiffEqNoiseProcess.save_noise!(W)
end

OrdinaryDiffEqCore.noise_curt(W::DiffEqNoiseProcess.AbstractNoiseProcess) = W.curt

OrdinaryDiffEqCore.is_noise_saveable(W::NoiseProcess) = true
OrdinaryDiffEqCore.is_noise_saveable(W::DiffEqNoiseProcess.AbstractNoiseProcess) = false
