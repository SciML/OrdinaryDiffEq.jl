# construct solver for fixed-point iterations
function build_fpsolver(
        alg, fpalg::Union{NLFunctional, NLAnderson}, u, uEltypeNoUnits,
        uBottomEltypeNoUnits, ::Val{true}
    )
    # no fixed-point iterations if the algorithm is constrained
    isconstrained(alg) && return

    # define unitless type
    uTolType = real(uBottomEltypeNoUnits)

    # could use integrator.cache.atmp if guaranteed to exist
    atmp = similar(u, uEltypeNoUnits)
    dz = similar(u)

    # build cache
    if fpalg isa NLFunctional
        fpcache = FPFunctionalCache(atmp, dz)
    elseif fpalg isa NLAnderson
        max_history = min(fpalg.max_history, fpalg.max_iter, length(u))
        Δz₊s = [zero(u) for i in 1:max_history]
        Q = Matrix{uEltypeNoUnits}(undef, length(u), max_history)
        R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
        γs = Vector{uEltypeNoUnits}(undef, max_history)

        dzold = zero(u)
        z₊old = zero(u)

        fpcache = FPAndersonCache(atmp, dz, dzold, z₊old, Δz₊s, Q, R, γs, 0, fpalg.aa_start, fpalg.droptol)
    end

    # build solver
    ηold = one(uTolType)

    return FPSolver{typeof(fpalg), true, uTolType, typeof(fpcache)}(fpalg, uTolType(fpalg.κ), uTolType(fpalg.fast_convergence_cutoff), ηold, 10000, fpalg.max_iter, SlowConvergence, fpcache, 0)
end

function build_fpsolver(
        alg, fpalg::Union{NLFunctional, NLAnderson}, u, uEltypeNoUnits,
        uBottomEltypeNoUnits, ::Val{false}
    )
    # no fixed-point iterations if the algorithm is constrained
    isconstrained(alg) && return

    # define unitless type
    uTolType = real(uBottomEltypeNoUnits)

    # build cache
    if fpalg isa NLFunctional
        fpcache = FPFunctionalConstantCache()
    elseif fpalg isa NLAnderson
        max_history = min(fpalg.max_history, fpalg.max_iter, length(u))
        Δz₊s = Vector{typeof(u)}(undef, max_history)
        Q = Matrix{uEltypeNoUnits}(undef, length(u), max_history)
        R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
        γs = Vector{uEltypeNoUnits}(undef, max_history)

        dz = u
        dzold = u
        z₊old = u

        fpcache = FPAndersonConstantCache(
            dz, dzold, z₊old, Δz₊s, Q, R, γs, 0,
            fpalg.aa_start,
            fpalg.droptol
        )
    end

    # build solver
    ηold = one(uTolType)

    return FPSolver{typeof(fpalg), false, uTolType, typeof(fpcache)}(
        fpalg, uTolType(fpalg.κ),
        uTolType(fpalg.fast_convergence_cutoff),
        ηold, 10_000,
        fpalg.max_iter,
        SlowConvergence, fpcache, 0
    )
end

## resize

function resize_fpsolver!(integrator::DDEIntegrator, i::Int)
    (; fpsolver) = integrator

    if fpsolver !== nothing
        resize!(fpsolver, integrator, i)
    end

    return nothing
end

function Base.resize!(fpsolver::FPSolver, integrator::DDEIntegrator, i::Int)
    return resize!(fpsolver.cache, fpsolver, integrator, i)
end
