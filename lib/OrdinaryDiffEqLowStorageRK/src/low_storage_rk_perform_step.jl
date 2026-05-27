function initialize!(integrator, cache::LowStorageRKTableau{TwoN})
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

function perform_step!(integrator, cache::LowStorageRKTableau{TwoN}, repeat_step = false)
    return _perform_step_oop!(integrator, cache)
end

get_fsalfirstlast(cache::LowStorageRK2NCache, u) = (nothing, nothing)

function initialize!(integrator, cache::LowStorageRK2NCache)
    (; k, tmp, williamson_condition) = cache
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = k
    integrator.f(k, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::LowStorageRK2NCache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end

function initialize!(integrator, cache::LowStorageRKTableau{TwoC})
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

function perform_step!(integrator, cache::LowStorageRKTableau{TwoC}, repeat_step = false)
    return _perform_step_oop!(integrator, cache)
end

function initialize!(integrator, cache::LowStorageRK2CCache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::LowStorageRK2CCache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end

# 3S low storage methods
function initialize!(integrator, cache::LowStorageRK3SConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

function perform_step!(integrator, cache::LowStorageRK3SConstantCache, repeat_step = false)
    return _perform_step_oop!(integrator, cache)
end

function initialize!(integrator, cache::LowStorageRK3SCache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::LowStorageRK3SCache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end

function initialize!(integrator, cache::LowStorageRK3SpConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

function perform_step!(integrator, cache::LowStorageRK3SpConstantCache, repeat_step = false)
    return _perform_step_oop!(integrator, cache)
end

function initialize!(integrator, cache::LowStorageRK3SpCache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    return integrator.k[1] = integrator.fsalfirst
end

function perform_step!(integrator, cache::LowStorageRK3SpCache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end

function initialize!(integrator, cache::LowStorageRK3SpFSALConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function perform_step!(
        integrator, cache::LowStorageRK3SpFSALConstantCache, repeat_step = false
    )
    return _perform_step_oop!(integrator, cache)
end

function initialize!(integrator, cache::LowStorageRK3SpFSALCache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::LowStorageRK3SpFSALCache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end

function initialize!(integrator, cache::LowStorageRK2RPConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

function perform_step!(integrator, cache::LowStorageRK2RPConstantCache, repeat_step = false)
    return _perform_step_oop!(integrator, cache)
end

function initialize!(integrator, cache::LowStorageRK2RPCache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::LowStorageRK2RPCache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end

function initialize!(integrator, cache::LowStorageRK3RPConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

function perform_step!(integrator, cache::LowStorageRK3RPConstantCache, repeat_step = false)
    return _perform_step_oop!(integrator, cache)
end

function initialize!(integrator, cache::LowStorageRK3RPCache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::LowStorageRK3RPCache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end


function initialize!(integrator, cache::LowStorageRK4RPConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

function perform_step!(integrator, cache::LowStorageRK4RPConstantCache, repeat_step = false)
    return _perform_step_oop!(integrator, cache)
end

function initialize!(integrator, cache::LowStorageRK4RPCache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::LowStorageRK4RPCache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end


function initialize!(integrator, cache::LowStorageRK5RPConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

function perform_step!(integrator, cache::LowStorageRK5RPConstantCache, repeat_step = false)
    return _perform_step_oop!(integrator, cache)
end

function initialize!(integrator, cache::LowStorageRK5RPCache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::LowStorageRK5RPCache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end

function initialize!(integrator, cache::RK46NLConstantCache)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.kshortsize = 1
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsallast = zero(integrator.fsalfirst)
    return integrator.k[1] = integrator.fsalfirst
end

function perform_step!(integrator, cache::RK46NLConstantCache, repeat_step = false)
    return _perform_step_oop!(integrator, cache)
end

function initialize!(integrator, cache::RK46NLCache)
    integrator.kshortsize = 1
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::RK46NLCache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end

function initialize!(integrator, cache::SHLDDRK52ConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::SHLDDRK52ConstantCache, repeat_step = false)
    return _perform_step_oop!(integrator, cache)
end

function initialize!(integrator, cache::SHLDDRK52Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::SHLDDRK52Cache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end

function initialize!(integrator, cache::SHLDDRK_2NConstantCache)
    integrator.kshortsize = 2
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
    OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    return integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator, cache::SHLDDRK_2NConstantCache, repeat_step = false)
    return _perform_step_oop!(integrator, cache)
end

function initialize!(integrator, cache::SHLDDRK_2NCache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    return OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
end

function perform_step!(integrator, cache::SHLDDRK_2NCache, repeat_step = false)
    return _perform_step_iip!(integrator, cache, cache.tab)
end
