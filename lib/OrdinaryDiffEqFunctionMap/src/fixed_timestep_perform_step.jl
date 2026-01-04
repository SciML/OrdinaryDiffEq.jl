function initialize!(integrator, cache::FunctionMapCache)
    integrator.kshortsize = 0
    return resize!(integrator.k, integrator.kshortsize)
end

function perform_step!(integrator, cache::FunctionMapCache, repeat_step = false)
    (; u, uprev, dt, t, f, p) = integrator
    alg = unwrap_alg(integrator, nothing)
    (; tmp) = cache
    return if integrator.f != DiffEqBase.DISCRETE_INPLACE_DEFAULT &&
            !(
            integrator.f isa DiffEqBase.EvalFunc &&
                integrator.f.f === DiffEqBase.DISCRETE_INPLACE_DEFAULT
        )
        if FunctionMap_scale_by_time(alg)
            f(tmp, uprev, p, t + dt)
            @muladd @.. broadcast = false u = uprev + dt * tmp
        else
            f(u, uprev, p, t + dt)
        end
        OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
    end
end

function perform_step!(integrator, cache::FunctionMapConstantCache, repeat_step = false)
    (; uprev, dt, t, f, p) = integrator
    alg = unwrap_alg(integrator, nothing)
    return if integrator.f != DiffEqBase.DISCRETE_OUTOFPLACE_DEFAULT &&
            !(
            integrator.f isa DiffEqBase.EvalFunc &&
                integrator.f.f === DiffEqBase.DISCRETE_OUTOFPLACE_DEFAULT
        )
        if FunctionMap_scale_by_time(alg)
            tmp = f(uprev, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
            @muladd integrator.u = @.. broadcast = false uprev + dt * tmp
        else
            integrator.u = f(uprev, p, t + dt)
            OrdinaryDiffEqCore.increment_nf!(integrator.stats, 1)
        end
    end
end

function initialize!(integrator, cache::FunctionMapConstantCache)
    integrator.kshortsize = 0
    return integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
end
