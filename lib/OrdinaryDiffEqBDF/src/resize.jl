function resize_bdf_buffers!(buffers, i)
    for buffer in buffers
        buffer isa AbstractArray && ArrayInterface.ismutable(buffer) && resize!(buffer, i)
    end
    return nothing
end

function resize_non_user_cache!(
        integrator::OrdinaryDiffEqCore.ODEIntegrator,
        cache::FBDFCache, i::Int
    )
    resize!(integrator.u, i)
    resize!(integrator.uprev, i)
    # The generic resize runs before integrator.u has its new size.
    resize_nlsolver!(integrator, i)
    resize_bdf_buffers!(cache.u_history, i)
    resize_bdf_buffers!(cache.u_corrector, i)
    resize_bdf_buffers!(cache.dense, i)
    cache.order = 1
    cache.prev_order = 1
    cache.nconsteps = 0
    cache.consfailcnt = 0
    cache.qwait = 3
    cache.iters_from_event = 0
    stald_reset!(cache.stald)
    # A resize invalidates multistep history, including when called outside a callback.
    integrator.derivative_discontinuity = true
    return nothing
end

function resize_non_user_cache!(
        integrator::OrdinaryDiffEqCore.ODEIntegrator, cache::QNDFCache, i::Int
    )
    resize_bdf_buffers!(cache.D, i)
    resize_bdf_buffers!(cache.Dtmp, i)
    resize_bdf_buffers!(cache.prevD, i)
    resize_bdf_buffers!(cache.dense, i)
    cache.order = 1
    cache.prevorder = 1
    cache.nconsteps = 0
    cache.consfailcnt = 0
    integrator.derivative_discontinuity = true
    return nothing
end
