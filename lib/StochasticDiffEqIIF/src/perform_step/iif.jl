mutable struct RHS_IIF1M_Scalar{F, CType, tType, P} <: Function
    f::F
    t::tType
    dt::tType
    tmp::CType
    p::P
end

function (f::RHS_IIF1M_Scalar)(resid, u)
    return resid[1] = u[1] - f.tmp - f.dt * f.f.f2(u[1], f.p, f.t + f.dt)[1]
end

mutable struct RHS_IIF2M_Scalar{F, CType, tType, P} <: Function
    f::F
    t::tType
    dt::tType
    tmp::CType
    p::P
end

function (f::RHS_IIF2M_Scalar)(resid, u)
    return resid[1] = u[1] - f.tmp - 0.5f.dt * f.f.f2(u[1], f.p, f.t + f.dt)[1]
end

@muladd function initialize!(
        integrator, cache::Union{
            IIF1MConstantCache, IIF2MConstantCache, IIF1MilConstantCache,
        }
    )
    cache.uhold[1] = integrator.uprev
end

@muladd function perform_step!(
        integrator, cache::Union{
            IIF1MConstantCache, IIF2MConstantCache, IIF1MilConstantCache,
        }
    )
    (; t, dt, uprev, u, W, p, f) = integrator
    (; uhold, rhs, nl_rhs) = cache
    alg = unwrap_alg(integrator, true)
    A = integrator.f.f1(u, p, t)
    if cache isa IIF1MilConstantCache
        error("Milstein correction does not work.")
    elseif cache isa IIF1MConstantCache
        tmp = exp(A * dt) * (uprev + integrator.f.g(uprev, p, t) * W.dW)
    elseif cache isa IIF2MConstantCache
        tmp = exp(A * dt) * (uprev + 0.5dt * integrator.f.f2(uprev, p, t) + integrator.f.g(uprev, p, t) * W.dW)
    end

    if integrator.iter > 1 && !integrator.u_modified
        uhold[1] = current_extrapolant(t + dt, integrator)
    end # else uhold is previous value.

    rhs.t = t
    rhs.dt = dt
    rhs.tmp = tmp
    nlres = alg.nlsolve(nl_rhs, uhold)

    u = nlres[1]
    integrator.u = u
end

mutable struct RHS_IIF1{F, uType, tType, DiffCacheType, SizeType, P} <: Function
    f::F
    tmp::uType
    t::tType
    dt::tType
    dual_cache::DiffCacheType
    sizeu::SizeType
    p::P
end
function (f::RHS_IIF1)(resid, u)
    _du = get_du(f.dual_cache, eltype(u))
    du = reinterpret(eltype(u), _du)
    f.f.f2(du, reshape(u, f.sizeu), f.p, f.t + f.dt)
    return @.. resid = u - f.tmp - f.dt * du
end

mutable struct RHS_IIF2{F, uType, tType, DiffCacheType, SizeType, P} <: Function
    f::F
    tmp::uType
    t::tType
    dt::tType
    dual_cache::DiffCacheType
    sizeu::SizeType
    p::P
end
function (f::RHS_IIF2)(resid, u)
    _du = get_du(f.dual_cache, eltype(u))
    du = reinterpret(eltype(u), _du)
    f.f.f2(du, reshape(u, f.sizeu), f.p, f.t + f.dt)
    return @.. resid = u - f.tmp - 0.5f.dt * du
end

@muladd function perform_step!(integrator, cache::Union{IIF1MCache, IIF2MCache})
    (; rtmp1, rtmp2, rtmp3, tmp, noise_tmp) = cache
    (; uhold, rhs, nl_rhs) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    alg = unwrap_alg(integrator, true)
    uidx = eachindex(u)

    integrator.f.g(rtmp2, uprev, p, t)

    if is_diagonal_noise(integrator.sol.prob)
        rmul!(rtmp2, W.dW) # rtmp2 === rtmp3
    else
        mul!(rtmp3, rtmp2, W.dW)
    end

    rtmp3 .+= uprev

    if cache isa IIF2MCache
        integrator.f.f2(rtmp1, uprev, p, t)
        dto2 = dt / 2
        @.. rtmp3 = dto2 * rtmp1 + rtmp3
    end

    A = integrator.f.f1(rtmp1, uprev, p, t)
    M = exp(A * dt)
    mul!(tmp, M, rtmp3)

    if integrator.iter > 1 && !integrator.u_modified
        current_extrapolant!(uhold, t + dt, integrator)
    end # else uhold is previous value.

    rhs.t = t
    rhs.dt = dt
    rhs.tmp = tmp
    rhs.sizeu = size(u)
    nlres = alg.nlsolve(nl_rhs, uhold)

    copyto!(uhold, nlres)
end

@muladd function perform_step!(integrator, cache::IIF1MilCache)
    (; rtmp1, rtmp2, rtmp3, tmp, noise_tmp) = cache
    (; uhold, rhs, nl_rhs) = cache
    (; t, dt, uprev, u, W, p, f) = integrator
    alg = unwrap_alg(integrator, true)

    dW = W.dW
    sqdt = integrator.sqdt
    f = integrator.f
    g = integrator.f.g

    A = integrator.f.f1(t, uprev, rtmp1)
    M = exp(A * dt)

    uidx = eachindex(u)
    integrator.f.g(rtmp2, uprev, p, t)
    if cache isa Union{IIF1MCache, IIF2MCache}
        if is_diagonal_noise(integrator.sol.prob)
            rmul!(rtmp2, W.dW) # rtmp2 === rtmp3
        else
            mul!(rtmp3, rtmp2, W.dW)
        end
    else #Milstein correction
        rtmp2 = M * rtmp2 # mul!(rtmp2,M,gtmp)
        (; gtmp, gtmp2) = cache
        #error("Milstein correction does not work.")
        mul!(rtmp3, rtmp2, W.dW)
        I = zeros(length(dW), length(dW))
        Dg = zeros(length(dW), length(dW))
        mil_correction = zeros(length(dW))
        mil_correction .= 0.0
        for i in 1:length(dW), j in 1:length(dW)

            I[j, i] = 0.5 * dW[i] * dW[j]
            j == i && (I[i, i] -= 0.5 * dt) # Ito correction
        end
        for j in 1:length(uprev)
            #Kj = uprev .+ dt.*du1 + sqdt*rtmp2[:,j] # This works too
            Kj = uprev .+ sqdt * rtmp2[:, j]
            g(gtmp, Kj, p, t)
            mul!(gtmp2, M, gtmp)
            Dgj = (gtmp2 - rtmp2) / sqdt
            mil_correction .+= Dgj * I[:, j]
        end
        rtmp3 .+= mil_correction
    end

    if cache isa IIF2MCache
        integrator.f.f2(t, uprev, rtmp1)
        dto2 = dt / 2
        @.. rtmp1 = dto2 * rtmp1 + uprev + rtmp3
        mul!(tmp, M, rtmp1)
    elseif !(cache isa IIF1MilCache)
        @.. rtmp1 = uprev + rtmp3
        mul!(tmp, M, rtmp1)
    else
        mul!(tmp, M, uprev)
        tmp .+= rtmp3
    end

    if integrator.iter > 1 && !integrator.u_modified
        current_extrapolant!(uhold, t + dt, integrator)
    end # else uhold is previous value.

    rhs.t = t
    rhs.dt = dt
    rhs.tmp = tmp
    rhs.sizeu = size(u)
    nlres = alg.nlsolve(nl_rhs, uhold)

    copyto!(uhold, nlres)
end
