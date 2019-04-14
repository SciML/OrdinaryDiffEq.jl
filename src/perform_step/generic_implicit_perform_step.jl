mutable struct ImplicitRHS_Scalar{F,uType,tType,P} <: Function
  f::F
  tmp::uType
  a::tType
  t::tType
  dt::tType
  p::P
end

function (p::ImplicitRHS_Scalar)(resid,u)
  resid[1] = first(u) - p.tmp - p.a*first(p.f(first(u),p.p,p.t+p.dt))
end

mutable struct ImplicitRHS{F,uType,tType,DiffCacheType,P} <: Function
  f::F
  tmp::uType
  a::tType
  t::tType
  dt::tType
  dual_cache::DiffCacheType
  p::P
end

function (p::ImplicitRHS)(resid,u)
  _du1 = get_du(p.dual_cache, eltype(u))
  du1 = reinterpret(eltype(u),_du1)
  p.f(du1,u,p.p,p.t+p.dt)
  @.. resid = u - p.tmp - p.a*du1
end

function initialize!(integrator,
                     cache::Union{GenericImplicitEulerConstantCache,GenericTrapezoidConstantCache})
  cache.uhold[1] = integrator.uprev
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::GenericImplicitEulerConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uhold,rhs,nl_rhs = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  rhs.tmp = uprev

  if integrator.success_iter > 0 && !integrator.reeval_fsal && alg.extrapolant == :interpolant
    uhold[1] = current_extrapolant(t+dt,integrator)
  elseif alg.extrapolant == :linear
    uhold[1] = uprev + dt*integrator.fsalfirst
  else # :constant
    uhold[1] = uprev
  end

  rhs.t = t
  rhs.dt = dt
  rhs.a = dt
  nlres = alg.nlsolve(nl_rhs,uhold)
  uhold[1] = nlres[1]
  integrator.fsallast = f(uhold[1],p,t+dt)
  integrator.destats.nf += 1
  u = uhold[1]

  if integrator.opts.adaptive && integrator.success_iter > 0
    # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
    # use 2nd divided differences (DD) a la SPICE and Shampine

    # TODO: check numerical stability
    uprev2 = integrator.uprev2
    tprev = integrator.tprev

    dt1 = dt*(t+dt-tprev)
    dt2 = (t-tprev)*(t+dt-tprev)
    c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
    r = c*dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

    tmp = r*abs((u - uprev)/dt1 - (uprev - uprev2)/dt2)
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  else
    integrator.EEst = 1
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,
                     cache::Union{GenericImplicitEulerCache,GenericTrapezoidCache})
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1

  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::GenericImplicitEulerCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack dual_cache,k,nl_rhs,rhs,tmp,atmp = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  copyto!(tmp,uprev)

  if integrator.success_iter > 0 && !integrator.reeval_fsal && alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
  elseif alg.extrapolant == :linear
    @.. u = uprev + dt*integrator.fsalfirst
  else
    copyto!(u,uprev)
  end

  rhs.t = t
  rhs.dt = dt
  rhs.a = dt
  nlres = alg.nlsolve(nl_rhs,u)
  copyto!(u,nlres)

  if integrator.opts.adaptive && integrator.success_iter > 0
    # local truncation error (LTE) bound by dt^2/2*max|y''(t)|
    # use 2nd divided differences (DD) a la SPICE and Shampine

    # TODO: check numerical stability
    uprev2 = integrator.uprev2
    tprev = integrator.tprev

    dt1 = dt*(t+dt-tprev)
    dt2 = (t-tprev)*(t+dt-tprev)
    c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
    r = c*dt^2 # by mean value theorem 2nd DD equals y''(s)/2 for some s

    @.. tmp = r*abs((u - uprev)/dt1 - (uprev - uprev2)/dt2)
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  else
    integrator.EEst = 1
  end

  f(k, u, p, t+dt)
  integrator.destats.nf += 1
end

function initialize!(integrator, cache::GenericTrapezoidConstantCache)
  cache.uhold[1] = integrator.uprev
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::GenericTrapezoidConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack uhold,rhs,nl_rhs = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  rhs.tmp = first(uprev) + (dt/2)*first(integrator.fsalfirst)

  if integrator.success_iter > 0 && !integrator.reeval_fsal && alg.extrapolant == :interpolant
    uhold[1] = current_extrapolant(t+dt,integrator)
  elseif alg.extrapolant == :linear
    uhold[1] = uprev + dt*integrator.fsalfirst
  else # :constant
    uhold[1] = uprev
  end

  rhs.t = t
  rhs.dt = dt
  rhs.a = dt/2
  nlres = alg.nlsolve(nl_rhs,uhold)
  uhold[1] = nlres[1]
  integrator.fsallast = f(uhold[1],p,t+dt)
  integrator.destats.nf += 1
  u = uhold[1]

  if integrator.opts.adaptive
    if integrator.iter > 2
      # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
      # use 3rd divided differences (DD) a la SPICE and Shampine

      # TODO: check numerical stability
      uprev2 = integrator.uprev2
      tprev = integrator.tprev
      uprev3 = cache.uprev3
      tprev2 = cache.tprev2

      dt1 = dt*(t+dt-tprev)
      dt2 = (t-tprev)*(t+dt-tprev)
      dt3 = (t-tprev)*(t-tprev2)
      dt4 = (tprev-tprev2)*(t-tprev2)
      dt5 = t+dt-tprev2
      c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
      r = c*dt^3/2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

      # tmp = r*abs(((u - uprev)/dt1 - (uprev - uprev2)/dt2) - ((uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4)/dt5)
      DD31 = (u - uprev)/dt1 - (uprev - uprev2)/dt2
      DD30 = (uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4
      tmp = r*abs((DD31 - DD30)/dt5)
      atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
      integrator.EEst = integrator.opts.internalnorm(atmp,t)
      if integrator.EEst <= 1
        cache.uprev3 = uprev2
        cache.tprev2 = tprev
      end
    elseif integrator.success_iter > 0
      integrator.EEst = 1
      cache.uprev3 = integrator.uprev2
      cache.tprev2 = integrator.tprev
    else
      integrator.EEst = 1
    end
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::GenericTrapezoidCache)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t)
  integrator.destats.nf += 1
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::GenericTrapezoidCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack dual_cache,k,rhs,nl_rhs,tmp,atmp = cache
  alg = typeof(integrator.alg) <: CompositeAlgorithm ? integrator.alg.algs[integrator.cache.current] : integrator.alg
  tmp .= uprev .+ (dt/2).*integrator.fsalfirst

  if integrator.success_iter > 0 && !integrator.reeval_fsal && alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
  elseif alg.extrapolant == :linear
    @.. u = uprev + dt*integrator.fsalfirst
  else
    copyto!(u,uprev)
  end

  # copyto!(rhs.fsalfirst,fsalfirst) Implicitly done by pointers: fsalfirst === fsalfirst == rhs.fsalfirst
  rhs.t = t
  rhs.dt = dt
  rhs.a = dt/2
  nlres = alg.nlsolve(nl_rhs,u)
  copyto!(u,nlres)

  if integrator.opts.adaptive
    if integrator.iter > 2
      # local truncation error (LTE) bound by dt^3/12*max|y'''(t)|
      # use 3rd divided differences (DD) a la SPICE and Shampine

      # TODO: check numerical stability
      uprev2 = integrator.uprev2
      tprev = integrator.tprev
      uprev3 = cache.uprev3
      tprev2 = cache.tprev2

      dt1 = dt*(t+dt-tprev)
      dt2 = (t-tprev)*(t+dt-tprev)
      dt3 = (t-tprev)*(t-tprev2)
      dt4 = (tprev-tprev2)*(t-tprev2)
      dt5 = t+dt-tprev2
      c = 7/12 # default correction factor in SPICE (LTE overestimated by DD)
      r = c*dt^3/2 # by mean value theorem 3rd DD equals y'''(s)/6 for some s

      # @.. tmp = r*abs(((u - uprev)/dt1 - (uprev - uprev2)/dt2) - ((uprev - uprev2)/dt3 - (uprev2 - uprev3)/dt4)/dt5)
      @inbounds for i in eachindex(u)
        DD31 = (u[i] - uprev[i])/dt1 - (uprev[i] - uprev2[i])/dt2
        DD30 = (uprev[i] - uprev2[i])/dt3 - (uprev2[i] - uprev3[i])/dt4
        tmp[i] = r*abs((DD31 - DD30)/dt5)
      end
      calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
      integrator.EEst = integrator.opts.internalnorm(atmp,t)
      if integrator.EEst <= 1
        copyto!(cache.uprev3,uprev2)
        cache.tprev2 = tprev
      end
    elseif integrator.success_iter > 0
      integrator.EEst = 1
      copyto!(cache.uprev3,integrator.uprev2)
      cache.tprev2 = integrator.tprev
    else
      integrator.EEst = 1
    end
  end

  f(k, u, p, t+dt)
  integrator.destats.nf += 1
end
