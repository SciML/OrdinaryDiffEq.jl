function initialize!(integrator,cache::AN5ConstantCache)
  integrator.kshortsize = 16
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::AN5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,l,m,tq,tau,vern9tab = cache
  # handle callbacks, rewind back to order one.
  if integrator.u_modified
    cache.step = 1
  end
  # Nordsieck form needs to build the history vector
  if cache.step == 1
    # Start the Nordsieck vector in one shot!
    perform_step!(integrator, vern9tab, repeat_step)
    cache.step = 5
    z[1] = integrator.uprev
    z[2] = integrator.k[1]*dt
    ode_addsteps!(integrator.k,t,integrator.uprev,integrator.u,dt,f,p,vern9tab)
    z[3] = ode_interpolant(t,dt,nothing,nothing,integrator.k,vern9tab,nothing,Val{2})*dt^2/2
    z[4] = ode_interpolant(t,dt,nothing,nothing,integrator.k,vern9tab,nothing,Val{3})*dt^3/6
    z[5] = ode_interpolant(t,dt,nothing,nothing,integrator.k,vern9tab,nothing,Val{4})*dt^4/24
    z[6] = zero(cache.z[6])
    fill!(tau, dt)
  end
  # Reset time
  for i in endof(tau):-1:2
    tau[i] = tau[i-1]
  end
  tau[1] = dt
  dt != tau[2] && nordsieck_rescale!(cache, false)
  integrator.k[1] = z[2]/dt
  # Perform 5th order Adams method in Nordsieck form
  perform_predict!(cache, false)
  calc_coeff!(cache)
  isucceed = nlsolve_functional!(integrator, cache)
  isucceed || return nothing

  ################################### Error estimation

  if integrator.opts.adaptive
    tmp = cache.Δ*cache.tq
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
    if integrator.EEst >= one(integrator.EEst)
      perform_predict!(cache, true)
      nordsieck_rescale!(cache, true)
      return nothing
    end
  end

  # Corrector
  perform_correct!(cache)

  ################################### Finalize

  integrator.k[2] = cache.z[2]/dt
  return nothing
end

function initialize!(integrator, cache::AN5Cache)
  integrator.kshortsize = 16
  integrator.fsalfirst = cache.vern9cache.k1; integrator.fsallast = cache.vern9cache.k7 # setup pointers
  resize!(integrator.k, integrator.kshortsize)
  # Setup k pointers
  integrator.k[1] = cache.vern9cache.k1
  integrator.k[2] = cache.vern9cache.k2
  integrator.k[3] = cache.vern9cache.k3
  integrator.k[4] = cache.vern9cache.k4
  integrator.k[5] = cache.vern9cache.k5
  integrator.k[6] = cache.vern9cache.k6
  integrator.k[7] = cache.vern9cache.k7
  integrator.k[8] = cache.vern9cache.k8
  integrator.k[9] = cache.vern9cache.k9
  integrator.k[10] = cache.vern9cache.k10
  integrator.k[11] = cache.vern9cache.k11
  integrator.k[12] = cache.vern9cache.k12
  integrator.k[13] = cache.vern9cache.k13
  integrator.k[14] = cache.vern9cache.k14
  integrator.k[15] = cache.vern9cache.k15
  integrator.k[16] = cache.vern9cache.k16
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
end

@muladd function perform_step!(integrator, cache::AN5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,uprev2 = integrator
  @unpack const_cache,utilde,tmp,ratetmp,atmp,vern9cache = cache
  @unpack z,l,m,tq,tau, = const_cache
  # handle callbacks, rewind back to order one.
  if integrator.u_modified
    const_cache.step = 1
  end
  # Nordsieck form needs to build the history vector
  if const_cache.step == 1
    ## Start the Nordsieck vector in two shots!
    perform_step!(integrator, vern9cache, repeat_step)
    const_cache.step = 5
    @. z[1] = integrator.uprev
    @. z[2] = integrator.k[1]*dt
    ode_addsteps!(integrator.k, t, uprev, u, dt, f, p, vern9cache)
    ode_interpolant!(z[3],t,dt,nothing,nothing,integrator.k,vern9cache,nothing,Val{2})
    ode_interpolant!(z[4],t,dt,nothing,nothing,integrator.k,vern9cache,nothing,Val{3})
    ode_interpolant!(z[5],t,dt,nothing,nothing,integrator.k,vern9cache,nothing,Val{4})
    @. z[3] = z[3]*dt^2/2
    @. z[4] = z[4]*dt^3/6
    @. z[5] = z[5]*dt^4/24
    fill!(z[6], zero(eltype(z[6])))
    fill!(tau, dt)
  end

  # Reset time
  for i in endof(tau):-1:2
    tau[i] = tau[i-1]
  end
  tau[1] = dt
  # Rescale
  dt != tau[2] && nordsieck_rescale!(cache)
  # Perform 5th order Adams method in Nordsieck form
  perform_predict!(cache, false)
  calc_coeff!(cache)
  isucceed = nlsolve_functional!(integrator, cache)
  isucceed || return nothing

  ################################### Error estimation

  if integrator.opts.adaptive
    @. tmp = const_cache.Δ*const_cache.tq
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
    if integrator.EEst >= one(integrator.EEst)
      perform_predict!(cache, true)
      nordsieck_rescale!(cache, true)
      return nothing
    end
  end

  # Corrector
  perform_correct!(cache)

  ################################### Finalize

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = const_cache.z[2]/dt
  return nothing
end
