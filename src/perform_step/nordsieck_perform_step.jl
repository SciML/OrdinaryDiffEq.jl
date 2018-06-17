function initialize!(integrator,cache::AN5ConstantCache)
  integrator.kshortsize = 7
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
  @unpack z,l,m,c_LTE,tau,tsit5tab = cache
  # handle callbacks, rewind back to order one.
  if integrator.u_modified
    cache.step = 1
  end
  # Nordsieck form needs to build the history vector
  if cache.step == 1
    # Start the Nordsieck vector in one shot!
    perform_step!(integrator, tsit5tab, repeat_step)
    cache.step = 4
    z[1] = integrator.uprev
    z[2] = integrator.k[1]*dt
    z[3] = ode_interpolant(t,dt,nothing,nothing,integrator.k,tsit5tab,nothing,Val{2})*dt^2/2
    z[4] = ode_interpolant(t,dt,nothing,nothing,integrator.k,tsit5tab,nothing,Val{3})*dt^3/6
    z[5] = ode_interpolant(t,dt,nothing,nothing,integrator.k,tsit5tab,nothing,Val{4})*dt^4/24
    z[6] = zero(cache.z[6])
    fill!(tau, dt)
    perform_predict!(cache)
    cache.Δ = integrator.u - integrator.uprev
    update_nordsieck_vector!(cache)
    if integrator.opts.adaptive && integrator.EEst >= one(integrator.EEst)
      cache.step = 1
    end
  else
    # Reset time
    for i in endof(tau):-1:2
      tau[i] = tau[i-1]
    end
    tau[1] = dt
    dt != tau[2] && nordsieck_rescale!(cache)
    integrator.k[1] = z[2]/dt
    # Perform 5th order Adams method in Nordsieck form
    perform_predict!(cache)
    calc_coeff!(cache)
    isucceed = nlsolve_functional!(integrator, cache)
    if !isucceed
      # rewind Nordsieck vector
      integrator.force_stepfail = true
      nordsieck_rewind!(cache)
      return nothing
    end

    ################################### Error estimation

    if integrator.opts.adaptive
      atmp = calculate_residuals(cache.Δ, uprev, integrator.u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp) * cache.c_LTE
    end

    # Correct Nordsieck vector
    cache.step = min(cache.step+1, 5)
    update_nordsieck_vector!(cache)

    ################################### Finalize

    integrator.k[2] = cache.z[2]/dt
  end
  return nothing
end

function initialize!(integrator, cache::AN5Cache)
  integrator.kshortsize = 7
  integrator.fsalfirst = cache.tsit5cache.k1; integrator.fsallast = cache.tsit5cache.k7 # setup pointers
  resize!(integrator.k, integrator.kshortsize)
  # Setup k pointers
  integrator.k[1] = cache.tsit5cache.k1
  integrator.k[2] = cache.tsit5cache.k2
  integrator.k[3] = cache.tsit5cache.k3
  integrator.k[4] = cache.tsit5cache.k4
  integrator.k[5] = cache.tsit5cache.k5
  integrator.k[6] = cache.tsit5cache.k6
  integrator.k[7] = cache.tsit5cache.k7
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
end

@muladd function perform_step!(integrator, cache::AN5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,uprev2 = integrator
  @unpack const_cache,utilde,tmp,ratetmp,atmp,tsit5cache = cache
  @unpack z,l,m,c_LTE,tau, = const_cache
  # handle callbacks, rewind back to order one.
  if integrator.u_modified
    const_cache.step = 1
  end
  # Nordsieck form needs to build the history vector
  if const_cache.step == 1
    ## Start the Nordsieck vector in two shots!
    perform_step!(integrator, tsit5cache, repeat_step)
    copy!(tmp, integrator.u)
    const_cache.step = 4
    @. z[1] = integrator.uprev
    @. z[2] = integrator.k[1]*dt
    ode_interpolant!(z[3],t,dt,nothing,nothing,integrator.k,tsit5cache,nothing,Val{2})
    ode_interpolant!(z[4],t,dt,nothing,nothing,integrator.k,tsit5cache,nothing,Val{3})
    ode_interpolant!(z[5],t,dt,nothing,nothing,integrator.k,tsit5cache,nothing,Val{4})
    @. z[3] = z[3]*dt^2/2
    @. z[4] = z[4]*dt^3/6
    @. z[5] = z[5]*dt^4/24
    fill!(z[6], 0)
    fill!(tau, dt)
    perform_predict!(cache)
    @. const_cache.Δ = integrator.u - integrator.uprev
    update_nordsieck_vector!(cache)
    if integrator.opts.adaptive && integrator.EEst >= one(integrator.EEst)
      const_cache.step = 1
    end
  else
    # Reset time
    for i in endof(tau):-1:2
      tau[i] = tau[i-1]
    end
    tau[1] = dt
    # Rescale
    dt != tau[2] && nordsieck_rescale!(cache)
    @. integrator.k[1] = z[2]/dt
    # Perform 5th order Adams method in Nordsieck form
    perform_predict!(cache)
    calc_coeff!(cache)
    isucceed = nlsolve_functional!(integrator, cache)
    if !isucceed
      integrator.force_stepfail = true
      # rewind Nordsieck vector
      nordsieck_rewind!(cache)
      return nothing
    end

    ################################### Error estimation

    if integrator.opts.adaptive
      calculate_residuals!(atmp, const_cache.Δ, uprev, integrator.u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp) * const_cache.c_LTE
    end

    # Correct Nordsieck vector
    const_cache.step = min(const_cache.step+1, 5)
    update_nordsieck_vector!(cache)

    ################################### Finalize

    @. integrator.k[2] = const_cache.z[2]/dt
  end
  return nothing
end

function initialize!(integrator,cache::JVODEConstantCache)
  integrator.kshortsize = 7
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

@muladd function perform_step!(integrator, cache::JVODEConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,l,m,c_LTE,tau,tsit5tab = cache
  # handle callbacks, rewind back to order one.
  if integrator.u_modified
    cache.step = 1
  end
  if cache.step == 1
    z[1] = integrator.uprev
    z[2] = f(uprev, p, t)*dt
    z[3] = zero(cache.z[3])
    fill!(tau, dt)
    perform_predict!(cache)
    calc_coeff!(cache)
    isucceed = nlsolve_functional!(integrator, cache)
    if !isucceed
      # rewind Nordsieck vector
      integrator.force_stepfail = true
      nordsieck_rewind!(cache)
      return nothing
    end
    update_nordsieck_vector!(cache)
  else
    # Nordsieck form needs to build the history vector
    # Reset time
    for i in endof(tau):-1:2
      tau[i] = tau[i-1]
    end
    tau[1] = dt
    dt != tau[2] && nordsieck_rescale!(cache)
    integrator.k[1] = z[2]/dt
    # Perform 5th order Adams method in Nordsieck form
    perform_predict!(cache)
    calc_coeff!(cache)
    isucceed = nlsolve_functional!(integrator, cache)
    if !isucceed
      # rewind Nordsieck vector
      integrator.force_stepfail = true
      nordsieck_rewind!(cache)
      return nothing
    end

    # Correct Nordsieck vector
    cache.step = min(cache.nextorder, 12)
    update_nordsieck_vector!(cache)

    ################################### Finalize

    integrator.k[2] = cache.z[2]/dt
  end
  ################################### Error estimation
  if integrator.opts.adaptive
    atmp = calculate_residuals(cache.Δ, uprev, integrator.u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp) * cache.c_LTE
  end
  cache.n_wait -= 1
  cache.prev_𝒟 = cache.c_𝒟
  return nothing
end

function initialize!(integrator, cache::JVODECache)
  integrator.kshortsize = 7
  integrator.fsalfirst = cache.tsit5cache.k1; integrator.fsallast = cache.tsit5cache.k7 # setup pointers
  resize!(integrator.k, integrator.kshortsize)
  # Setup k pointers
  integrator.k[1] = cache.tsit5cache.k1
  integrator.k[2] = cache.tsit5cache.k2
  integrator.k[3] = cache.tsit5cache.k3
  integrator.k[4] = cache.tsit5cache.k4
  integrator.k[5] = cache.tsit5cache.k5
  integrator.k[6] = cache.tsit5cache.k6
  integrator.k[7] = cache.tsit5cache.k7
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
end

@muladd function perform_step!(integrator, cache::JVODECache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,uprev2 = integrator
  @unpack const_cache,utilde,tmp,ratetmp,atmp,tsit5cache = cache
  @unpack z,l,m,c_LTE,tau, = const_cache
  # handle callbacks, rewind back to order one.
  if integrator.u_modified
    const_cache.step = 1
  end
  if const_cache.step == 1
    perform_step!(integrator, EulerCache(u, uprev, tmp, z[2], integrator.fsalfirst), repeat_step)
    @. z[1] = integrator.uprev
    @. z[2] = integrator.k[1]*dt
    @. z[3] = 0
    fill!(tau, dt)
    perform_predict!(cache)
    const_cache.Δ = integrator.u - integrator.uprev
    update_nordsieck_vector!(cache)
  else
    # Reset time
    for i in endof(tau):-1:2
      tau[i] = tau[i-1]
    end
    tau[1] = dt
    # Rescale
    dt != tau[2] && nordsieck_rescale!(cache)
    @. integrator.k[1] = z[2]/dt
    # Perform 5th order Adams method in Nordsieck form
    perform_predict!(cache)
    calc_coeff!(cache)
    isucceed = nlsolve_functional!(integrator, cache)
    if !isucceed
      integrator.force_stepfail = true
      # rewind Nordsieck vector
      nordsieck_rewind!(cache)
      return nothing
    end

    ################################### Error estimation

    if integrator.opts.adaptive
      calculate_residuals!(atmp, const_cache.Δ, uprev, integrator.u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
      integrator.EEst = integrator.opts.internalnorm(atmp) * const_cache.c_LTE
    end

    # Correct Nordsieck vector
    const_cache.step = min(const_cache.step+1, 12)
    update_nordsieck_vector!(cache)

    ################################### Finalize

    @. integrator.k[2] = const_cache.z[2]/dt
  end
  const_cache.n_wait -= 1
  return nothing
end
