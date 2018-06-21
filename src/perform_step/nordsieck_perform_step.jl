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
    cache.order = 1
  end
  # Nordsieck form needs to build the history vector
  if cache.order == 1
    # Start the Nordsieck vector in one shot!
    perform_step!(integrator, tsit5tab, repeat_step)
    cache.order = 4
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
      cache.order = 1
    end
  else
    # Reset time
    tmp = tau[6]
    for i in 5:-1:1
      tau[i+1] = tau[i]
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
      if integrator.EEst > one(integrator.EEst)
        for i in 1:5
          tau[i] = tau[i+1]
        end
        tau[6] = tmp
      end
    end

    # Correct Nordsieck vector
    cache.order = 5
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
    const_cache.order = 1
  end
  # Nordsieck form needs to build the history vector
  if const_cache.order == 1
    ## Start the Nordsieck vector in two shots!
    perform_step!(integrator, tsit5cache, repeat_step)
    copy!(tmp, integrator.u)
    const_cache.order = 4
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
      const_cache.order = 1
    end
  else
    # Reset time
    tmp = tau[6]
    for i in 5:-1:1
      tau[i+1] = tau[i]
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
      if integrator.EEst > one(integrator.EEst)
        for i in 1:5
          tau[i] = tau[i+1]
        end
        tau[6] = tmp
      end
    end

    # Correct Nordsieck vector
    const_cache.order = 5
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
  if integrator.u_modified || integrator.iter == 1
    cache.order = 1
    z[1] = integrator.uprev
    z[2] = f(uprev, p, t)*dt
    tau[1] = dt
  end
  # Reset time
  tmp = tau[13]
  for i in 12:-1:1
    tau[i+1] = tau[i]
  end
  tau[1] = dt
  dt != tau[2] && nordsieck_rescale!(cache)
  integrator.k[1] = z[2]/dt
  # Perform 5th order Adams method in Nordsieck form
  perform_predict!(cache)
  cache.order = min(cache.nextorder, 12)
  calc_coeff!(cache)
  isucceed = nlsolve_functional!(integrator, cache)
  if !isucceed
    # rewind Nordsieck vector
    integrator.force_stepfail = true
    nordsieck_rewind!(cache)
    return nothing
  end

  # Correct Nordsieck vector
  update_nordsieck_vector!(cache)

  ################################### Finalize
  cache.n_wait -= 1
  if nordsieck_change_order(cache, 1) && cache.order != 12
    cache.z[end] = cache.Δ
    cache.prev_𝒟 = cache.c_𝒟
  end

  integrator.k[2] = cache.z[2]/dt
  ################################### Error estimation
  if integrator.opts.adaptive
    atmp = calculate_residuals(cache.Δ, uprev, integrator.u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp) * cache.c_LTE
    if integrator.EEst > one(integrator.EEst)
      for i in 1:12
        tau[i] = tau[i+1]
      end
      tau[13] = tmp
    end
  end
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
  if integrator.u_modified || integrator.iter == 1
    const_cache.order = 1
    @. z[1] = integrator.uprev
    f(z[2], uprev, p, t)
    @. z[2] = z[2]*dt
    tau[1] = dt
  end
  tmp = tau[13]
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
  const_cache.order = min(const_cache.nextorder, 12)
  calc_coeff!(cache)
  isucceed = nlsolve_functional!(integrator, cache)
  if !isucceed
    integrator.force_stepfail = true
    # rewind Nordsieck vector
    nordsieck_rewind!(cache)
    return nothing
  end

  # Correct Nordsieck vector
  update_nordsieck_vector!(cache)
  const_cache.n_wait -= 1

  ################################### Finalize
  if nordsieck_change_order(cache, 1) && const_cache.order != 12
    const_cache.z[end] = const_cache.Δ
    const_cache.prev_𝒟 = const_cache.c_𝒟
  end

  ################################### Error estimation

  if integrator.opts.adaptive
    calculate_residuals!(atmp, const_cache.Δ, uprev, integrator.u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp) * const_cache.c_LTE
    if integrator.EEst > one(integrator.EEst)
      for i in 1:12
        tau[i] = tau[i+1]
      end
      tau[13] = tmp
    end
  end
  return nothing
end
