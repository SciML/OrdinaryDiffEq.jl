function initialize!(integrator,cache::AN5ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::AN5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack z,l,m,tq,tau,tsit5tab = cache
  z[1] = uprev
  # handle callbacks, rewind back to order one.
  if integrator.u_modified
    cache.step = 1
  end
  # Nordsieck form needs to build the history vector, so one has to start from order one.
  if cache.step <= 4
    if cache.step == 1
      # Grow the Nordsieck vector to length 2+1
      perform_step!(integrator, tsit5tab, repeat_step)
      # The first history vector is `h ydot`
      cache.z[2] = integrator.fsalfirst*dt
      cache.step += 1
    elseif cache.step == 2
      # Grow the Nordsieck vector to length 3+1
      cache.step += 1
    elseif cache.step == 3
      # Grow the Nordsieck vector to length 4+1
      cache.step += 1
    elseif cache.step == 4
      # Grow the Nordsieck vector to length 5+1
      cache.step += 1
    end
  else
  # Perform 5th order Adams method in Nordsieck form

  end

  ################################### Error estimation

  if integrator.opts.adaptive
    utilde =
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  ################################### Finalize

  # Reset time
  for i in endof(tau):-1:2
    tau[i] = tau[i-1]
  end
  tau[1] = dt

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = uₙ
  return nothing
end

# Work on the out-of-place version first
#=
function initialize!(integrator, cache::AN5Cache)
  integrator.kshortsize = 2
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # For the interpolation, needs k at the updated point
end

@muladd function perform_step!(integrator, cache::AN5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack hist1,hist2,hist3,hist4,hist5,l,m,tq,tau = cache
  # handle callbacks, rewind back to order one.
  if integrator.u_modified
    cache.step = 1
  end
  # Nordsieck form needs to build the history vector, so one has to start from order one.
  if cache.step <= 4
    if cache.step == 1
      # Grow the Nordsieck vector to length 2+1
      perform_step!(integrator, EulerConstantCache(), repeat_step)
      # The first history vector is `h ydot`
      cache.hist1 = integrator.fsalfirst*dt
      cache.step += 1
      return nothing
    elseif cache.step == 2
      # Grow the Nordsieck vector to length 3+1
      cache.step += 1
      return nothing
    elseif cache.step == 3
      # Grow the Nordsieck vector to length 4+1
      cache.step += 1
      return nothing
    elseif cache.step == 4
      # Grow the Nordsieck vector to length 5+1
      cache.step += 1
      return nothing
    end
  end

  # WIP

  ################################### Finalize

  # Reset time
  for i in endof(tau):-1:2
    tau[i] = tau[i-1]
  end
  tau[1] = dt

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = uₙ
  return nothing
end
=#
