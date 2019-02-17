
# 2N low storage methods
function initialize!(integrator,cache::LowStorageRK2NConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::LowStorageRK2NConstantCache,repeat_step=false)
  @unpack t,dt,u,f,p = integrator
  @unpack A2end,B1,B2end,c2end = cache

  # u1
  tmp = dt*integrator.fsalfirst
  u   = u + B1*tmp

  # other stages
  for i in eachindex(A2end)
    k = f(u, p, t+c2end[i]*dt)
    tmp = A2end[i]*tmp + dt*k
    u   = u + B2end[i]*tmp
  end

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::LowStorageRK2NCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::LowStorageRK2NCache,repeat_step=false)
  @unpack t,dt,u,f,p = integrator
  @unpack k,fsalfirst,tmp = cache
  @unpack A2end,B1,B2end,c2end = cache.tab

  # u1
  @. tmp = dt*fsalfirst
  @. u   = u + B1*tmp

  # other stages
  for i in eachindex(A2end)
    f(k, u, p, t+c2end[i]*dt)
    @. tmp = A2end[i]*tmp + dt*k
    @. u   = u + B2end[i]*tmp
  end

  f(k, u, p, t+dt)
end



# 2C low storage methods
function initialize!(integrator,cache::LowStorageRK2CConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::LowStorageRK2CConstantCache,repeat_step=false)
  @unpack t,dt,u,f,p = integrator
  @unpack A2end,B1,B2end,c2end = cache

  # u1
  k = integrator.fsalfirst
  u = u + B1*dt*k

  # other stages
  for i in eachindex(A2end)
    tmp = u + A2end[i]*dt*k
    k   = f(tmp, p, t+c2end[i]*dt)
    u   = u + B2end[i]*dt*k
  end

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::LowStorageRK2CCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::LowStorageRK2CCache,repeat_step=false)
  @unpack t,dt,u,f,p = integrator
  @unpack k,fsalfirst,tmp = cache
  @unpack A2end,B1,B2end,c2end = cache.tab

  # u1
  @. k = integrator.fsalfirst
  @. u = u + B1*dt*k

  # other stages
  for i in eachindex(A2end)
    @. tmp = u + A2end[i]*dt*k
    f(k, tmp, p, t+c2end[i]*dt)
    @. u   = u + B2end[i]*dt*k
  end

  f(k, u, p, t+dt)
end



# 3S low storage methods
function initialize!(integrator,cache::LowStorageRK3SConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::LowStorageRK3SConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end = cache

  # u1
  tmp = u
  u   = tmp + β1*dt*integrator.fsalfirst

  # other stages
  for i in eachindex(γ12end)
    k   = f(u, p, t+c2end[i]*dt)
    tmp = tmp + δ2end[i]*u
    u   = γ12end[i]*u + γ22end[i]*tmp + γ32end[i]*uprev + β2end[i]*dt*k
  end

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::LowStorageRK3SCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::LowStorageRK3SCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,tmp = cache
  @unpack γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end = cache.tab

  # u1
  @. tmp = u
  @. u   = tmp + β1*dt*integrator.fsalfirst

  # other stages
  for i in eachindex(γ12end)
    f(k, u, p, t+c2end[i]*dt)
    @. tmp = tmp + δ2end[i]*u
    @. u   = γ12end[i]*u + γ22end[i]*tmp + γ32end[i]*uprev + β2end[i]*dt*k
  end

  f(k, u, p, t+dt)
end



# 2R+ low storage methods
function initialize!(integrator,cache::LowStorageRK2RPConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::LowStorageRK2RPConstantCache,repeat_step=false)
  @unpack t,dt,u,uprev,f,fsalfirst,p = integrator
  @unpack A1nm1,Bend,Bhend,B1nm1,Bh1nm1,C1nm1 = cache

  k   = fsalfirst
  integrator.opts.adaptive && (tmp = zero(uprev))

  #stages 1 to s-1
  for i in eachindex(A1nm1)
    integrator.opts.adaptive && (tmp = tmp + (B1nm1[i] - Bh1nm1[i])*dt*k)
    u = u + B1nm1[i]*dt*k
    gprev = u + (A1nm1[i] - B1nm1[i])*dt*k
    k = f(u + (A1nm1[i] - B1nm1[i])*dt*k, p, t + C1nm1[i]*dt)
  end

  #last stage
  integrator.opts.adaptive && (tmp = tmp + (Bend - Bhend)*dt*k)
  u   = u  + Bend*dt*k

  #Error estimate
  if integrator.opts.adaptive
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.u = u
end

function initialize!(integrator,cache::LowStorageRK2RPCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t)
end

@muladd function perform_step!(integrator,cache::LowStorageRK2RPCache,repeat_step=false)
  @unpack t,dt,u,uprev,f,fsalfirst,p = integrator
  @unpack k,gprev,tmp,atmp = cache
  @unpack A1nm1,Bend,Bhend,B1nm1,Bh1nm1,C1nm1 = cache.tab

  @. k   = fsalfirst
  integrator.opts.adaptive && (@. tmp = zero(uprev))

  #stages 1 to s-1
  for i in eachindex(A1nm1)
    integrator.opts.adaptive && (@. tmp = tmp + (B1nm1[i] - Bh1nm1[i])*dt*k)
    @. u     = u + B1nm1[i]*dt*k
    @. gprev = u + (A1nm1[i] - B1nm1[i])*dt*k
    f(k, gprev, p, t + C1nm1[i]*dt)
  end

  #last stage
  integrator.opts.adaptive && (@. tmp = tmp + (Bend - Bhend)*dt*k)
  @. u   = u  + Bend*dt*k

  #Error estimate
  if integrator.opts.adaptive
    calculate_residuals!(atmp,tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  f(k, u, p, t+dt)
end
