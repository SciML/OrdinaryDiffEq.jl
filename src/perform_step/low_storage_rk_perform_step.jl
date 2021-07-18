
# 2N low storage methods
function initialize!(integrator,cache::LowStorageRK2NConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
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
    integrator.destats.nf += 1
    tmp = A2end[i]*tmp + dt*k
    u   = u + B2end[i]*tmp
  end

  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  integrator.fsalfirst = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.u = u
end

function initialize!(integrator,cache::LowStorageRK2NCache)
  @unpack k, tmp, williamson_condition = cache
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = k
  integrator.fsalfirst = k # used for get_du
  integrator.fsallast = k
  integrator.f(k ,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::LowStorageRK2NCache,repeat_step=false)
  @unpack t,dt,u,f,p = integrator
  @unpack k,tmp,williamson_condition,stage_limiter!,step_limiter! = cache
  @unpack A2end,B1,B2end,c2end = cache.tab

  # u1
  f(k, u, p, t)
  integrator.destats.nf += 1
  @.. tmp = dt*k
  @.. u   = u + B1*tmp
  # other stages
  for i in eachindex(A2end)
    if williamson_condition
      f(ArrayFuse(tmp, u, (A2end[i], dt, B2end[i])), u, p, t+c2end[i]*dt)
    else
      @.. tmp = A2end[i]*tmp
      stage_limiter!(u, integrator, p, t+c2end[i]*dt)
      f(k, u, p, t+c2end[i]*dt)
      @.. tmp += dt * k
      @.. u   = u + B2end[i]*tmp
    end
    integrator.destats.nf += 1
  end
  step_limiter!(u, integrator, p, t+dt)
end

# 2C low storage methods
function initialize!(integrator,cache::LowStorageRK2CConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
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
  k = integrator.fsalfirst = f(u, p, t)
  integrator.k[1] = integrator.fsalfirst
  integrator.destats.nf += 1
  u = u + B1*dt*k

  # other stages
  for i in eachindex(A2end)
    tmp = u + A2end[i]*dt*k
    k   = f(tmp, p, t+c2end[i]*dt)
    integrator.destats.nf += 1
    u   = u + B2end[i]*dt*k
  end

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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::LowStorageRK2CCache,repeat_step=false)
  @unpack t,dt,u,f,p = integrator
  @unpack k,fsalfirst,tmp = cache
  @unpack A2end,B1,B2end,c2end = cache.tab

  # u1
  @.. k = integrator.fsalfirst
  @.. u = u + B1*dt*k

  # other stages
  for i in eachindex(A2end)
    @.. tmp = u + A2end[i]*dt*k
    f(k, tmp, p, t+c2end[i]*dt)
    integrator.destats.nf += 1
    @.. u   = u + B2end[i]*dt*k
  end

  f(k, u, p, t+dt)
  integrator.destats.nf += 1
end



# 3S low storage methods
function initialize!(integrator,cache::LowStorageRK3SConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
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
    integrator.destats.nf += 1
    tmp = tmp + δ2end[i]*u
    u   = γ12end[i]*u + γ22end[i]*tmp + γ32end[i]*uprev + β2end[i]*dt*k
  end

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.destats.nf += 1
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::LowStorageRK3SCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,tmp = cache
  @unpack γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end = cache.tab

  # u1
  @.. tmp = u
  @.. u   = tmp + β1*dt*integrator.fsalfirst

  # other stages
  for i in eachindex(γ12end)
    f(k, u, p, t+c2end[i]*dt)
    integrator.destats.nf += 1
    @.. tmp = tmp + δ2end[i]*u
    @.. u   = γ12end[i]*u + γ22end[i]*tmp + γ32end[i]*uprev + β2end[i]*dt*k
  end

  f(k, u, p, t+dt)
  integrator.destats.nf += 1
end



# 3S+ low storage methods: 3S methods adding another memory location for the embedded method (non-FSAL version)
function initialize!(integrator,cache::LowStorageRK3SpConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::LowStorageRK3SpConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end = cache

  # u1
  integrator.fsalfirst = f(uprev, p, t)
  integrator.destats.nf += 1
  integrator.k[1] = integrator.fsalfirst
  tmp = uprev
  u   = tmp + β1*dt*integrator.fsalfirst
  if integrator.opts.adaptive
    utilde = bhat1*dt*integrator.fsalfirst
  end

  # other stages
  for i in eachindex(γ12end)
    k   = f(u, p, t+c2end[i]*dt)
    integrator.destats.nf += 1
    tmp = tmp + δ2end[i]*u
    u   = γ12end[i]*u + γ22end[i]*tmp + γ32end[i]*uprev + β2end[i]*dt*k
    if integrator.opts.adaptive
      utilde = utilde + bhat2end[i]*dt*k
    end
  end

  if integrator.opts.adaptive
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  integrator.u = u
end

function initialize!(integrator,cache::LowStorageRK3SpCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::LowStorageRK3SpCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,utilde,atmp = cache
  @unpack γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end = cache.tab

  # u1
  f(integrator.fsalfirst, uprev, p, t)
  integrator.destats.nf += 1
  @.. tmp = uprev
  @.. u   = tmp + β1*dt*integrator.fsalfirst
  if integrator.opts.adaptive
    @.. utilde = bhat1*dt*integrator.fsalfirst
  end

  # other stages
  for i in eachindex(γ12end)
    f(k, u, p, t+c2end[i]*dt)
    integrator.destats.nf += 1
    @.. tmp = tmp + δ2end[i]*u
    @.. u   = γ12end[i]*u + γ22end[i]*tmp + γ32end[i]*uprev + β2end[i]*dt*k
    if integrator.opts.adaptive
      @.. utilde = utilde + bhat2end[i]*dt*k
    end
  end

  if integrator.opts.adaptive
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
end



# 3S+ FSAL low storage methods: 3S methods adding another memory location for the embedded method (FSAL version)
function initialize!(integrator,cache::LowStorageRK3SpFSALConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::LowStorageRK3SpFSALConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end, bhatfsal = cache

  # u1
  tmp = uprev
  u   = tmp + β1*dt*integrator.fsalfirst
  if integrator.opts.adaptive
    utilde = bhat1*dt*integrator.fsalfirst
  end

  # other stages
  for i in eachindex(γ12end)
    k   = f(u, p, t+c2end[i]*dt)
    integrator.destats.nf += 1
    tmp = tmp + δ2end[i]*u
    u   = γ12end[i]*u + γ22end[i]*tmp + γ32end[i]*uprev + β2end[i]*dt*k
    if integrator.opts.adaptive
      utilde = utilde + bhat2end[i]*dt*k
    end
  end

  # FSAL
  integrator.fsallast = f(u, p, t+dt)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    utilde = utilde + bhatfsal*dt*integrator.fsallast
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::LowStorageRK3SpFSALCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::LowStorageRK3SpFSALCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,tmp,utilde,atmp = cache
  @unpack γ12end, γ22end, γ32end, δ2end, β1, β2end, c2end, bhat1, bhat2end, bhatfsal = cache.tab

  # u1
  @.. tmp = uprev
  @.. u   = tmp + β1*dt*integrator.fsalfirst
  if integrator.opts.adaptive
    @.. utilde = bhat1*dt*integrator.fsalfirst
  end

  # other stages
  for i in eachindex(γ12end)
    f(k, u, p, t+c2end[i]*dt)
    integrator.destats.nf += 1
    @.. tmp = tmp + δ2end[i]*u
    @.. u   = γ12end[i]*u + γ22end[i]*tmp + γ32end[i]*uprev + β2end[i]*dt*k
    if integrator.opts.adaptive
      @.. utilde = utilde + bhat2end[i]*dt*k
    end
  end

  # FSAL
  f(k, u, p, t+dt)
  integrator.destats.nf += 1

  if integrator.opts.adaptive
    @.. utilde = utilde + bhatfsal*dt*k
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
end



# 2R+ low storage methods
function initialize!(integrator,cache::LowStorageRK2RPConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::LowStorageRK2RPConstantCache,repeat_step=false)
  @unpack t,dt,u,uprev,f,fsalfirst,p = integrator
  @unpack Aᵢ,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ = cache

  k   = fsalfirst
  integrator.opts.adaptive && (tmp = zero(uprev))

  #stages 1 to s-1
  for i in eachindex(Aᵢ)
    integrator.opts.adaptive && (tmp = tmp + (Bᵢ[i] - B̂ᵢ[i])*dt*k)
    gprev = u + Aᵢ[i]*dt*k
    u = u + Bᵢ[i]*dt*k
    k = f(gprev, p, t + Cᵢ[i]*dt)
    integrator.destats.nf += 1
  end

  #last stage
  integrator.opts.adaptive && (tmp = tmp + (Bₗ - B̂ₗ)*dt*k)
  u   = u  + Bₗ*dt*k

  #Error estimate
  if integrator.opts.adaptive
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.destats.nf += 1
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::LowStorageRK2RPCache,repeat_step=false)
  @unpack t,dt,u,uprev,f,fsalfirst,p = integrator
  @unpack k,gprev,tmp,atmp = cache
  @unpack Aᵢ,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ = cache.tab

  @.. k   = fsalfirst
  integrator.opts.adaptive && (@.. tmp = zero(uprev))

  #stages 1 to s-1
  for i in eachindex(Aᵢ)
    integrator.opts.adaptive && (@.. tmp = tmp + (Bᵢ[i] - B̂ᵢ[i])*dt*k)
    @.. gprev = u + Aᵢ[i]*dt*k
    @.. u     = u + Bᵢ[i]*dt*k
    f(k, gprev, p, t + Cᵢ[i]*dt)
    integrator.destats.nf += 1
  end

  #last stage
  integrator.opts.adaptive && (@.. tmp = tmp + (Bₗ - B̂ₗ)*dt*k)
  @.. u   = u  + Bₗ*dt*k

  #Error estimate
  if integrator.opts.adaptive
    calculate_residuals!(atmp,tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  f(k, u, p, t+dt)
  integrator.destats.nf += 1
end



# 3R+ low storage methods
function initialize!(integrator,cache::LowStorageRK3RPConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::LowStorageRK3RPConstantCache,repeat_step=false)
  @unpack t,dt,u,uprev,f,fsalfirst,p = integrator
  @unpack Aᵢ₁,Aᵢ₂,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ = cache

  fᵢ₋₂  = zero(fsalfirst)
  k     = fsalfirst
  uᵢ₋₁  = uprev
  uᵢ₋₂  = uprev
  integrator.opts.adaptive && (tmp = zero(uprev))

  #stages 1 to s-1
  for i in eachindex(Aᵢ₁)
    integrator.opts.adaptive && (tmp = tmp + (Bᵢ[i] - B̂ᵢ[i])*dt*k)
    gprev = uᵢ₋₂ + (Aᵢ₁[i]*k+Aᵢ₂[i]*fᵢ₋₂)*dt
    u     = u + Bᵢ[i]*dt*k
    fᵢ₋₂  = k
    uᵢ₋₂  = uᵢ₋₁
    uᵢ₋₁  = u
    k     = f(gprev, p, t + Cᵢ[i]*dt)
    integrator.destats.nf += 1
  end

  #last stage
  integrator.opts.adaptive && (tmp = tmp + (Bₗ - B̂ₗ)*dt*k)
  u   = u  + Bₗ*dt*k

  #Error estimate
  if integrator.opts.adaptive
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.destats.nf += 1
  integrator.u = u
end

function initialize!(integrator,cache::LowStorageRK3RPCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t)
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::LowStorageRK3RPCache,repeat_step=false)
  @unpack t,dt,u,uprev,f,fsalfirst,p = integrator
  @unpack k,uᵢ₋₁,uᵢ₋₂,gprev,fᵢ₋₂,tmp,atmp = cache
  @unpack Aᵢ₁,Aᵢ₂,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ = cache.tab

  @.. fᵢ₋₂  = zero(fsalfirst)
  @.. k     = fsalfirst
  integrator.opts.adaptive && (@.. tmp = zero(uprev))
  @.. uᵢ₋₁  = uprev
  @.. uᵢ₋₂  = uprev

  #stages 1 to s-1
  for i in eachindex(Aᵢ₁)
    integrator.opts.adaptive && (@.. tmp = tmp + (Bᵢ[i] - B̂ᵢ[i])*dt*k)
    @.. gprev = uᵢ₋₂ + (Aᵢ₁[i]*k+Aᵢ₂[i]*fᵢ₋₂)*dt
    @.. u     = u + Bᵢ[i]*dt*k
    @.. fᵢ₋₂  = k
    @.. uᵢ₋₂  = uᵢ₋₁
    @.. uᵢ₋₁  = u
    f(k, gprev, p, t + Cᵢ[i]*dt)
    integrator.destats.nf += 1
  end

  #last stage
  integrator.opts.adaptive && (@.. tmp = tmp + (Bₗ - B̂ₗ)*dt*k)
  @.. u   = u  + Bₗ*dt*k

  #Error estimate
  if integrator.opts.adaptive
    calculate_residuals!(atmp,tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  f(k, u, p, t+dt)
  integrator.destats.nf += 1
end


# 4R+ low storage methods
function initialize!(integrator,cache::LowStorageRK4RPConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::LowStorageRK4RPConstantCache,repeat_step=false)
  @unpack t,dt,u,uprev,f,fsalfirst,p = integrator
  @unpack Aᵢ₁,Aᵢ₂,Aᵢ₃,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ = cache

  fᵢ₋₂  = zero(fsalfirst)
  fᵢ₋₃  = zero(fsalfirst)
  k     = fsalfirst
  uᵢ₋₁  = uprev
  uᵢ₋₂  = uprev
  uᵢ₋₃  = uprev
  integrator.opts.adaptive && (tmp = zero(uprev))

  #stages 1 to s-1
  for i in eachindex(Aᵢ₁)
    integrator.opts.adaptive && (tmp = tmp + (Bᵢ[i] - B̂ᵢ[i])*dt*k)
    gprev = uᵢ₋₃ + (Aᵢ₁[i]*k+Aᵢ₂[i]*fᵢ₋₂+Aᵢ₃[i]*fᵢ₋₃)*dt
    u     = u + Bᵢ[i]*dt*k
    fᵢ₋₃  = fᵢ₋₂
    fᵢ₋₂  = k
    uᵢ₋₃  = uᵢ₋₂
    uᵢ₋₂  = uᵢ₋₁
    uᵢ₋₁  = u
    k     = f(gprev, p, t + Cᵢ[i]*dt)
    integrator.destats.nf += 1
  end

  #last stage
  integrator.opts.adaptive && (tmp = tmp + (Bₗ - B̂ₗ)*dt*k)
  u   = u  + Bₗ*dt*k

  #Error estimate
  if integrator.opts.adaptive
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.destats.nf += 1
  integrator.u = u
end

function initialize!(integrator,cache::LowStorageRK4RPCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t)
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::LowStorageRK4RPCache,repeat_step=false)
  @unpack t,dt,u,uprev,f,fsalfirst,p = integrator
  @unpack k,uᵢ₋₁,uᵢ₋₂,uᵢ₋₃,gprev,fᵢ₋₂,fᵢ₋₃,tmp,atmp = cache
  @unpack Aᵢ₁,Aᵢ₂,Aᵢ₃,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ = cache.tab

  @.. fᵢ₋₂  = zero(fsalfirst)
  @.. fᵢ₋₃  = zero(fsalfirst)
  @.. k     = fsalfirst
  integrator.opts.adaptive && (@.. tmp = zero(uprev))
  @.. uᵢ₋₁  = uprev
  @.. uᵢ₋₂  = uprev
  @.. uᵢ₋₃  = uprev


  #stages 1 to s-1
  for i in eachindex(Aᵢ₁)
    integrator.opts.adaptive && (@.. tmp = tmp + (Bᵢ[i] - B̂ᵢ[i])*dt*k)
    @.. gprev = uᵢ₋₃ + (Aᵢ₁[i]*k+Aᵢ₂[i]*fᵢ₋₂+Aᵢ₃[i]*fᵢ₋₃)*dt
    @.. u     = u + Bᵢ[i]*dt*k
    @.. fᵢ₋₃  = fᵢ₋₂
    @.. fᵢ₋₂  = k
    @.. uᵢ₋₃  = uᵢ₋₂
    @.. uᵢ₋₂  = uᵢ₋₁
    @.. uᵢ₋₁  = u
    f(k, gprev, p, t + Cᵢ[i]*dt)
    integrator.destats.nf += 1
  end

  #last stage
  integrator.opts.adaptive && (@.. tmp = tmp + (Bₗ - B̂ₗ)*dt*k)
  @.. u   = u  + Bₗ*dt*k

  #Error estimate
  if integrator.opts.adaptive
    calculate_residuals!(atmp,tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  f(k, u, p, t+dt)
  integrator.destats.nf += 1
end


# 5R+ low storage methods
function initialize!(integrator,cache::LowStorageRK5RPConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::LowStorageRK5RPConstantCache,repeat_step=false)
  @unpack t,dt,u,uprev,f,fsalfirst,p = integrator
  @unpack Aᵢ₁,Aᵢ₂,Aᵢ₃,Aᵢ₄,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ = cache

  fᵢ₋₂  = zero(fsalfirst)
  fᵢ₋₃  = zero(fsalfirst)
  fᵢ₋₄  = zero(fsalfirst)
  k     = fsalfirst
  uᵢ₋₁  = uprev
  uᵢ₋₂  = uprev
  uᵢ₋₃  = uprev
  uᵢ₋₄  = uprev
  integrator.opts.adaptive && (tmp = zero(uprev))

  #stages 1 to s-1
  for i in eachindex(Aᵢ₁)
    integrator.opts.adaptive && (tmp = tmp + (Bᵢ[i] - B̂ᵢ[i])*dt*k)
    gprev = uᵢ₋₄ + (Aᵢ₁[i]*k+Aᵢ₂[i]*fᵢ₋₂+Aᵢ₃[i]*fᵢ₋₃+Aᵢ₄[i]*fᵢ₋₄)*dt
    u     = u + Bᵢ[i]*dt*k
    fᵢ₋₄  = fᵢ₋₃
    fᵢ₋₃  = fᵢ₋₂
    fᵢ₋₂  = k
    uᵢ₋₄  = uᵢ₋₃
    uᵢ₋₃  = uᵢ₋₂
    uᵢ₋₂  = uᵢ₋₁
    uᵢ₋₁  = u
    k     = f(gprev, p, t + Cᵢ[i]*dt)
    integrator.destats.nf += 1
  end

  #last stage
  integrator.opts.adaptive && (tmp = tmp + (Bₗ - B̂ₗ)*dt*k)
  u   = u  + Bₗ*dt*k

  #Error estimate
  if integrator.opts.adaptive
    atmp = calculate_residuals(tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  integrator.k[1] = integrator.fsalfirst
  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.destats.nf += 1
  integrator.u = u
end

function initialize!(integrator,cache::LowStorageRK5RPCache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t)
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::LowStorageRK5RPCache,repeat_step=false)
  @unpack t,dt,u,uprev,f,fsalfirst,p = integrator
  @unpack k,uᵢ₋₁,uᵢ₋₂,uᵢ₋₃,uᵢ₋₄,gprev,fᵢ₋₂,fᵢ₋₃,fᵢ₋₄,tmp,atmp = cache
  @unpack Aᵢ₁,Aᵢ₂,Aᵢ₃,Aᵢ₄,Bₗ,B̂ₗ,Bᵢ,B̂ᵢ,Cᵢ = cache.tab

  @.. fᵢ₋₂  = zero(fsalfirst)
  @.. fᵢ₋₃  = zero(fsalfirst)
  @.. fᵢ₋₄  = zero(fsalfirst)
  @.. k     = fsalfirst
  integrator.opts.adaptive && (@.. tmp = zero(uprev))
  @.. uᵢ₋₁  = uprev
  @.. uᵢ₋₂  = uprev
  @.. uᵢ₋₃  = uprev
  @.. uᵢ₋₄  = uprev


  #stages 1 to s-1
  for i in eachindex(Aᵢ₁)
    integrator.opts.adaptive && (@.. tmp = tmp + (Bᵢ[i] - B̂ᵢ[i])*dt*k)
    @.. gprev = uᵢ₋₄ + (Aᵢ₁[i]*k+Aᵢ₂[i]*fᵢ₋₂+Aᵢ₃[i]*fᵢ₋₃+Aᵢ₄[i]*fᵢ₋₄)*dt
    @.. u     = u + Bᵢ[i]*dt*k
    @.. fᵢ₋₄  = fᵢ₋₃
    @.. fᵢ₋₃  = fᵢ₋₂
    @.. fᵢ₋₂  = k
    @.. uᵢ₋₄  = uᵢ₋₃
    @.. uᵢ₋₃  = uᵢ₋₂
    @.. uᵢ₋₂  = uᵢ₋₁
    @.. uᵢ₋₁  = u
    f(k, gprev, p, t + Cᵢ[i]*dt)
    integrator.destats.nf += 1
  end

  #last stage
  integrator.opts.adaptive && (@.. tmp = tmp + (Bₗ - B̂ₗ)*dt*k)
  @.. u   = u  + Bₗ*dt*k

  #Error estimate
  if integrator.opts.adaptive
    calculate_residuals!(atmp,tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end

  f(k, u, p, t+dt)
  integrator.destats.nf += 1
end
